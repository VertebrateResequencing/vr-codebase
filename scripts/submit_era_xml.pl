#!/usr/bin/env perl
use strict;
use warnings;
no warnings 'uninitialized';
use XML::Simple;
use Data::Dumper;
use VertRes::Utils::VRTrackFactory;
use Getopt::Long;
 
my ($file_prefix, $help, $spp,$verbose,$skip_test, $skip_live);

GetOptions(
    'p|prefix=s'  =>  \$file_prefix,
    's|spp=s'     =>  \$spp,
    'v|verbose'   =>  \$verbose,
    'n|no_test'   =>  \$skip_test,
    'l|no_live'   =>  \$skip_live,
    'h|help'	  =>  \$help,
    );


my %db_for_spp = ( 'mouse'  => 'mouse_reseq_track',
		    'g1k'   => 'g1k_track',
                 );

# For testing
#my %db_for_spp = ( 'mouse'  => 'mouse_reseq_track_test',
#		    'g1k'   => 'g1k_track_test',
#	  );

my $db = $db_for_spp{$spp};

my $subfile = "$file_prefix.submission.xml";
my $runfile = "$file_prefix.run.xml";
my $expfile = "$file_prefix.experiment.xml";

(-f $subfile && -f $runfile && -f $expfile && $db && !$help) or die <<USAGE;
    Usage: $0   
                --prefix    <file prefix to submission files>
                --spp       <species, i.e. g1k or mouse>
                [--no_test  <skip test submission, just submit live>
                [--no_live  <skip live submission, just test>
                [--verbose  <verbose output>]
                [--help     <this message>]

e.g.
$0 --prefix \$G1K/SRA/20100607/g1k-srp000546-sc-20100607 --spp g1k

Submits submission, experiment and run XMLs to ERA and updates tracking
database with the accessions returned.

First submits to the test instance to check if the submission is valid, and
only if so does it submit to the live instance and update the database.  This 
behaviour can be disabled by --no_test; you probably only want to do this if
you have already sent the submission to the test server and now it will fail
due to the submission existing already.

USAGE

print "Species: $spp\n";
print "Database: $db\n";
my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $db,
                                                          mode => 'rw');
unless ($vrtrack){
    die "Can't connect to tracking db $db\n";
}

# environment for curl
$ENV{HTTP_PROXY} = "http://webcache.sanger.ac.uk:3128";
$ENV{HTTPS_PROXY} = "http://webcache.sanger.ac.uk:3128";
my $test_url = 'https://wwwdev.ebi.ac.uk/ena/submit/drop-box/submit?auth=ERA+era-drop-2+zVFfPuFNDgZIfIhnNI6tx95Lp8o%3D';
my $live_url = 'https://www.ebi.ac.uk/ena/submit/drop-box/submit?auth=ERA+era-drop-2+zVFfPuFNDgZIfIhnNI6tx95Lp8o%3D';

my $xmlhash;

# Make test submission
if ($skip_test){
    print "Skipping test submission.  Making real submission\n";
}
else {
    $xmlhash = make_submission($subfile, $expfile, $runfile, $test_url);

    if (submission_was_successful($xmlhash)){
        print "Test submission successful.  Making real submission\n";
    }
    else {
        print "Test submission failed\n";
        print join "\n", @{$xmlhash->{MESSAGES}->[0]->{ERROR}},"\n";
        exit 1;
    }
}

# Make live submission
if ($skip_live){
    print "Skipping real submission.\n";
    exit;
}
else {
    $xmlhash = make_submission($subfile, $expfile, $runfile, $live_url, "$file_prefix.receipt.xml");

    unless (submission_was_successful($xmlhash)){
        print "ERA submission failed\n";
        print join "\n", @{$xmlhash->{MESSAGES}->[0]->{ERROR}},"\n";
        exit 1;
    }
}

my($submission_alias, $submission_acc) = get_submission_alias_acc_from_xml($xmlhash);
my $lanehash = get_lanes_from_xml($xmlhash);

# get all the lanes for this submission alias
# Note this is both accessioned and unaccessioned. We'll check both.
my $sql = qq[select l.hierarchy_name from submission sub, latest_lane l where l.submission_id = sub.submission_id and sub.name="$submission_alias"];
my %vr_lanes;
my $sth = $vrtrack->{_dbh}->prepare($sql);

if ($sth->execute()){
    foreach(@{$sth->fetchall_arrayref()}){
        $vr_lanes{$_->[0]} = "";
    }
}
else {
    die(sprintf('Cannot retrieve vr_lanes: %s', $DBI::errstr));
}

# check two sets of lanes are congruent
foreach my $era_lane (keys %$lanehash){
    unless (exists $vr_lanes{$era_lane}){
        die "$era_lane is not in vr_track for submission $submission_alias\n";
    }
}

foreach my $vr_lane (keys %vr_lanes){
    unless (exists $lanehash->{$vr_lane}){
        die "$vr_lane is not in ERA for submission $submission_alias\n";
    }
}
        
# right, now update the accessions
foreach my $era_lane (keys %$lanehash){
    my $lane_acc = $lanehash->{$era_lane};
    my $vr_lane = VRTrack::Lane->new_by_hierarchy_name($vrtrack,$era_lane);
    unless ($vr_lane){
       print  "Can't find lane $era_lane in database\n";
       next;
    }
    if ($vr_lane->submission){
        if ($vr_lane->submission->name ne $submission_alias){
            # major problem somewhere
            die "Existing submission for $era_lane is ",$vr_lane->submission->name, " not $submission_alias";
        }
    }
    else {
        warn "$era_lane has no associated submission.  Attempting to add.\n";
        # should never happen, as we just checked that all vr_lanes are in this submission
    }

    my $submission = $vr_lane->submission($submission_alias);
    unless ($submission){
        $submission = $vr_lane->add_submission($submission_alias);
        $submission_alias =~ /sc-(20\d{6})/;   # get date from name
        my $sub_date = $1;
        $submission->acc($submission_acc);
        $submission->date($sub_date) if $sub_date;
        $submission->update;
    }
    # update existing submission acc
    unless ($submission->acc){
        $submission->acc($submission_acc);
        $submission->update;
    }

    if ($vr_lane->acc){
        if ($vr_lane->acc ne $lane_acc){
            die "Lane has acc of ",$vr_lane->acc," rather than $lane_acc\n";
        }
    }
    $vr_lane->acc($lane_acc);

    $vr_lane->update;
}

print "ERA submission was successful\n";

###############################################################################

sub make_submission {
    my ($subfile, $expfile, $runfile, $url, $dumpfile) = @_;
    my $cmd = qq(curl -s -k -F "SUBMISSION=\@$subfile" -F "EXPERIMENT=\@$expfile" -F "RUN=\@$runfile" '$url');
    print "$cmd\n" if $verbose;
    my $DUMPFILE;
    if ($dumpfile){
        open $DUMPFILE, ">$dumpfile" or die "Can't open $dumpfile for writing: $!\n";
    }

    my $xmlhash;
    eval{
        my $content = `$cmd`;
        if ($dumpfile){
            print $DUMPFILE "$content\n";
            close $DUMPFILE;
        }
        my $xml = XML::Simple->new(KeyAttr => [], ForceArray=>1);
        $xmlhash = $xml->XMLin($content);
    };
    if($@){
        die "[[XML ERROR]] $cmd: $@\n";
    }
    return $xmlhash;

}

sub submission_was_successful {
    my $xml = shift;
    my $success = $xml->{'success'};
    die "Can't retrieve submission success\n" unless defined $success;
    return $success eq 'true' ? 1 : 0;
}

sub get_submission_alias_acc_from_xml {
    my $xml = shift;
    my $submission_alias = $xml->{SUBMISSION}->[0]->{'alias'};
    my $submission_acc = $xml->{SUBMISSION}->[0]->{'accession'};
    die "Can't retrieve submission info\n" unless ($submission_alias && $submission_acc);
    return ($submission_alias, $submission_acc);
}

sub get_lanes_from_xml {
    my $xml = shift;
    my %lanes;
    foreach my $run (@{$xml->{'RUN'}}){
        $lanes{$run->{alias}} = $run->{accession};
    }

    return \%lanes;
}
