#!/usr/bin/env perl

use strict;
use warnings;
no warnings 'uninitialized';
use Getopt::Long;
use VertRes::Utils::VRTrackFactory;
use VertRes::Wrapper::iRODS;
use VRTrack::Lane;
use VRTrack::Library;
use VRTrack::Sample;
use VRTrack::Individual;

my ($spp, $single_proj, $submit_as, $nosample, $nostudy, $libsource, $libselection, $libstrategy, $machine, $dryrun, $lanes, $holdDate, $do_swaps, $pilot2main, $help);

GetOptions(
    'd|dryrun'      =>  \$dryrun,
    'project=s'     =>  \$single_proj,
    'submit_as=s'   =>  \$submit_as,
    'species=s'     =>  \$spp,
    'nosample'	    =>  \$nosample,
    'nostudy'	    =>  \$nostudy,
    'source=s'	    =>  \$libsource,
    'selection=s'   =>  \$libselection,
    'strategy=s'    =>  \$libstrategy,
    'machine=s'     =>  \$machine,
    'doswaps'       =>  \$do_swaps,
    'p2m'           =>  \$pilot2main, # submit all g1k pilot lanes as appropriate main project by population
    'holduntil=s'   =>  \$holdDate, #must be in YYYY-MM-DD format
    'l|lanes=s'	    =>  \$lanes, #an optional list of lane names to submit
    'h|help'        =>  \$help,
    );

$libsource ||= 'GENOMIC';
$libselection ||= 'RANDOM';
$libstrategy ||= 'WGS';
$machine ||= 'Illumina Genome Analyzer II';

$spp = lc($spp);
my %db_for_spp = ( 'mouse'  => 'mouse_reseq_track',
		    'g1k'   => 'g1k_track',
		  );

my $db = $db_for_spp{$spp};

# only do pilot2main on g1k lanes
$pilot2main = undef unless $spp eq 'g1k';
$submit_as = undef if $pilot2main;  # p2m overrides submit_as

if ($submit_as && !$single_proj){
    warn "Need project set if submit_as is used\n";
    $help++;
}

($db && !$help) or die <<USAGE;
    Usage: $0
        --species    <VRTrack to work with, i.e. g1k or mouse>
        [--project   <limit submission to a single project name>]
        [--submit_as <submit project as this project instead. Requires --project to be set to source project>]
        [--source    <library source, default: GENOMIC>]
        [--selection <library selection method, default: RANDOM>]
        [--strategy  <library strategy, default: WGS>]
        [--machine   <sequencing machine, e.g. "Illumina Genome Analyzer II"]
        [--nosample  <don't make sample XML - for existing projects>]
        [--nostudy   <don't make study XML - for existing projects>]
        [--dryrun    <make files, but don't update database]
        [--doswaps   <submit any lanes with the wrong genotype as the sample they appear to be>]
        [--holduntil <specifies a future date (YYYY-MM-DD) to the archive to keep the files private until>]
        [--lanes     <list of lane srf to limit submission to, in same format as this script produces as output>]  
        --help	     <this message>

Generates a minimal set of submission files for the ERA/SRA.

Works directly from the VRTrack databases, and writes submission info back to
them.

--project and --submit_as work together (--submit_as requires that --project is
set).  This is used for, e.g. submitting g1k pilot lanes as main project
lanes.  If --submit_as is omitted, --project is submitted to --project.

--doswaps assumes that every lane which has a genotype of wrong (and is qc passed) is a sample swap, and so should be submitted to the sample that it was found to be, rather than the sample which its library comes from.  Note that this also associates the library of that lane with the new sample.

Writes XMLs for every lane that is qc passed and has not been submitted
previously.  Also prints to STDOUT a list of srfs for all the lanes involved.

Builds one set of submission files for each study.

Note that this is extremely simple, and assumes paired-end Illumina
sequencing.  It should produce the basic & minimal study, sample, experiment,
run, and submission XMLs required to submit to the ERA.

USAGE

$|++;

our $CENTER = 'SC';
our $center = 'sc';
our $DOMAIN      = 'sanger.ac.uk';
our %MACHINECODE = ('Illumina Genome Analyzer'	    => '1',
		    'Illumina Genome Analyzer II'   => '2',
	            'unspecified'		    => '0',
		  );

my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $db,
                                                          mode => 'rw');

unless ($vrtrack){
    die "Can't connect to tracking database\n";
}

#sanity check the hold date
if( $holdDate )
{
	die("Invalid holduntil data provided (YYYY-MM-DD): $holdDate") unless $holdDate =~ /^(\d\d\d\d)-(\d\d)-(\d\d)$/;
	my ($year, $month, $day) = split( /-/, $holdDate );
	die('Invalid holduntil data provided (YYYY-MM-DD)') unless $day > 0 && $day < 31 && $month > 0 && $month < 13 && $year > 2009;
}

# Build list of lanes to submit
# If we have a list of lanes, limit to those

my %limit_to_lanes;
if ($lanes){
    open (my $LANES, $lanes) or die "Can't open $lanes: $!\n";
    while (<$LANES>){
        chomp;
        $limit_to_lanes{$_}++;
    }
    close $LANES;
}


my %laneinfo;
my %expinfo;
my %sampleinfo;
my %studyinfo;

my @projects;
if ($single_proj){
    my $proj = $vrtrack->get_project_by_name($single_proj);
    @projects = ($proj);
}
else {
    @projects = @{$vrtrack->projects};
}

my $submit_proj;
if ($submit_as){
    $submit_proj = $vrtrack->get_project_by_name($submit_as);
    unless ($submit_proj){
        die "Can't get submit_as project $submit_as\n";
    }
}
    
my $subdate = dateStr(time);

# 2010-08-11 jws
# Right, change this all around.  Rather than descend each project and store
# info on the way down, we're going to get all the lanes that need submitting
# and ascend as required.  This allows us to handle sample swaps, i.e. submit
# as a sample that is not the sample that the library is attached to.

my $dbh = $vrtrack->{_dbh};
my $sql = q[SELECT lane_id FROM latest_lane WHERE qc_status="passed" AND submission_id IS NULL];
my $lane_ids = $dbh->selectall_arrayref($sql);

foreach my $lane_id(map {$_->[0]} @$lane_ids){
    my $lane = VRTrack::Lane->new($vrtrack, $lane_id);
    $lane or die "No VRTrack::Lane for $lane_id\n";
    
    my ($file, $type, $md5) = get_lane_file($lane->name);
    unless ($file){
        warn "No file for $lane, skipping\n";
        next;
    }
    my $cycles = $lane->read_len;
    my $slxcode = $MACHINECODE{$machine};
    die "Unrecognised platform: $machine" unless defined $slxcode;
    # if we've been passed a list of lanes, only do those
    if (%limit_to_lanes){
        next unless $limit_to_lanes{$file};
    }

    my $library = VRTrack::Library->new($vrtrack, $lane->library_id);
    $library or die "No VRTrack::Library off lane ",$lane->name,"\n";

    # Detect swaps.  To tell if it's a sample swap, we'd need to check 
    # all lanes off this library, but for now let's assume all swaps with
    # qc_passed are sample swaps.  Library swaps would be difficult to handle
    # (you don't know what library the lane actually came from, so would
    # probably have to fake a library on the right sample).

    my $gt_status = $lane->genotype_status;
    my $mapstats = $lane->latest_mapping;
    my $gt_exp = $mapstats->genotype_expected;
    my $gt_found = $mapstats->genotype_found;
    my $individual; # set this appropriately
    my $study;    # ditto
    if ($gt_status eq 'wrong'){
        unless ($do_swaps){
            die "lane ",$lane->name," is from sample $gt_found not $gt_exp.  If this is a sample swap, you can submit this as $gt_found using --doswaps\n";
        }
        # right, make the swap.
        $individual = VRTrack::Individual->new_by_name($vrtrack, $gt_found);

        my $sql = qq[select distinct p.study_id as study_id from latest_project p, allocation a, seq_centre sc where a.study_id = p.study_id and a.individual_id=? and a.seq_centre_id=sc.seq_centre_id and sc.name="$CENTER"];
        my $study_ids = $dbh->selectall_arrayref($sql,{Slice => {}},$individual->id);
        if (scalar @$study_ids > 1){
            die "More than one matching study for individual $gt_found\n";
        }
        $study = VRTrack::Study->new($vrtrack, $study_ids->[0]{'study_id'});
        $study or die "No VRTrack::Study off individual $gt_found\n";
    }
    else {
        # do a quick check that the genotype and sample match
        my $sample = VRTrack::Sample->new($vrtrack, $library->sample_id);
        $sample or die "No VRTrack::Sample off library ",$library->name,"\n";
        $individual = $sample->individual;
        if ($gt_exp ne 'none' && ($individual->alias ne $gt_exp)){
            die "lane ",$lane->name," is on sample $individual but the genotype check expected $gt_exp.  NB this is not a swap - the genotype check has been performed wrongly\n";
        }
        my $project = VRTrack::Project->new($vrtrack, $sample->project_id);
        $project or die "No VRTrack::Project off sample ",$sample->name,"\n";
        $study=$project->study;
    }

    # mess about with study to submit to
    my $orig_acc;
    if ($pilot2main){
        if ($study->acc eq 'SRP000031'){    # Lowcoverage pilot 1
            my $pop = $individual->population->name;
            $pop = 'TOS' if $pop eq 'TSI'; # stupid name change
            my $main_pname = "1000Genomes-B1-$pop";
            my $main_proj = $vrtrack->get_project_by_name($main_pname);
            unless ($main_proj){
                die "Can't get pilot2main project $main_pname\n";
            }
            $orig_acc = $study->acc;
            $study = $main_proj->study;
        }
    }
    elsif ($submit_proj){
        $orig_acc = $study->acc;
        $study = $submit_proj->study;
    }
    # OK, at this point we have a lane, a library, an individual and a study object
    # construct data we want to pass to XML builds    
    # print join "\t",($lane->name,$individual->name, $study->acc,"\n");

    # STUDY
    my $study_acc = $study->acc;
    my $projname;
    if ($orig_acc){
        $projname = join '-',($spp,$study_acc,$orig_acc);
    }
    else {
        $projname = join '-',($spp,$study_acc);
    }
    $studyinfo{$study_acc} = $projname;

    # EXPERIMENT 
    my $sample_name = $individual->name;
    my $sample_acc = $individual->acc;
    my $species = $individual->species->name;

    my $lib_name = $library->name;
    my $lib_insert = $library->insert_size;
    # need library name to match existing G1K libraries via David
    # Carter's script,
    $lib_name =~ s/[_ ]/-/g;
    $lib_name = sprintf("%s-%s-%s",$spp,$center,$lib_name);

    my $run = $lane->name;
    $run =~ s/\_.*//;
    my $experiment_id = sprintf("%s-%s-%s", $lib_name,$run,$subdate);

    $sampleinfo{$study_acc}{$sample_name}{species} = $species;
    
    $laneinfo{$study_acc}{$file} = {'filetype'  => $type,
                                    'md5'       => $md5,
                                    'experiment'=> $experiment_id,
                                    'lane'      => $lane,
                                    };


    unless (defined $expinfo{$study_acc}{$experiment_id}){
        $expinfo{$study_acc}{$experiment_id} = {	
                                        'sample'    => $sample_name,
                                        'sampleacc' => $sample_acc,
                                        'species'   => $species,
                                        'library'   => $lib_name,
                                        'insert_size' => $lib_insert,
                                        'cycles'    => $cycles,
                                        'machine'   => $machine,
                                      };
    }
}

foreach my $study(keys %studyinfo){
    # Submission name:
    my $projname = $studyinfo{$study};
    my $subname = lc(sprintf("%s-%s-%s", $projname, $center ,$subdate));

    createStudyXML($study, $projname, $subname) unless $nostudy;
    createSampleXML($sampleinfo{$study}, $subname) unless $nosample;
    createExperimentXML($expinfo{$study}, $study, $subname, $libsource, $libstrategy, $libselection);
    createRunXML($laneinfo{$study}, $subname);
    createSubXML($laneinfo{$study}, $subname, $holdDate );

    my @files_to_submit;
    foreach my $file(keys %{$laneinfo{$study}}){
        my $lane = $laneinfo{$study}{$file}{lane};
        my $submission = $lane->submission;
        if ($submission){
            # double-check here - shouldn't happen as the lanes are already
            # selected as not having a submission
            die $lane->name, " has a submission already!: ",$submission->name;
        }
        else {
            push @files_to_submit, $file;
            if (! $dryrun){
                $submission = $lane->submission($subname);
                unless($submission){
                    $submission = $lane->add_submission($subname);
                    $subname =~ /sc-(20\d{6})/;   # get date from name
                    my $sub_date = $1;
                    $submission->date($sub_date) if $sub_date;
                    $submission->update;
                }
                $lane->update
            }
        }
    }
    print join "\n", @files_to_submit;
    print "\n";
}


sub createStudyXML {
    my ($studyname, $projname, $subname) = @_;
    open (my $STUDY, ">$subname.study.xml") or die "Can't open $subname.study.xml : $!\n";
    my $studyxml = <<'EOXML';
<?xml version="1.0" encoding="UTF-8"?>
<STUDY_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
     <STUDY alias="%s">
          <DESCRIPTOR>
               <STUDY_TITLE></STUDY_TITLE>
               <STUDY_TYPE existing_study_type="Whole Genome Sequencing"/>
               <CENTER_NAME>%s</CENTER_NAME>
               <CENTER_PROJECT_NAME>%s</CENTER_PROJECT_NAME>
               <PROJECT_ID>0</PROJECT_ID>
          </DESCRIPTOR>
     </STUDY>
</STUDY_SET>
EOXML

    printf $STUDY $studyxml, ($studyname, $CENTER,$projname);
    close ($STUDY);
}

sub createSampleXML {
    my ($sampleref,  $subname) = @_;
    open (my $SAMPLE, ">$subname.sample.xml") or die "Can't open $subname.sample.xml : $!\n";
    print $SAMPLE qq(<?xml version="1.0" encoding="UTF-8"?>\n);
    print $SAMPLE qq(<SAMPLE_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n);

    my $sampxml = <<'EOXML';
    <SAMPLE alias="%s">
          <SAMPLE_NAME>
               <COMMON_NAME>%s</COMMON_NAME>
          </SAMPLE_NAME>
     </SAMPLE>
EOXML

    foreach my $samplename (sort keys %{$sampleref}){
	my $species = $sampleref->{$samplename}{species};
	printf $SAMPLE $sampxml, ($samplename,$species);
    }
    
    print $SAMPLE "</SAMPLE_SET>\n";

    close ($SAMPLE);
}

sub createExperimentXML {
    my ($expref, $studyname, $subname, $libsource, $libstrategy, $libselection) = @_;
    open (my $EXP, ">$subname.experiment.xml") or die "Can't open $subname.experiment.xml : $!\n";
    print $EXP qq(<?xml version="1.0" encoding="UTF-8"?>\n);
    print $EXP qq(<EXPERIMENT_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n);
    my $expxml = <<'EOXML';
     <EXPERIMENT alias="%s">
          <STUDY_REF %s="%s"/>
          <DESIGN>
               <DESIGN_DESCRIPTION>Illumina sequencing of %s %s library</DESIGN_DESCRIPTION>
               <SAMPLE_DESCRIPTOR %s="%s"/>
               <LIBRARY_DESCRIPTOR>
                    <LIBRARY_NAME>%s</LIBRARY_NAME>
                    <LIBRARY_STRATEGY>%s</LIBRARY_STRATEGY>
                    <LIBRARY_SOURCE>%s</LIBRARY_SOURCE>
		    <LIBRARY_SELECTION>%s</LIBRARY_SELECTION>
                    <LIBRARY_LAYOUT>
                         <PAIRED NOMINAL_LENGTH="%s"/>
                    </LIBRARY_LAYOUT>
               </LIBRARY_DESCRIPTOR>
               <SPOT_DESCRIPTOR>
                    <SPOT_DECODE_SPEC>
			 <NUMBER_OF_READS_PER_SPOT>2</NUMBER_OF_READS_PER_SPOT>
			 <READ_SPEC>
                              <READ_INDEX>0</READ_INDEX>
                              <READ_CLASS>Application Read</READ_CLASS>
                              <READ_TYPE>Forward</READ_TYPE>
                              <BASE_COORD>1</BASE_COORD>
                         </READ_SPEC>
                         <READ_SPEC>
                              <READ_INDEX>1</READ_INDEX>
                              <READ_CLASS>Application Read</READ_CLASS>
                              <READ_TYPE>Reverse</READ_TYPE>
                              <BASE_COORD>%s</BASE_COORD>
                         </READ_SPEC>
                    </SPOT_DECODE_SPEC>
               </SPOT_DESCRIPTOR>
          </DESIGN>
          <PLATFORM>
               <ILLUMINA>
                    <INSTRUMENT_MODEL>%s</INSTRUMENT_MODEL>
                    <CYCLE_SEQUENCE/>
                    <CYCLE_COUNT>%s</CYCLE_COUNT>
               </ILLUMINA>
          </PLATFORM>
          <PROCESSING>
               <BASE_CALLS>
                    <SEQUENCE_SPACE>Base Space</SEQUENCE_SPACE>
                    <BASE_CALLER>Solexa primary analysis</BASE_CALLER>
               </BASE_CALLS>
               <QUALITY_SCORES qtype="other">
                    <QUALITY_SCORER>Solexa primary analysis</QUALITY_SCORER>
                    <NUMBER_OF_LEVELS>80</NUMBER_OF_LEVELS>
                    <MULTIPLIER>1</MULTIPLIER>
               </QUALITY_SCORES>
          </PROCESSING>
     </EXPERIMENT>
EOXML

    foreach my $exp(sort keys %{$expref}){
	my $sample;
        my $samplelink;
        if ($expref->{$exp}{'sampleacc'}){
            $sample = $expref->{$exp}{'sampleacc'};
            $samplelink = 'accession';
        }
        else {
            $sample = $expref->{$exp}{'sample'};
            $samplelink = 'refname';
        }
	my $species = $expref->{$exp}{'species'};
	my $cycles = $expref->{$exp}{'cycles'};
	my $library = $expref->{$exp}{'library'};
	my $lib_insert = $expref->{$exp}{'insert_size'};
	my $machine = $expref->{$exp}{'machine'};
	my $studylink;	# link by name or by accession?
	if ($studyname =~ /^[ES]RP\d+$/){
	    $studylink = 'accession';
	}
	else {
	    $studylink = 'refname';
	}

	printf $EXP $expxml, ($exp,$studylink,$studyname,$species,$library,$samplelink,$sample,$library,$libstrategy, $libsource, $libselection, $lib_insert, $cycles + 1, $machine, $cycles);
    }

    print $EXP "</EXPERIMENT_SET>";
    close ($EXP);
}


sub createRunXML {
    my ($laneref, $subname) = @_;
    open (my $RUN, ">$subname.run.xml") or die "Can't open $subname.run.xml : $!\n";
    print $RUN qq(<?xml version="1.0" encoding="UTF-8"?>\n);
    print $RUN qq(<RUN_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n);
    my $runxml = <<'EOXML';
     <RUN alias="%s" total_data_blocks="1">
          <EXPERIMENT_REF refname="%s"/>
          <DATA_BLOCK name="%s">
               <FILES>
                    <FILE filename="%s" filetype="%s"/>
               </FILES>
          </DATA_BLOCK>
     </RUN>
EOXML

    my %exp_count;
    foreach my $lane(sort keys %{$laneref}){
	my $filetype= $laneref->{$lane}{'filetype'};
	my $blockname = $lane;
	$blockname =~ s/\..*//;	# trim off type
	my $exp  = $laneref->{$lane}{'experiment'};
	$exp_count{$exp}++; 
	my $index = $exp_count{$exp}; 
        # 2010-02-22 jws 
        # Need to change run name as experiments could collide.  A good run name
        # is simply the 1234_1 name as this is unique.  So, blockname.
	# use to be: my $runname = "${exp}_$index";
        my $runname=$blockname;

	printf $RUN $runxml, ($runname,$exp,$blockname,$lane,$filetype);
    }
    print $RUN "</RUN_SET>";
    close ($RUN);
}


sub createSubXML {
    my ($laneref, $subname,$holdDate) = @_;
    open (my $SUB, ">$subname.submission.xml") or die "Can't open $subname.submission.xml : $!\n";

    # header
    print $SUB qq(<?xml version="1.0" encoding="UTF-8"?>\n);
    print $SUB sprintf (qq(<SUBMISSION center_name="%s" submission_date="%s" submission_id="%s" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n), $CENTER,isoDate(time),$subname);

    # Contact info
    my $conxml = <<'EOXML';
  <CONTACTS>
    <CONTACT inform_on_error="%s" inform_on_status="%s" name="%s" />
  </CONTACTS>
EOXML

    my $userdata = getUserData();
    my $email = $userdata->{mail};
    my $name = $userdata->{name};
    printf $SUB $conxml, ($email, $email, $name);

    # reference the other XML files
    my $actxml =  <<'EOXML';
    <ACTION>
      <ADD schema="%s" source="%s" />
    </ACTION>
EOXML

    print $SUB "  <ACTIONS>\n";
    unless($nostudy){
	printf $SUB $actxml, "study", "$subname.study.xml";
    }
    unless ($nosample){
	printf $SUB $actxml, "sample", "$subname.sample.xml";
    }
    foreach $_ qw(experiment run){
	printf $SUB $actxml, $_, "$subname.$_.xml";
    }
	
	if( $holdDate )
	{
		print $SUB qq[
      <ACTION>
        <HOLD HoldUntilDate="$holdDate"/> 
      </ACTION>
	  ];
	}
    print $SUB "  </ACTIONS>\n";
	
    # Files
    my $filexml =  <<'EOXML';
    <FILE filename="%s" checksum="%s" checksum_method="MD5" />
EOXML

    print $SUB "  <FILES>\n";
    foreach my $lane(sort keys %{$laneref}){
	my $md5= $laneref->{$lane}{'md5'};
	printf $SUB $filexml, ($lane,$md5);
    }
    print $SUB "  </FILES>\n";

    print $SUB "</SUBMISSION>\n";
    close $SUB;
}


sub dateStr {
    my ($t) = @_;
    my ($sec,$min,$hour,$mday,$mon,$year) = localtime($t);
    my $date = sprintf("%4d%.2d%.2d", $year+1900,$mon+1,$mday);
    return $date;
}


sub isoDate {
    my ($t) = @_;
    my ($sec,$min,$hour,$mday,$mon,$year) = localtime($t);
    my $date = sprintf("%4d-%2d-%2dT%2d:%2d:%2dZ",
                      $year+1900,$mon+1,$mday,$hour,$min,$sec);
    $date =~ tr/ /0/;
    return $date;
}


sub getUserData {
    my $u = $ENV{USER};
    my @a = getpwnam($u);
    my $name = $a[6];
    $name =~ s/,.*//;
    return {uid => $u,
            name => $name,
            mail => "$u\@$DOMAIN"};
}

sub get_lane_file {
    # check for srf, if not found, check for bam
    my $lane = shift;
    my $type;
    my $file = "$lane.srf";
    my $md5 = `/software/solexa/bin/mpsa_download -m -f $file`;
    unless ($md5){
        $file =~ s/_/_s_/;
        $md5 = `/software/solexa/bin/mpsa_download -m -f $file`;
    }

    if ($md5){
        ($md5) = split /\s+/,$md5;
        $type = "srf";
    }
    else {
        #die "Can't get md5 for $lane\n";
        warn "no srf for $lane, try bam\n";
        $lane =~ /^(\d+)_(\d+)/;
        my ($run, $pos) = ($1,$2);
        my $irods = VertRes::Wrapper::iRODS->new();
        my @bam = @{$irods->find_files_by_run_lane($run,$pos)};
        if (scalar @bam > 1){
            # TODO: handle multiplex here
            die "More than one bam for $lane\n";
        }
        elsif (scalar @bam  == 1){
            my $bam = $bam[0];
            $file = $bam;
            $file =~ s|.*/||;   # trim off irods path
            $type = "bam";
            $md5 = $irods->get_file_md5($bam);
        }
    }
    unless ($md5){  # return undef if we haven't found it.
        $file = $type = undef;
    }
    return ($file,$type, $md5);
}

