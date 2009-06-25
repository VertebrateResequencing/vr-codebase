#!/usr/local/bin/perl

use strict;
use warnings;
no warnings 'uninitialized';
use Getopt::Long;
use DBI;

$|++;

my $g1k_dbh = DBI->connect("DBI:mysql:host=mcs4a:port=3306;database=g1k_track", "vreseq_ro",undef,
                             {'RaiseError' => 1, 'PrintError'=>0});
my ($lane_file, $md5_file, $spp, $help);

GetOptions(
    'm|md5|md5s=s'  =>  \$md5_file,
    's|spp|species=s'=>  \$spp,
    'h|help'	    =>  \$help,
    );

(-s $md5_file && $spp && !$help) or die <<USAGE;
    Usage: $0   
                --md5s   <file of md5s for srfs>
                --spp    <species (e.g. 'Homo sapiens') >
		--help   <this message>

Build a tab-delimited submission info file for minimalsubmission.pl to use

**** G1K Main Project specific!!! ****

USAGE

my %accession_from_study = ('TSI' => 'SRP000540',
			    'YRI' => 'SRP000542',
			    'LWK' => 'SRP000543',
			    'JPT' => 'SRP000544',
			    'CHB' => 'SRP000546',
			    'CEU' => 'SRP000547',
			    );

open (my $LANE, $md5_file) or die "Can't open $md5_file: $!\n";

while (<$LANE>){
    chomp;
    my ($md5, $lane) = split (/ +/, $_);
    $lane =~ s/\*//;
    $lane =~ s/_s_/_/;
    $lane =~ s/.srf$//;

    my $sql = qq[select project.name as project_name, sample.name as sample_name, library.name as library_name, readlen from project, sample, library, lane where lane.name=? and project.project_id=sample.project_id and sample.sample_id = library.sample_id and library.library_id = lane.library_id;];
    my $ref = $g1k_dbh->selectrow_hashref($sql, undef, ($lane));
    my ($study, $sample, $library, $cycles);
    if ($ref){
	$study = $ref->{'project_name'};
        $sample = $ref->{'sample_name'};
        $library = $ref->{'library_name'};
        $cycles = $ref->{'readlen'};
    }
    else {
	die "No data for $lane\n";
    }
    $study =~ s/1000.*_//;	# should be left with CEU, TOS, etc
    $study = 'TSI' if $study eq 'TOS';  #TuScans from Italy.  Not Toscans.
    my $project = "g1k-$study";
    my $study_acc = $accession_from_study{$study};
    $study_acc or die "Can't find accession for study $study\n";
    my $srf = "${lane}.srf";

    print join "\t", ($study_acc, $project, $srf,'srf',$md5,$sample,$library,$spp,$cycles,'Illumina Genome Analyzer II');
    print "\n";

}

close $LANE;
