#!/usr/bin/env perl

use strict;
use warnings;
no warnings 'uninitialized';
use Getopt::Long;
use VertRes::Utils::VRTrackFactory;
use DBI;

$|++;

my ($lane_file, $md5_file, $species, $help);

my %spp_config = ('human' => {  'spp' => 'Homo sapiens',
                                'db'  => 'g1k_track',
                                'proj_pref' => 'g1k',
                                'study_acc' => {'TSI' => 'SRP000540',
                                                'YRI' => 'SRP000542',
                                                'LWK' => 'SRP000543',
                                                'JPT' => 'SRP000544',
                                                'CHB' => 'SRP000546',
                                                'CEU' => 'SRP000547',
                                                }
                            },
                 'mouse' => {   'spp'  => 'Mus musculus',
                                'db'  => 'mouse_reseq_track',
                                'proj_pref' => 'mouse',
                                'study_acc' => {
                                               'NZO_Mouse_Genome' => 'ERP000047',
                                                'PWK_Ph_Mouse_Genome' => 'ERP000048',
                                                'BALBc_J_Mouse_Genome' => 'ERP000039',
                                                'DBA_2J_Mouse_Genome' => 'ERP000044',
                                                'WSB_Ei_Mouse_Genome' => 'ERP000050',
                                                '129P2_Mouse_Genome' => 'ERP000034',
                                                'Spretus_Ei_Mouse_Genome' => 'ERP000049',
                                                'C57BL_6N_Mouse_Genome' => 'ERP000041',
                                                'NOD_Mouse_Genome' => 'ERP000046',
                                                '129S1_SvImJ_Mouse_Genome' => 'ERP000035',
                                                'CBA_J_Mouse_Genome' => 'ERP000043',
                                                'A_J_Mouse_Genome' => 'ERP000038',
                                                'C3H_HeJ_Mouse_Genome' => 'ERP000040',
                                                '129S5_Mouse_Genome' => 'ERP000036',
                                                'CAST_Ei_Mouse_Genome' => 'ERP000042',
                                                'AKR_J_Mouse_Genome' => 'ERP000037',
                                                'LP_J_Mouse_Genome' => 'ERP000045',
                                                }
                            }
                 );

GetOptions(
    'm|md5|md5s=s'  =>  \$md5_file,
    's|spp|species=s'=>  \$species,
    'h|help'	    =>  \$help,
    );

$species = lc($species);
my $spp = $spp_config{$species}{'spp'};

(-s $md5_file && $spp && !$help) or die <<USAGE;
    Usage: $0   
                --md5s   <file of md5s for srfs>
                --spp    <species (i.e. 'human' or 'mouse') >
		--help   <this message>

Build a tab-delimited submission info file for minimalsubmission.pl to use

USAGE

my $db = $spp_config{$species}{'db'};

my %cd = VertRes::Utils::VRTrackFactory->connection_details('r');

my $dbh = DBI->connect("DBI:mysql:host=$cd{host};database=$db",
			$cd{user}, undef,  
                      {'RaiseError' => 0, 'PrintError'=>0}
		      );
if ($DBI::err){
      die(sprintf('DB connection failed: %s', $DBI::errstr));
}
open (my $LANE, $md5_file) or die "Can't open $md5_file: $!\n";

while (<$LANE>){
    chomp;
    my ($md5, $lane) = split (/ +/, $_);
    $lane =~ s/\*//;
    $lane =~ s/_s_/_/;
    $lane =~ s/.srf$//;

    my $sql = qq[select project.name as project_name, sample.name as sample_name, library.name as library_name, readlen from project, sample, library, lane where lane.name=? and project.project_id=sample.project_id and sample.sample_id = library.sample_id and library.library_id = lane.library_id;];
    my $ref = $dbh->selectrow_hashref($sql, undef, ($lane));
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
    $study =~ s/1000.*_//;	# should be left with CEU, TOS, etc for G1K
    $study = 'TSI' if $study eq 'TOS';  #TuScans from Italy.  Not Toscans.
    my $project = $spp_config{$species}{'proj_pref'}."-$study";
    my $study_acc = $spp_config{$species}{'study_acc'}{$study};
    $study_acc or die "Can't find accession for study $study\n";
    my $srf = "${lane}.srf";

    print join "\t", ($study_acc, $project, $srf,'srf',$md5,$sample,$library,$spp,$cycles,'Illumina Genome Analyzer II');
    print "\n";

}

close $LANE;
