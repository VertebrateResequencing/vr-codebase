#!/usr/bin/env perl
use strict;
use warnings;
no warnings 'uninitialized';
use Getopt::Long;
use VertRes::Utils::VRTrackFactory;
use VRTrack::Lane;
use Sfind::Fastq;
use Date::Manip;
use Data::Dumper;

my ($mindate,$spp, $help);

GetOptions(
    's|spp=s'       =>  \$spp,
    'd|date=s'      =>  \$mindate,
    'h|help'	    =>  \$help,
    );

$mindate ||= "1000-01-01";  # very early min date...
my $since = ParseDate($mindate);

my %db_for_spp = ( 'mouse'  => 'mouse_reseq_track',
		    'g1k'   => 'g1k_track',
	        );

my $db = $db_for_spp{$spp};

($db && $since && !$help) or die <<USAGE;
    Usage: $0   
                --spp       <species, i.e. g1k or mouse>
                [--date     <only include lanes after this date>]
                --help      <this message>

Summarises VRTrack lanes in monthly bins: Gb and mean Q, by read-length.

USAGE

warn "Species: $spp\n";
warn "Database: $db\n";
my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $db,
                                                          mode => 'r');

unless ($vrtrack){
    die "Can't connect to tracking database\n";
}


my $lanes = $vrtrack->qc_filtered_lane_names();
my %lanes_in_bin;
my %readlength;
foreach my $lanename (@$lanes){
    $lanename =~ /^(\d+)_(\d)$/ or die "Can't match $lanename";
    my $run = $1;
    my $lane = $2;
    my $vrlane = VRTrack::Lane->new_by_hierarchy_name($vrtrack,$lanename);
    die "No lane for $lanename" unless $vrlane;
    my $dateflag = Date_Cmp($vrlane->run_date,$since);
    if ($dateflag < 1){ 
        # $run_date is earlier than $since
        next;
    }
    my $rundate = ParseDate($vrlane->run_date());
    my $bin = UnixDate($rundate, '%Y-%m');
    # get readlen for hashkey, add raw bases
    my $readlen = $vrlane->read_len();
    my $bp = $vrlane->raw_bases;

    $readlength{$readlen}++;
    $lanes_in_bin{$bin}{$readlen}{lanes}++;
    $lanes_in_bin{$bin}{$readlen}{bp} += $bp;

    foreach my $fastq (@{$vrlane->files}){
        $lanes_in_bin{$bin}{$readlen}{meanq} += $fastq->mean_q;
        $lanes_in_bin{$bin}{$readlen}{fastq}++;
    }
}

my @bins = sort keys %lanes_in_bin;

my @readlengths = sort {$a <=> $b}  keys %readlength;
# header
print "run";
foreach my $rl(@readlengths,"total"){
    print "\t";
    print join "\t",("$rl lanes","$rl Gb","$rl meanq");
}
print "\n";
foreach my $bin (@bins){
    my $bintotal_lanes = 0;
    my $bintotal_bp = 0;
    my $bintotal_q_total = 0;   # total of q means
    my $bintotal_q_n = 0;   # number of q means to average over

    print $bin;
    foreach my $rl(@readlengths){
        print "\t";
        my $mean_mean_q = 0;
        my $total_bp = 0;
        my $total_lanes = 0;
        if ($lanes_in_bin{$bin}{$rl}){
            $mean_mean_q = sprintf("%.1f",$lanes_in_bin{$bin}{$rl}{meanq}/$lanes_in_bin{$bin}{$rl}{fastq});
            $total_bp = sprintf("%.1f",$lanes_in_bin{$bin}{$rl}{bp}/1e9);
            $total_lanes = $lanes_in_bin{$bin}{$rl}{lanes};
        }
        print join "\t", (  $total_lanes,
                            $total_bp,
                            $mean_mean_q);
        $bintotal_lanes += $lanes_in_bin{$bin}{$rl}{lanes};
        $bintotal_bp += $lanes_in_bin{$bin}{$rl}{bp};
        $bintotal_q_total += $lanes_in_bin{$bin}{$rl}{meanq};
        $bintotal_q_n += $lanes_in_bin{$bin}{$rl}{fastq};
    }
    # totals
    print "\t";
    print join "\t", (  $bintotal_lanes,
                        sprintf("%.1f",$bintotal_bp/1e9),
                        sprintf("%.1f",$bintotal_q_total/$bintotal_q_n)
                        );
    print "\n";
}
