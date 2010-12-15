#!/usr/bin/env perl
use strict;
use warnings;
no warnings 'uninitialized';
use Getopt::Long;
use VertRes::Utils::VRTrackFactory;
use VRTrack::Lane;
use File::Path;

my ($spp, $db, $hroot, $help);

GetOptions(
    's|spp=s'       =>  \$spp,
    'd|db=s'        =>  \$db,
    'r|root=s'      =>  \$hroot,
    'h|help'	    =>  \$help,
    );

my %db_for_spp = ( 'mouse'  => 'mouse_reseq_track',
		    'g1k'   => 'g1k_track',
	  );

# For testing
#my %db_for_spp = ( 'mouse'  => 'mouse_cancer_track',
#		    'g1k'   => 'mouse_cancer_track',
#		  );

$db ||= $db_for_spp{$spp};

($db && !$help) or die <<USAGE;
    Usage: $0 <file of lanes to delete>  
                --spp       <species, i.e. g1k or mouse>
                --db        <specify db name directly, instead of via spp>
               [--root      <root of the hierarchy the lane dir is under>]
                --help      <this message>

Deletes lanes (and associated files and mapstats) from tracking database, and,
if hierarchy root is specified, deletes the lane directory from a file
hierarchy.

USAGE

warn "Species: $spp\n" if $spp;
warn "Database: $db\n";
$db || die "a known --spp or a --db must be specified\n";
my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $db,
                                                          mode => 'rw');

unless ($vrtrack){
    die "Can't connect to tracking database\n";
}

# get lanes from fofn or stdin
while (<>){
    my $lanename = $_;
    chomp $lanename;
    my $lane = VRTrack::Lane->new_by_name($vrtrack,$lanename);
    unless ($lane){
        print "Can't get lane $lanename\n";
        next;
    }
    if ($hroot){
        my $laneh = $vrtrack->hierarchy_path_of_lane_name($lane->name);
        if (-d "$hroot/$laneh"){
            print "Deleting $hroot/$laneh\n";
            rmtree("$hroot/$laneh");
        }
        else {
            print "Can't find lane directory $hroot/$laneh\n";
        }
    }
    if ($lane->delete){
        print "$lanename deleted from database\n";
    }
    else {
        print "Problem deleting $lanename from database\n";
    }
}
