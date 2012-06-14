#!/usr/bin/env perl

#
# delete_lane.pl --root /path/to/TRACKING --db dbname <file-of-lane-names>
# (optional --verbose flag can be set to output actions of the script)
#
# For example:
# delete_lane.pl --root /lustre/scratch106/projects/uk10k/TRACKING 
# --db vrtrack_uk10k_cohort --verbose lanes.txt
#
# Sets up the lanes that need to have the genotype checking performed again. It deletes the lane directory on disk
# and removes the lane and all descendant objects from the tracking database.
#
# Author: Jim Stalker (jws)
#

use strict;
use warnings;
no warnings 'uninitialized';

use Getopt::Long;
use VertRes::Utils::VRTrackFactory;
use VRTrack::Lane;
use File::Path qw(remove_tree);

my ($db, $root, $help, $verbose);

GetOptions(
    'd|db=s'        =>  \$db,
    'r|root=s'		=>  \$root,
    'v|verbose'     =>  \$verbose,
    'h|help'	    =>  \$help,
    );

($db && !$help) or die <<USAGE;
    Usage: $0 <file of lane names to delete>  
                --db        <specify db name>
                [--root      <root directory for the analyses>]
                --verbose   [be extra verbose about what it is doing]
                --help      <this message>

Deletes the lanes specified. Deletes the lane and all descendant objects
(files, mapstats, etc) from the tracking database and, if root is specified,
removes the lane directory from disk.

USAGE

if ($root){
    if ('/' ne substr $root,-1,1) {
            $root = $root.'/';
    }
    die "Root directory $root does not exist\n" unless (-e $root);
}

my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $db,
                                                          mode => 'rw');
unless ($vrtrack){
    die "Can't connect to tracking database\n";
}

# get lanes from file of lane names or stdin
while (<>){
    my $lanename = $_;
    chomp $lanename;
    #Check lane actually exists or bail out
    my $lane = VRTrack::Lane->new_by_name($vrtrack,$lanename);
    unless ($lane){
       print "Can't get lane $lanename\n";
       next;
    }
    if ($root){
        #Delete files first
        #Get full path to lane directory
        my $lanedir = $root.$vrtrack->hierarchy_path_of_lane_name($lane->name);
        print "Deleting: \n" if $verbose;
        remove_tree($lanedir, {verbose => $verbose, safe => 1});
    }

    #update database
    if ($lane->delete){
        print "$lanename deleted from database\n" if $verbose;
    }
    else {
        print "Problem deleting $lanename from database\n";
    }
}
