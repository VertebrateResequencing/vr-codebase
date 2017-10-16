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
use VertRes::Utils::FileSystem;
use File::Path qw(remove_tree);
use Cwd 'abs_path';

my ($db, $root, $clean, $help, $verbose);

GetOptions(
    'd|db=s'        =>  \$db,
    'r|root=s'		=>  \$root,
    'c|clean'		=>  \$clean,
    'v|verbose'     =>  \$verbose,
    'h|help'	    =>  \$help,
    );

($db && !$help) or die <<USAGE;
    Usage: $0 <file of lane names to delete>  
                --db        <specify db name>
                [--root      <root directory for the analyses>]
                --clean     [delete stored folder - in nexsan]
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
	
	next if(length($lanename) < 3 );
	
    #Check lane actually exists or bail out
    my $lane = VRTrack::Lane->new_by_name($vrtrack,$lanename);
    unless ($lane){
       print "Can't get lane $lanename\n";
       next;
    }
    if ($root){
        # We first get the hierarchy structure for this lane. This is usually a symlink. 
        # If -c is specified, we delete any entries in the fsu_file_exists table for this lane
        # and also the folder which the symlink points to. We then delete the symlink itself.
        # We don't rely on the stored_path in the lanes table because it may not be the
        # one that was used to actually store the data.
        # Further down in this script, the record for this lane is deleted from the DB
        # so we do not have to do any updating of the processed flag etc.
        # 5 June 2013
        
        #Get full path to lane directory
		
		my $lane_suffix_dir = $vrtrack->hierarchy_path_of_lane_name($lane->name);
		# If you dont check this exists then you end up deleting the root directory of the pipeline 
		next if(! defined($lane_suffix_dir) || (length($lane_suffix_dir) < 10) );
        my $lanedir = $root.$lane_suffix_dir;
        
        if($clean){
        	
        	#Clear up any files in the FSU FILE EXISTS tables
        	my $fsu = VertRes::Utils::FileSystem->new(reconnect_db=>1);
        	print "Deleting: \n Files in the fsu_file_exists table that have a path like $lanedir \n" if $verbose;
    		$fsu->file_exists($lanedir, recurse =>1, wipe_out=>1);
 
 	        my $stored_path = readlink $lanedir;  #Get the folder pointed to by the symlink

         	print "Deleting: \n" if $verbose;
        	remove_tree($stored_path, {verbose => $verbose, safe => 0}); #Delete folder pointed to by symlink
        }
               
        print "Deleting: \n" if $verbose;
        remove_tree($lanedir, {verbose => $verbose, safe => 1}); #Delete symlink
        
    }

    #update database
    if ($lane->delete){
        print "$lanename deleted from database\n" if $verbose;
    }
    else {
        print "Problem deleting $lanename from database\n";
    }
}
