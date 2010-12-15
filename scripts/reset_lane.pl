#!/usr/bin/env perl
# 
# WARNING: Be very careful using this script. Files from lane directories 
# 	will be deleted
# 
# Given file of lane_ids and the path to the META directory will
# reset the directory, i.e. set all processed flags, except
# 'stored', to zero. Then removes all files in the directory,
# if it exists.
# 
# Usage:
# 	reset_lane.pl --meta /path/to/META/ --ids file.of.lane.ids --g1k
# 
# Author: Shane McCarthy, sm15@sanger.ac.uk
# 

use strict;
use warnings;

use VertRes::Utils::VRTrackFactory;
use VRTrack::Lane;
use VertRes::Utils::FileSystem;
use File::Path;
use File::Spec;
use Getopt::Long;

# get user input
my ($help, $g1k_mode, $meta_path, $file_of_ids, $database, $verbose);
GetOptions('g1k'         => \$g1k_mode,
           'meta=s'     => \$meta_path,
           'database=s'  => \$database,
           'ids=s' => \$file_of_ids,
           'verbose'     => \$verbose,
           'h|help'      => \$help);

if ($g1k_mode) {
    $database ||= 'g1k_meta';
}

my $missing_opts = 0;
unless ($file_of_ids && $database && $meta_path) {
    $missing_opts = 1;
}

($help || $missing_opts) and die <<USAGE;

WARNING: Be very careful using this script. Files from lane directories 
	will be deleted

Given file of lane_ids and the path to the META directory will
reset the directory, i.e. set all processed flags, except
'stored', to zero. Then removes all files in the directory,
if it exists.

Usage: $0 --ids file.of.lane.ids
        --meta        path to meta directory
        --g1k         sets database to g1k_meta
        --database    <database name>
        --verbose     be extra verbose about what it is doing
        --help        <this message>

USAGE

my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $database, mode => "rw"); 

unless ($vrtrack) {
    die "DB connection failed: ".$DBI::errstr."\n";
}

die "META directory $meta_path does not exist\n" unless (-e $meta_path);

open my $fh, "<$file_of_ids" || die "Could not open filehandle for file $file_of_ids\n";;

while (<$fh>) {
	next if /^(\s)*$/; # skip blank lines
	chomp;
	
	# Get the lane 
	my $lane_id = $_;
	my $lane = VRTrack::Lane->new_by_hierarchy_name($vrtrack, $lane_id) || die "Could not get lane $lane_id\n";
	my $hier_path = $vrtrack->hierarchy_path_of_lane($lane);
	my $lane_path = File::Spec->catdir($meta_path, $hier_path);
	
	# Warn if lane path does not exist
	die "$lane_path does not exist\n" unless ( -e $lane_path );
	
	my %flags = $lane->allowed_processed_flags();
    
	# Reset the processed flags
	my $number_unprocessed = 0;
	foreach my $flag (keys %flags) {
		next if $flag eq 'stored';
		$lane->is_processed($flag, 0);
		$number_unprocessed++ unless ( $lane->is_processed($flag) );
	}
	$lane->update;
	
	# Check to see that all flags but the 'stored' flag have been reset
	my $reset = $number_unprocessed == (scalar keys %flags) - 1 ? 1 : 0;
		
	# If database has been reset, remove the lane directory
	# and remake an empty directory in its place.
	if ($reset) {
		print "Processed flags for $lane_id have been reset. " if ($verbose);
		print "Removing files from lane directory $lane_path...\n" if ($verbose);
		
		my $start_dir = File::Spec->curdir();
		chdir($lane_path) || die "Could not change in to directory $lane_path\n";
		my $lane_dir = File::Spec->curdir();
		opendir (my $dh, $lane_dir) || die "Could not open directory $lane_dir\n";
		while (my $item = readdir $dh) {
			next if ($item =~ /^\.{1,2}$/); # skip the ./ and ../ directories
			if (-d $item) {
				my $fsu = VertRes::Utils::FileSystem->new();
				my $dir = File::Spec->catdir($lane_dir, $item);
				$fsu->rmtree( $dir );
			} else {
				unlink $item;
			}
		}
		closedir $dh;
		chdir($start_dir);
	}
}
close $fh;
