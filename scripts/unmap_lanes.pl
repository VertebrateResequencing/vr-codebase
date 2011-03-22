#!/usr/bin/env perl
#
# unmap_lanes.pl --lanes lanes.fofn --database g1k
#
# Prepares lanes that have previously been mapped to be mapped again. Typically
# used prior to mapping with a different mapper.
#
# Author: Sendu Bala <bix@sendu.me.uk>
#

use strict;
use warnings;
use Getopt::Long;
use VertRes::IO;
use VertRes::Utils::VRTrackFactory;
use VertRes::Utils::Sam;
use File::Basename;
use File::Spec;

# get user input
my ($help, $fofn, $single_lane, $lanes_file, $database, $reset_improvement);
my $spinner = 0;
GetOptions('fofn=s'     => \$fofn,
           'lane=s'     => \$single_lane,
           'lanes=s'     => \$lanes_file,
           'database=s' => \$database,
           'reset_improvement' => \$reset_improvement,
           'h|help'     => \$help);

my $missing_opts = 0;
unless (($fofn || $single_lane || $lanes_file) && $database) {
    $missing_opts = 1;
}

($help || $missing_opts) and die <<USAGE;
Prepares lanes that have previously been mapped to be mapped again. Typically
used prior to mapping with a different mapper.

Contrary to the name of this script, no actual "unmapping" occurs. It simply
unsets the mapped status on lanes in the database. After running this script,
run the mapping pipeline as normal.

Usage: $0 --fofn lanes.fofn --database my_vrtrack_meta
        --fofn        <path> A file with absolute paths to mapped lane dirs
        --lane        <name> The name of a single name
        --lanes       <path> A file with a list of lane names, one per line
        --database    <database name>
        --reset_improvement  if supplied, a lane will become unimproved as well
        --help        <this message>

USAGE


# connect to the database
my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $database,
                                                          mode => 'rw');
unless ($vrtrack) {
    die "DB connection failed: ".$DBI::errstr."\n";
}

my @lane_names;

# read the fofn
if ($fofn) {
    my $vio = VertRes::IO->new;
    my @lanes = $vio->parse_fod($fofn);
    
    my $su = VertRes::Utils::Sam->new;
    
    foreach my $lane_dir (@lanes) {
        my $lane_name = basename($lane_dir);
        push(@lane_names, $lane_name);
    }
}
if ($single_lane) {
    push(@lane_names, $single_lane);
}
if ($lanes_file) {
    open(my $fh, $lanes_file) || die "Could not open lanes file '$lanes_file'\n";
    while (<$fh>) {
        chomp;
        push(@lane_names, $_);
    }
}

# for each lane
foreach my $lane_name (@lane_names) {
    # update the database
    my $lane = VRTrack::Lane->new_by_name($vrtrack, $lane_name);
    unless ($lane) {
        warn "lane '$lane_name' was not found in database '$database'\n";
        next;
    }
    
    $lane->is_processed(mapped => 0);
    $lane->is_processed(improved => 0) if $reset_improvement;
    $lane->update;
}

exit;
