#!/usr/bin/env perl
#
# release_lanes.pl --root /abs/META --db g1k_meta --all_samples --ignore_platform SOLID > lanes.fod
#
# Prints out the paths to all the lane directories that would be needed to
# build sample-level release bams. This file can be used as input to the
# release pipeline (VertRes::Pipelines::Release). It will only output lanes for
# a sample if all that sample's lanes have been mapped (ignoring lanes matching
# an --ignore_platform option).
#
# Author: Sendu Bala <bix@sendu.me.uk>
#

use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use File::Spec;
use VertRes::Utils::Hierarchy;
use VertRes::Utils::VRTrackFactory;
use VRTrack::Lane;
use VRTrack::History;

# get user input
my ($help, $database, $all_samples, @ignore_platforms, @samples, $root, $qc, $date);
my $spinner = 0;
GetOptions('db=s'              => \$database,
           'all_samples'       => \$all_samples,
           'samples:s{1,}'     => \@samples,
           'root=s'            => \$root,
           'ignore_platform=s' => \@ignore_platforms,
           'qc'                => \$qc,
           'date=s'        => \$date,
           'h|help'            => \$help);

my $missing_opts = 0;
unless ($database && ($all_samples || @samples) && $root) {
    $missing_opts = 1;
}

($help || $missing_opts) and die <<USAGE;
Prints out the paths to all the lane directories that would be needed to build
sample-level release bams. This file can be used as input to the release
pipeline (VertRes::Pipelines::Release). It will only output lanes for a sample
if all that sample's lanes have been mapped (ignoring lanes matching an
--ignore_platform option).

Usage: $0 --db my_vrtrack_meta --all_samples --root /abs/path/to/META
        --db               <database name>
        --all_samples      get all lane paths in the database
        --root             root directory of the tracking database
        --qc               get only lanes that are set as qc passed
        --samples          <sample name> [<sample name> ...]
        --ignore_platform  <SLX|454|SOLID> ignore lanes sequenced with this
                                           platform (this whole option can be
                                           specified more than once)
        --date             '2010-01-04 10:49:10' view the database as it was
                                                 immediately prior to this date
        --help             this message

USAGE

if ($all_samples && @samples) {
    warn "Both --all_samples and --samples were set; ignoring the --all_samples request\n";
    $all_samples = 0;
}

# travel back in time?
my $hist = VRTrack::History->new();
$hist->time_travel($date) if $date;

# get all the lanes for our desired samples, excluding ignored platforms
my @platforms;
my %ignore_platforms = map { $_ => 1 } @ignore_platforms;
foreach my $platform ('SLX', '454', 'SOLID') {
    next if exists $ignore_platforms{$platform};
    push(@platforms, $platform);
}

my %cd = VertRes::Utils::VRTrackFactory->connection_details('r');
$cd{database} = $database;

my $hu = VertRes::Utils::Hierarchy->new();
my @lanes = $hu->get_lanes(db => { %cd },
                           $all_samples ? () : (sample => [@samples]),
                           platform => [@platforms]);

# get the paths of all lanes, and determine which samples are fully mapped
# (now, not at $date)
$hist->time_travel('latest');
my $vrtrack = $lanes[0]->vrtrack;

my %bad_samples;
my %lanes_by_sample;
my %seen_paths;
foreach my $lane (@lanes) {
    if ($date) {
        $lane = VRTrack::Lane->new($vrtrack, $lane->id);
    }
    my %objs = $hu->lane_hierarchy_objects($lane);
    my $sample = $objs{sample}->name;
    
    next if exists $bad_samples{$sample};
    next if( $qc && ( !$lane->qc_status() || $lane->qc_status() ne 'passed' ) );
    
    # are we mapped?
    unless ($lane->is_processed('mapped')) {
	my $lane_name = $lane->name;
	warn "$lane_name was not mapped!\n";
        $bad_samples{$sample} = 1;
        delete $lanes_by_sample{$sample};
        warn "Not all the lanes under sample '$sample' were mapped; they will all be excluded\n";
        next;
    }
    
    my $path = File::Spec->catdir($root, $vrtrack->hierarchy_path_of_lane($lane));
    if (exists $seen_paths{$path}) {
        warn "already saw $path...\n";
        next;
    }
    $seen_paths{$path} = 1;
    push(@{$lanes_by_sample{$sample}}, $path);
}

# output good lanes
while (my ($sample, $lanes) = each %lanes_by_sample) {
    foreach my $path (@{$lanes}) {
        print $path, "\n";
    }
}

exit;
