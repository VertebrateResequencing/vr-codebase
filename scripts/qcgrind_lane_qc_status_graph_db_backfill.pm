#!/usr/bin/env perl
use strict;
use warnings;
use VRPipe;
use VRPipe::Schema;
use VRTrack::Factory;

my $database = shift;
chomp($database);

my $vrtrack = VRTrack::Factory->instantiate(database => $database, mode => 'r');
my $schema = VRPipe::Schema->create("VRTrack");

foreach my $name (@{$vrtrack->qc_filtered_lane_names}) {
    my $lane = VRTrack::Lane->new_by_name($vrtrack, $name);
    my $status = $lane->qc_status;
    next if $status eq 'pending';
    
    my $node = $schema->get("Lane", { unique => $name });
    if ($node) {
        $node->qcgrind_qc_status($status);
        print "$name => ", $lane->qc_status, "\n";
    }
}

exit;
