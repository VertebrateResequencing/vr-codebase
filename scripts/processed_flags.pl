#!/usr/bin/env perl
use strict;
use VRTrack::Core_obj;

my $val = shift;

unless (defined $val){
    die "Usage: $0 processed_value\nOutputs processed states for a given VRTrack processed value\ne.g. $0 135\n";
}
my %flags = VRTrack::Core_obj->allowed_processed_flags();
my @flags = sort keys %flags;
my $out = join "\n", grep {$val & $flags{$_}} @flags;
print "$out\n" if $out;
