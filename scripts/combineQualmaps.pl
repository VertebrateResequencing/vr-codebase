#!/usr/bin/env perl

# Combine qualmapBayesian.pl files into an averaged one where we recalculate
# each output quality simple-mindedly, just to get an impression of how
# (e.g.) a chunk of BGI data is recalibrating. The output is in qualmap format.

use strict;

my (%N,%E);
while (<>) {
    chomp;
    my ($r,$p,$q1,undef,$n,$e) = split;
    if ($n) {
	my $key = sprintf("%d %2d %2d",$r,$p,$q1);
	$N{$key} += $n;
	$E{$key} += $e;
    }
}
foreach my $key (sort keys %N) {
    my $n = $N{$key} || 0;
    my $e = $E{$key} || 0;
    my $q=20;
    if ($n>0 && $e>0) {
	my $ep = $e/$n;
	$q = int(-10*log($ep)/log(10)+0.5);
    }
    printf("%s %2d %7d %7d\n",$key,$q,$n,$e);
}
