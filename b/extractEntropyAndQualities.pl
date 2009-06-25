#!/usr/bin/env perl

# Create entquals.txt for every lane in the v1 hierarchy.
# These files can be checked for problems with the data,
# i.e. if we have low qualities or (discounted) entropies much below 2.
#
# This script ought really to be a g1k.pl action...

use Utility;
use File::Basename;
use strict;

foreach my $d (glob("$ENV{G1K}/v1/*/*/lib*/lane*")) {
    my $tgt = "$d/entquals.txt";
    (my $srcD = $d) =~ s/v1/DATA/;
    my ($src1) = glob("$srcD/*.entropy");
    my ($src2) = glob("$srcD/*.readlengths");
    my $src3 = "$d/qualmapBayesian.txt";
    if (makeStatus([$src1,$src2,$src3],[$tgt]) > 0 && stamp($tgt)) {
	report("$tgt: creating");
	extractEntropyAndQualities($src1,$src2,$src3,$tgt);
	report();
	unstamp($tgt);
    }
}

sub extractEntropyAndQualities {
    my ($fEnt,$fRL,$fQual,$tgt) = @_;
    my $fh = openToRead($fRL);
    chomp ($_ = <$fh>);
    my ($len1,$len2) = split;
    $fh = openToRead($fEnt);
    my @ent;
    while (<$fh>) {
	chomp;
	my ($i,$h) = split;
	$i++; # to make it a read position
	if ($len2 == 0) {
	    $ent[0][$i] = $h;
	} elsif ($i <= $len1) {
	    $ent[1][$i] = $h;
	} else {
	    $ent[2][$i-$len1] = $h;
	}
    }
    $fh = openToRead($fQual);
    my (@uncal,@cal);
    while (<$fh>) {
	my ($r,$p,$q1,$q2,$n) = split;
	$uncal[$r][$p][0] += $n;
	$uncal[$r][$p][1] += $n*$q1;
	$cal[$r][$p][0] += $n;
	$cal[$r][$p][1] += $n*$q2;
    }
    my $gh = openToWrite($tgt);
    foreach my $r (0 .. $#ent) {
	next unless ($ent[$r]);
	my $pEnt = $ent[$r];
	foreach my $p (1 .. $#$pEnt) {
	    my $mu = $uncal[$r][$p][0] ? $uncal[$r][$p][1]/$uncal[$r][$p][0] : 0;
	    my $mc = $cal[$r][$p][0] ? $cal[$r][$p][1]/$cal[$r][$p][0] : 0;
	    printf($gh "%d %2d %6.4f %6.4f %6.4f\n",
		   $r,$p,$pEnt->[$p],$mu,$mc);
	}
    }
    $gh->close();
}

    
