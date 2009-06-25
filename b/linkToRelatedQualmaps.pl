#!/usr/bin/env perl

# This script can be run to rescue the situation when no qualmap file
# can be generated for a lane, usually because there isn't sufficient
# good data to train one on. On the (nearly always accurate) assumption
# that lanes for the same run will have very similar quality behaviours,
# it looks for such lanes, and links to the first one it finds.

# The script doesn't actually create the links; it outputs "ln" commands
# for you to check and run.

# When there is a qualmap file but it's been created on the basis of too
# few reads (should no longer happen) it's moved out of the way.

use G1K::Calibrate;
use File::Basename;
use strict;

my %dirsForRun;
my $DIR = "$ENV{G1K}/v1";
chdir($DIR);
print "cd $DIR\n";
foreach my $d (glob("*/NA?????/lib*/lane*")) {
    #unless ($d =~ /libBI/) { # postpone while still doing BI
	(my $r = basename($d)) =~ s/_[1-8]$//;
	push @{$dirsForRun{$r}},$d;
    #}
}

foreach my $pA (sort {$a->[0] cmp $b->[0]} values %dirsForRun) {
    my (@good,@bad);
    foreach my $d (@$pA) {
	if (G1K::Calibrate::enoughReadsMapped("$d/sample/aln.map")) {
	    push @good,[$d,basename(dirname($d))] if (-s "$d/qualmapBayesian.txt");
	} else {
	    push @bad,[$d,basename(dirname($d))];
	    if (-f "$d/qualmapBayesian.txt") {
		print "sv -m $d/qualmapBayesian.txt $d/recalibrated*\n";
	    }
	}
    }
    my $done;
    foreach my $pPair1 (@bad) {
	my ($d,$lib) = @$pPair1;
	my $done = 0;
      LOOP:
	foreach my $match (1,0) {
	    foreach my $pPair2 (@good) {
		unless ($match && $pPair2->[1] ne $lib) {
		    print "ln $pPair2->[0]/qualmapBayesian.txt $d/qualmapBayesian.txt\n";
		    $done = 1;
		    last LOOP;
		}
	    }
	}
	unless ($done) {
	    print "# CANNOT LINK FOR $d/qualmapBayesian.txt\n";
	}
    }
}
    
