#!/usr/bin/env perl

use Utility;
use strict;

autoflush STDOUT 1;

# Create *.fastq.entropy files to evaluate data quality. If the
# discounted entropy (see below) is much below 2, the data is likely to be
# problematic. See $RECAL variable below for what is actually analysed.

# Running extractEntropyAndQualities.pl after this script will give
# a fuller picture of what's going on. 

# This information is particularly useful when recalibration fails
# owing to (a) an excess of reads with "N"s in them or (b) not enough
# reads mapping with high quality. If the discounted entropy dips too
# much, you can be reasonably sure it's a data problem not a code one.

# The "discounted entropy" Hdisc of a set of nucleotide values (from
# {A,C,G,T,N}) with counts c(A), ..., c(N), typically all nucleotides at
# a given position in a given lane, is defined as follows:
# 
#  n = c(A)+c(C)+c(G)+c(T)
# 
#  f(x) = c(x)/n, for x=A,C,G,T
# 
#  H = - sum{x=A,C,G,T} f(x) log_2 f(x)
# 
#  Hdisc = H*n/(n+c(N))
# 
# In other words, Hdisc is H multiplied by the proportion of non-N's
# in the set.  When everything is working well, there are very few
# N's, and Hdisc ~ H ~ 1.95 (H would be 2 if A, C G and T were equally
# frequent, but the 40:60 CG:AT ratio pushes it down a bit). When
# there is a major error or systematic bias arising from sequencing
# problems, either one genuine nucleotide will predominate, or there
# will be a lot of N's. Whatever the cause, when Hdisc drops below
# about 1.5 for any position in a lane, or below about 1.85 for
# several positions, there are likely to be severe alignment problems.
#
# Note: really instead of H we should use something that allows
# for G/C being rarer than A/T, i.e. some sort of Kullback-Liebler
# distance-based measure...but it seems to make little 
# practical difference.

my ($dir,$base);

# This should be 1 or 0; or you could introduce a switch to allow the user to
# choose without editing the script...
my $RECAL = 0;

if ($RECAL) {
  # Calculate from recalibrated data in v1:
  $dir = "$ENV{G1K}/v1";
  $base = "recalibrated.fastq.gz";
} else {
  # Calculate from (unrecalibrated) repository data:
  $dir = "$ENV{G1K}/DATA";
  $base = "*f*q*";
}

chdir($dir);

foreach my $d1 (glob("*/NA?????")) {
    foreach my $d2 (glob("$d1/lib*/lane*")) {
	my @files;
	foreach my $f (glob("$d2/$base")) {
	    push @files,$f unless ($f =~ /\.entropy/);
	}
	if (@files != 1) {
	    report(sprintf("%s: found %d files matching %s",$d2,scalar(@files),$base),1);
	    next;
	}
	my $f = $files[0];
	my $fh = openToRead($f);
	my $g = "$f.entropy";
	$g =~ s/gz.entropy$/entropy/;
	if (stamp($g)) {
	    report("$g: creating");
	    my $n=0;
	    my (@c2,@c4);
	    while (<$fh>) {
		$n++;
		if ($n == 2) {
		    addToCounts($_,\@c2);
		} elsif ($n == 4) {
		    addToCounts($_,\@c4);
		    $n = 0;
		}
	    }
	    my $gh = openToWrite($g);
	    foreach my $i (0 .. $#c2) {
		printf($gh "%d %8.6f %8.6f\n",$i,entropy2($c2[$i],1),entropy4($c4[$i],0));
	    }
	    $gh->close();
	    unstamp($g);
	    report();
	}
    }
}

sub addToCounts {
    my ($line,$pC) = @_;
    chomp $line;
    my $len = length($line);
    foreach my $i (0 .. $len-1) {
	$pC->[$i]{substr($line,$i,1)} ++;
    }
}

sub entropy2 {
    my ($pH) = @_;
    my ($n,$h);
    foreach my $nuc (keys %$pH) {
	my $c = $pH->{$nuc};
	if ($nuc ne 'N' && $c > 1) {
	    $n += $c;
	    $h += $c*log($c);
	}
    }
    $h = $n > 0 ? ($n/($n+$pH->{N}))*(log($n) - $h/$n)/log(2) : 0;
    return $h;
}

sub entropy4 {
    my ($pH) = @_;
    my ($n,$h);
    foreach my $nuc (keys %$pH) {
	my $c = $pH->{$nuc};
	if ($c > 1) {
	    $n += $c;
	    $h += $c*log($c);
	}
    }
    $h = $n > 0 ? (log($n) - $h/$n)/log(2) : 0;
    return $h;
}

