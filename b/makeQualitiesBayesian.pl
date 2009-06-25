#!/usr/bin/env perl

use File::Basename;

use lib dirname($0);

use Utility;
use Getopt::Std;
use strict;

my %OPT = qw(m 40); # max of uniform prior on qualities; 40 seems a good default.
getopts('m:',\%OPT);

autoflush STDOUT 1;

my $QMAX = $OPT{m};

makeQualitiesBayesian();

sub makeQualitiesBayesian {
    my @qm;
    while (<>) {
	chomp;
	# Read, Position, UncalibratedQuality, CalibratedQuality, NBases, NErrors
	# "Read" is 0 for single end, 1 for paired end read 1, 2 for paired end read 2.
	my @a = split;
	# Check for conditions that occasionally (used to) arise from splitting error in mapcheck file.
	# If OK, re-estimate the quality, i.e. derive a new CalibratedQuality from NBases and NErrors.
	unless (/\./ || $a[5]>$a[4]) {
	    $qm[$a[0]][$a[1]][$a[2]] = [estimateQuality(@a[4,5]),@a[4,5]];
	}
    }
    foreach my $r (0 .. $#qm) {
	my $pQMR = $qm[$r];
	next unless ($pQMR);
	# Find min and max positions for which qualities are estimated for this Read number
	my $pMax = $#$pQMR;
	my $pMin = 0;
	while (!$pQMR->[$pMin]) {
	    $pMin ++;
	}
	foreach my $p ($pMin .. $pMax) {
	    if (defined $pQMR->[$p]) {
		# Use one in from the end for each end (pMin and pMax)
		my $pp = $p == $pMin ? $pMin+1 : $p == $pMax ? $pMax-1 : $p;
		foreach my $q (0 .. $#{$pQMR->[$p]}) {
		    my $pA = $pQMR->[$pp][$q];
		    if (defined $pA) {
			# Print the results, but now $pA->[0] replaces the input CalibratedQuality.
			print join("\t",$r,$p,$q,@$pA) . "\n";
		    }
		}
	    }
	}
    }
}

# Arguments: n bases, of which m are (we assume) errors.
# Prior distribution is uniform over qualities 0 ... $QMAX.
sub estimateQuality {
    my ($n,$m) = @_;
    # If we have no data or bad data, posterior = prior so we return its mean.
    return $QMAX/2 if ($n == 0 || $m > $n);
    # $lch is log(n Choose m).
    my $lch = 0;
    foreach my $i (0 .. $m-1) {
	$lch += log(($n-$i)/($m-$i));
    }
    my ($top,$bot);
    # For each quality in the prior range...
    foreach my $q (0 .. $QMAX) {
	# $p is error probability derived from $q. We
	# set an upper limit of 0.5 (very low qualities are
	# not very meaningful anyway).
	my $p = min(0.5,exp(-log(10)*0.1*$q));
	# Log prob (log likelihood) of m given n and p
	# is log(p^m * (1-p)^(n-m)), times n Choose m...
	my $lpm = $m*log($p)+($n-$m)*log(1-$p);
	my $lp = $lpm+$lch;
	# Posterior prob is exp($lp) (times the prior which is
	# constant and can be ignored)
	my $post = exp($lp);
	# Work towards estimating the posterior mean...
	$top += $post*$q;
	$bot += $post;
    }
    # If we didn't gather any data, return middle of range.
    return $QMAX/2 if ($bot == 0);
    # Find nearest integer to posterior mean $top/$bot.
    my $q = int(0.5+$top/$bot);
    return $q;
}
