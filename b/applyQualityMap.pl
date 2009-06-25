#!/usr/bin/env perl

# Args: input fastq file, quality map file, and output fastq file

# Each line in quality-map file is
#   read position fromQual toQual
# where:
#   read is 1 or 2 (or 0 for single-end)
#   position is position in read (from 1)
#   fromQual and toQuals are integers for qualities

use File::Basename;

use lib dirname($0);

use G1K::G1K;
use Utility;
use strict;

if (@ARGV == 3) {
    applyQualityMap(@ARGV);
} else {
    die("Usage: $0 in.fastq qualmap.txt out.fastq");
}

sub applyQualityMap {
    my ($fqIn,$qMapFile,$fqOut) = @_;

    die("Empty file: $qMapFile") unless (-s $qMapFile);

    my ($len1,$isSingle) = G1K::G1K::getReadLength(dirname($fqIn));

    my @qMap;
    my $mh = openToRead($qMapFile);
    while (<$mh>) {
	chomp;
	my ($r,$p,$m,$n) = split;
	if ($n < 0 || $n > 93) {
	    die("Bad line in $qMapFile: $_");
	}
	$qMap[$r][$p-1][$m] = $n;
    }
    $fqIn = readlink($fqIn) if (-l $fqIn);
    my $fh = openToRead($fqIn);
    report("$fqOut: creating");
    my $gh = openToWrite($fqOut);
    my $i=0;
    while (<$fh>) {
	chomp;
	$i ++;
	if ($i == 4) {
	    $_ = applyQualityMap1(\@qMap,$len1,$isSingle,$_);
	    $i = 0;
	}
	print $gh "$_\n";
    }
    $gh->close();
    my $tf = "$fqOut.touch";
    unlink($tf) if (-f $tf);
    report();
}
	
sub applyQualityMap1 {
    my ($pMap,$len1,$isSingle,$line) = @_;
    my $iLim = length($line);
    my $iHalf = $isSingle ? undef : $len1;
    for (my $i=0; $i<$iLim; $i++) {
	my ($r,$p);
	if ($isSingle) {
	    $r = 0;
	    $p = $i;
	} elsif ($i < $iHalf) {
	    $r = 1;
	    $p = $i;
	} else {
	    $r = 2;
	    $p = $i-$iHalf;
	}
	my $fqc = substr($line,$i,1);
	my $fq = ord($fqc)-33;
	my $pRP = $pMap->[$r][$p];
	my $tq;
	if (defined $pRP) {
	    $tq = $pRP->[$fq];
	    unless (defined $tq) {
		if ($fq > $#$pRP) {
		    $tq = $pRP->[$#$pRP];
		}
		my ($lo,$hi)=($fq,$fq);
		while (!defined $tq) {
		    if (defined $pRP->[$lo]) {
			if (defined $pRP->[$hi]) {
			    $tq = int($pRP->[$lo]+$pRP->[$hi])/2;
			} else {
			    $tq = $pRP->[$lo];
			}
		    } elsif (defined $pRP->[$hi]) {
			$tq = $pRP->[$hi];
		    }
		    if ($lo > 0) {
			$lo--;
			$hi++ if ($hi < $#$pRP);
		    } elsif ($hi < $#$pRP) {
			$hi--;
		    } else {
			last; # give up
		    }
		}
	    }
	}
	unless (defined $tq) {
	    die(sprintf("No mapped quality for %d %d %d",$r,$p+1,$fq))
	}
	my $tqc = chr($tq+33);
	substr($line,$i,1) = $tqc;
    }
    return $line;
}
