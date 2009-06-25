#!/usr/bin/env perl

# Find all the fastq files (in whatever guise) in repository directiories.
# Collect the output in a file and hand it to updateDataDirectory.pl for
# greater speed in the latter.

use Utility;
use strict;

autoflush STDOUT 1;

foreach my $n (0 .. 99) {
    my $n4 = sprintf("%4d",$n);
    $n4 =~ tr/ /0/;
    my $dir = "/nfs/repository/d$n4";
    #print "# $dir\n";
    foreach my $d (glob("$dir/SLX_*")) {
	foreach my $suf (qw(fastq fastq.gz fq fq.gz)) {
	    foreach my $f (glob("$d/*.$suf")) {
		print "$f\n";
	    }
	}
    }
}
