#!/usr/bin/env perl
use strict;
use warnings;

# Author: Kim Wong kw10@sanger.ac.uk

# This script is part of the mutect_runner package
# Reformat mutect_summary.txt to create a list of
# sites with samples listed on separate lines
# ie: one sample per line; sites and annotations are duplicated
# chr site1 annotations sample1
# chr site1 annotations sample2
# chr site2 annotations sample4
# chr site3 annotations sample1

my $help = <<END;

  This tool is part of the mutect_runner package.
  Input a mutect summary table and reformat to 
  create a file with one sample per line, eg:

    chr site1 annotations sample1,sample2

  becomes:

    chr site1 annotations sample1
    chr site1 annotations sample2

  Assumes the first 14 columns are annotations, sample info in 15+.

  Usage: summary_bysample.pl mutect_all_summary.txt > mutect_all_summary.bysample.txt 

END


my $file = shift @ARGV or die $help;
if (! -e $file) {
	die "No such file $file\n";
}
open F, $file or die $!;

my @samples;
while (<F>) {
	chomp;
	if (/^#/ ) {
		my @c = split "\t",$_;
		@samples = @c[14..$#c];
		chomp @samples;
		print join ("\t",@c[0..13])."\tSample\n";
		next;
	}
	my @line = split "\t", $_;
	my @cell = @line[14..$#line];
	foreach my $i (0..$#cell) {
		if ($cell[$i]==1) {
			print join ("\t",@line[0..13])."\t$samples[$i]\n";
		}
	}
}

close F;

