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

  Assumes the first 14 columns are annotations, sample info in 15+ unless specified.

  Usage: summary_bysample.pl mutect_all_summary.txt [column_number] > mutect_all_summary.bysample.txt 

END


my $file = shift @ARGV or die $help;
my $colnum = shift @ARGV;
$colnum = '15' if !$colnum; # 1-based assumed
$colnum--;

if (! -e $file) {
	die "No such file $file\n";
}
open F, $file or die $!;

my @samples;
my $recurrency;

while (<F>) {
	chomp;
	if (/^##/) {
		print $_."\n"; next;
	}
	elsif (/^#/ ) {
		my @c = split "\t",$_;
		$recurrency=1 if $c[-1] eq 'recurrency';
		@samples = $recurrency ? @c[($colnum)..$#c-1] : @c[($colnum)..$#c];
		chomp @samples;
		print join ("\t",@c[0..($colnum-1)])."\tSample\n";
		next;
	}
	my @line = split "\t", $_;
	pop @line if $recurrency ;
	my @cell = @line[$colnum..$#line];
	foreach my $i (0..$#cell) {
		if ($cell[$i]==1) {
			print join ("\t",@line[0..($colnum-1)])."\t$samples[$i]\n";
		}
	}
}

close F;

