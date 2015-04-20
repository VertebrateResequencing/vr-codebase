#!/usr/bin/env perl
use strict;
use warnings;

# Author: Kim Wong kw10@sanger.ac.uk

# This script is part of the mutect_runner package.
# Create a summary file from reformatted mutect files; 
# creates a table with 1 or 0 for each sample/site

my $file = shift @ARGV;

if (!$file) {
  my $help = <<END;

  This script is part of the mutect_runner package.
  This is a tool to create a summary table from multiple
  MuTect tab files created by 'reformatVCF.pl'.
  Input is a list of all tab files followed by a file with 
  sites to be included.

  Usage: checksamples.pl mutect_tabs.list all_sites.txt

  Output is written to STDOUT.

  Format of mutect_tabs.list file (files end in -Mutect.txt)
    1_cell_line-vs-1_spleen-MuTect.txt
    1_tumour-vs-1_spleen-MuTect.txt

  Format of all_sites.txt:

    # create all_sites.txt file (if tables have DP4T column):
    cut -f 1-5,12- 1_cell_line-vs-1_spleen-MuTect.txt | grep ^#CHR > all_sites.txt
    cat < mutect_tabs.list | grep -v ^# | cut -f 1-5,12- | sort -uV >> all_sites.txt

  Example usage:

   checksamples.pl mutect_tabs.list all_sites.txt > mutect_all_summary.txt

   # create a summary file with one sample per line:
   summary_bysample.pl mutect_all_summary.txt > mutect_all_summary.bysample.txt

END

	die $help;
}

elsif (! -e $file) {
	die "No such file $file\n";
}

my @files = `cat $file`;
chomp @files;
my @samples;
foreach my $f (@files) {
	my $s = `basename $f "-MuTect.txt"`;
	chomp $s;
	push @samples, $s;
	print STDERR "$s\n";
}

##print join ("\t",@samples)."\n";
##exit;

my $dp4check = `grep ^# $files[0] | grep DP4T`;
print STDERR "DP found\n" if $dp4check;

# input list of sites
my %inlines;

foreach my $in (@files) {
	my @lines = $dp4check ? `cut -f 1,2,4,5,13,14,16,18 $in | sed 's/\t/_/g'` :  `cut -f 2,4,5,12,13,15,17 $in | sed 's/\t/_/g'`; 
	for (@lines) {
		chomp;
		die "Error: already seen $_\n" if $inlines{$in}{$_};
		$inlines{$in}{$_}=1;
	}
}

while (<>) {
	my $line = $_;
	chomp $line;
	if (/^#/) {
		print "$line\t".join ("\t",@samples)."\n";
		next;
	}
	my @c = split "\t", $line;
	my $match = join ("_", $c[0],$c[1],$c[3],$c[4],$c[6],$c[7],$c[9],$c[11]);
	my @tally;
	foreach my $in (@files) {
		if ($inlines{$in}{$match}){
			push @tally, '1';
		}
		else {
			push @tally, '0';
		}
	}
	print "$line\t", join ("\t", @tally)."\n";
}
