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
    xargs cat < mutect_tabs.list | grep -v ^# | cut -f 1-5,12- | sort -uV >> all_sites.txt

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
my %seen;
my $full = 0;
foreach my $f (@files) {
	my $line = `grep ^#CHR $f`;
	my $s = $1 if $line =~/NormalAltBases\s\S+\t(\S+)/;
	push @samples, $s;
	print STDERR "$s\n";
	if ($seen{s}) {
		print STDERR "Sample name not unique $s: using names in filename\n";
		$full=1;
		last;
	}
	$seen{$s}++ if !$seen{$s};
}
if ($full==1) {
	@samples = ();
	foreach my $f (@files) {
		my $s = `basename $f "-MuTect.txt"`;
		chomp $s;
		push @samples, $s;
		print STDERR "$s\n";
	}
}
##print join ("\t",@samples)."\n";
##exit;

my $dp4check = `grep ^# $files[0] | grep DP4T`;
print STDERR "DP found\n" if $dp4check;
my $pick = `grep ^# $files[0] | grep PICK`;
print STDERR "PICK found\n" if $pick;
my $checkheader = `grep "#CHR" $files[0] | sed 's/^#//'`;

my @csqorder = split /\t/, $checkheader;
chomp @csqorder;
# record the index position of each annotation
my %csqindex;
foreach my $i (0..$#csqorder) {
	# add 1 to get the column num instead of the index number
	$csqindex{$csqorder[$i]}=$i+1;
}




# input list of sites
my %inlines;
#CHROM  POS     REF     ALT     GENE_SYMBOL     CONSEQUENCE     PROT_POS        TRANSCRIPTS
my $cols = join (",", $csqindex{CHROM},$csqindex{POS},$csqindex{REF},$csqindex{ALT},$csqindex{GENE_SYMBOL},$csqindex{CONSEQUENCE},$csqindex{PROT_POS},$csqindex{TRANSCRIPTS});

foreach my $in (@files) {
	my @lines = `cut -f $cols $in | sed 's/\t/_/g'`;
	for (@lines) {
		chomp;
		die "Error: already seen $_\n" if $inlines{$in}{$_};
		$inlines{$in}{$_}=1;
	}
}
my %csqindex2;
while (<>) {
	my $line = $_;
	chomp $line;
	if ($line=~/^##/) {
		print "$line\n";
		next;
	}
	elsif ($line=~/^#/) {
		print "$line\t".join ("\t",@samples)."\n";
		$line =~ s/^#//;
		my @csqorder = split /\t/, $line;
		chomp @csqorder;
		# record the index position of each annotation
		foreach my $i (0..$#csqorder) {
			$csqindex2{$csqorder[$i]}=$i+1;
		}
		next;
	}
	my @c = split "\t", $line;
	unshift @c, "NA";
	my $match = join ("_", $c[$csqindex2{CHROM}],$c[$csqindex2{POS}],$c[$csqindex2{REF}],$c[$csqindex2{ALT}],$c[$csqindex2{GENE_SYMBOL}],$c[$csqindex2{CONSEQUENCE}],$c[$csqindex2{PROT_POS}],$c[$csqindex2{TRANSCRIPTS}]);
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
