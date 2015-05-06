#!/usr/bin/env perl

use strict;
use warnings;

# Author: Kim Wong kw10@sanger.ac.uk
# This script adds double quotes to colums that should be designated
# as text when importing the tab-delimited file into Excel. This is a safegaurd
# against Excel converting gene names into dates. (Trust me, this has happened)

my $help = <<END;

 This script adds double quotes to columns in *-Mutect.txt and MuTect summary files 
 that should be designated as text when importing the tab-delimited file into Excel.
 This is a safegaurd against Excel converting gene names into dates. (Trust me, this has happened)

 Usage: convert_spreadsheet_friendly.pl mutect_tab_file.txt

 Prints to STDOUT.

END

my $f = shift @ARGV or die $help;


my @text = ( 'ID','ID_orig','Cosmic','Cosmic_annot','REF','ALT','MUTECT_FILTER','DP4T','TumourAltBases (Ref,Alt)','NormalAltBases (Ref,Alt)','ENS_GENEID','GENE_SYMBOL','CONSEQUENCE','AA_CHANGE','TRANSCRIPTS','TRANS_BIOTYPE','ENS_ID','Cell_line','Sample');
my %text = map {$_=>1} @text;
my @ind;
open F, $f or die $!;
while (<F>) {
	if (/^##/) {
		print; next;
	}
	elsif (/^#CHROM/) {
		@ind = ();
		my @c = split "\t", $_;
		chomp @c;
		# find out which columns are text
		for (0..$#c) {
			push @ind, $_ if $text{$c[$_]};
		}
		print;
	}
	elsif (/^\s*$/ || /^#CASE/) {
		print;
	}
	else {
		my @c = split "\t", $_;
		chomp @c;
		for (@ind) {
			$c[$_]="\"$c[$_]\"";
		}
		print join ("\t", @c) ."\n";
	}
}
close F;
