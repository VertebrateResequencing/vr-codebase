#!/usr/local/bin/perl

use strict;
use warnings;
no warnings 'uninitialized';
use Getopt::Long;
use Utility;
use File::Basename;

$|++;

my ($recal_dir, $help);

chdir($ENV{G1K});

GetOptions(
    'analysis_dir=s'	=>  \$recal_dir,
    'h|help'    =>  \$help,
    );

(-d $recal_dir && !$help) or die <<USAGE;
    Usage: $0   
                --analysis_dir  <name of analysis dir in \$G1K>
		--help		<this message>

This script generates "cp -p" commands to copy (newly-created) recalibrated.fastq.gz
files, resulting from "g1k.pl cal", to the repository. It doesn't actually do the
copying; you should redirect its output to a file and then (after checking) execute
that.

Assumes recalibrated files are generated under \$G1K/[analysis_dir], and that the
original fastq are symlinked into \$G1K/DATA

e.g. $0 --analysis_dir recal_2008_08 > recal_copy_commands

USAGE


print "cd $ENV{G1K}\n";

foreach my $f (glob("$recal_dir/*/NA*/lib*/lane*/recalibrated.fastq.gz")) {
    if (0 && -s "$f.touch") {
	print "# $f still being made\n";
	next;
    }
    (my $dd = dirname($f)) =~ s/$recal_dir/DATA/;
    my ($pfx,$sfx);
  LOOP:
    foreach my $suf (qw(fastq fastq.gz fq.gz)) {
	foreach my $sl (glob("$dd/*.$suf")) {
	    if (-l $sl) {
		$pfx = readlink($sl);
		if (-s $pfx) {
		    $pfx = substr($pfx,0,length($pfx)-length($suf)-1);
		    $sfx = $suf;
		    last LOOP;
		} else {
		    $pfx = undef;
		}
	    }
	}
    }
    if ($pfx) {
	my $tgt = "$pfx.recal.$sfx";
	$tgt .= ".gz" unless ($tgt =~ /\.gz/);
	my $st = -s $tgt;
	unless (-s $f == $st && !newer($f,$tgt)) {
	    if ($st) {
		print "# Overwriting $tgt\n";
	    }
	    print "cp -p $f $tgt\n";
	}
    }
}

	
