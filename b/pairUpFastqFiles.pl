#!/usr/bin/env perl

use Utility;
use File::Basename;
use Getopt::Long;
use strict;

my ($source, $target, $help);

GetOptions(
    'source=s'	=>  \$source,
    'target=s'	=>  \$target,
    'h|help'    =>  \$help,
    );

(-d $source && -d $target && !$help) or die <<USAGE;
    Usage: $0   
                --source    <source directory of [name]_1 and [name]_2 files
                --target    <directory to write the joined files>
		--help <this message>

    e.g. $0 --source /DATA/BGI/2008_08/bi_separate \
	    --target /DATA/BGI/2008_08

    Run this on directories containing fastq files (e.g. from Broad) which are
    for paired end data but where the two ends are in two files foo_1.fastq.gz
    and foo_2.fastq.gz.

    This script creates foo.fastq.gz in the target directory. The results
    should then be copied into the repository for use by
    updateDataDirectory.pl -A.

USAGE

pairUpFastqFiles($source, $target);

sub pairUpFastqFiles {
    my ($srcD,$tgtD) = @_;
    foreach my $f1 (glob("$srcD/*_1.fastq.gz")) {
	(my $f2 = $f1) =~ s/_1.fastq/_2.fastq/;
	die("Not found: $f2") unless (-s $f2);
	my $g = join('/',$tgtD,basename($f1));
	$g =~ s/_1.fastq/.fastq/;
	if (stamp($g)) {
	    report("$g: creating");
	    my $fh1 = openToRead($f1);
	    my $fh2 = openToRead($f2);
	    my $gh = openToWrite($g);
	  LOOP:
	    while (1) {
		my (@quad1,@quad2);
		foreach my $i (0 .. 3) {
		    my $line1 = <$fh1>;
		    last LOOP unless ($line1);
		    push @quad1,$line1;
		    my $line2 = <$fh2>;
		    push @quad2,$line2;
		}
		unless (@quad1 == 4 && @quad2 == 4) {
		    print "Mismatch quad:\n@quad1\n@quad2\n";
		    die();
		}
		foreach my $i (0 .. 3) {
		    if ($i == 0 || $i == 2) {
			my $id1 = $quad1[$i];
			$id1 =~ s/ .*//;    # only compare headers up to first space
			my $id2 = $quad2[$i];
			$id2 =~ s/ .*//;

			if ($id1 eq $id2) {
			    print $gh $quad1[$i] or die "Can't print: $!\n";
			} else {
			    die("Mismatch: $quad1[$i]$quad2[$i]");
			}
		    } else {
			chomp (my $seq1 = $quad1[$i]);
			chomp (my $seq2 = $quad2[$i]);
			# We truncate both reads to the length of the shorter one to avoid confusion.
			# The extra data is unlikely to be usable.
			my $len = min(length($seq1),length($seq2));
			printf($gh "%s%s\n",substr($seq1,0,$len),substr($seq2,0,$len)) or die "Can't print: $!\n";
		    }
		}
	    }
	    report();
	    unstamp($g);
	}
    }
}

# Edit this as required:
#pairUpFastqFiles("/lustre/thougen/thougen/EXTERNAL_DATA/BI/2008_08/bi_separate","/lustre/thougen/thougen/EXTERNAL_DATA/BI/2008_08");
#pairUpFastqFiles("/lustre/thougen/thougen/EXTERNAL_DATA/WUGSC/2008_08/bi_separate","/lustre/thougen/thougen/EXTERNAL_DATA/WUGSC/2008_08");
#pairUpFastqFiles("$ENV{G1K}/EXTERNAL_DATA/BI/test/bi_separate","$ENV{G1K}/EXTERNAL_DATA/BI/test");
#pairUpFastqFiles('/nfs/repository/d0037/SLX_1000_EXT_1/bi_new',"$ENV{G1K}/bi_tmp1");

#pairUpFastqFiles('/nfs/repository/d0037/SLX_1000_EXT_1/bi_separate',"$ENV{G1K}/bi_tmp1");
#pairUpFastqFiles('/nfs/repository/d0039/SLX_1000_EXT_4/bi_separate',"$ENV{G1K}/bi_tmp2");

