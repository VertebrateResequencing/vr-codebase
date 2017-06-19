#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

BEGIN { unshift( @INC, '/software/pathogen/internal/pathdev/vr-codebase/modules/' ) }
use Pathogens::QC::HetSNPCalculator;


my %options = (
    samtools => 'samtools-1.3',
    bcftools => 'bcftools-1.3',
    min_total_depth => 4,
    min_second_depth => 2,
    max_allele_freq => 0.9
);

my $options_ok = GetOptions(\%options,
    'samtools=s',
    'bctfools=s',
    'min_total_depth=i',
    'min_second_depth=i',
    'max_allele_freq=f',
);


unless (@ARGV == 3 and $options_ok) {
    print STDERR "usage: $0 [options] <sorted indexed BAM> <reference fasta> <outfiles prefix>

This script runs the heterozygous SNP calculation, using the same method as the
Pathogens pipeline. The defaults match the pipeline, but you can change them.

Whether or not a position in the reference is counted as heterozygous
is decided as follows.

The depth on each strand is considered independently. First, the
total read depth on each strand must be >= min_total_depth.
Then, on each strand we require:
  number of reads supporting variant >= min_second_depth;
  (number of reads supporting variant) / (total depth) <= max_allele_freq;
If a variant satisfies these requirements on both strands, then it is
counted. If a given position in the genome has at least 2 such variants,
it is counted as heterozygous.

The script options (see above explanation) are:

-min_total_depth
    minimum total read depth on each strand [default: 4]

-min_second_depth
    see explanation above [default: 2]

-max_allele_freq
    see explanation above [default: 0.9]
";

    exit(1);
}


my $het_snp_calc = Pathogens::QC::HetSNPCalculator->new(
    samtools => $options{samtools},
    bcftools => $options{bcftools},
    min_total_depth => $options{min_total_depth},
    min_second_depth => $options{min_second_depth},
    max_allele_freq => $options{max_allele_freq},
    fa_ref => $ARGV[1],
    bam => $ARGV[0],
    outprefix => $ARGV[2],
);

$het_snp_calc->run();
