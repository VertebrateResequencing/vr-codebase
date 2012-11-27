#!/usr/bin/env perl
use strict;
use warnings;

BEGIN { unshift(@INC, './modules') }
BEGIN {
    use Test::Most;
    use_ok('Pathogens::RNASeq::SequenceFile');
}

# get the total mapped reads
ok my $rna_seq_bam = Pathogens::RNASeq::SequenceFile->new(filename => 't/data/rna_seq_bitwise_flags_set.bam'), 'initialise the rna seq bam file container';
is $rna_seq_bam->total_mapped_reads, 2, 'extracted total mapped reads';
done_testing();