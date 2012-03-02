#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;

BEGIN { unshift(@INC, './modules') }
BEGIN {
    use Test::Most tests => 15;
    use_ok('Pathogens::RNASeq::GFF');
}

ok my $rna_seq_gff = Pathogens::RNASeq::GFF->new(filename => 't/data/Citrobacter_rodentium_ICC168_v1_test.gff'), 'Initialise valid GFF file';
ok $rna_seq_gff->features(), 'build features';


# Should extract feature in a valid manner
is $rna_seq_gff->features->{"continuous_feature_locus_tag"}->gene_start, 166, 'find gene start of CDS for continuous feature';
is $rna_seq_gff->features->{"continuous_feature_locus_tag"}->gene_end,   231, 'find gene end of CDS for continuous feature';
is $rna_seq_gff->features->{"continuous_feature_locus_tag"}->exon_length, 66, 'exon length for continuous feature';
is $rna_seq_gff->features->{"continuous_feature_locus_tag"}->gene_strand,  1, 'strand of continuous feature';
is $rna_seq_gff->features->{"continuous_feature_locus_tag"}->seq_id,  "FN543502", 'seq_id of continuous feature';

is $rna_seq_gff->features->{"non_cds_feature_id"}, undef, 'Dont parse non CDS features';

# discontinuous feature should be linked up
is $rna_seq_gff->features->{"discontinuous_feature_locus_tag"}->gene_start, 313, 'find gene start of CDS for discontinuous feature';
is $rna_seq_gff->features->{"discontinuous_feature_locus_tag"}->gene_end,   7420, 'find gene end of CDS for discontinuous feature';
is $rna_seq_gff->features->{"discontinuous_feature_locus_tag"}->exon_length, 3894, 'exon length for discontinuous feature';

# reverse strand
is $rna_seq_gff->features->{"reverse_strand_locus_tag"}->gene_strand,  -1, 'reverse strand';
is $rna_seq_gff->features->{"no_strand_identifier_locus_tag"}->gene_strand,  0, 'no strand identifier';

my @expected_gene_ids = ('ROD_00021','ROD_00031','ROD_00071','continuous_feature_locus_tag','discontinuous_feature_locus_tag','no_strand_identifier_locus_tag','reverse_strand_locus_tag');

is_deeply $rna_seq_gff->sorted_gene_ids, \@expected_gene_ids, 'sorting of gene ids';