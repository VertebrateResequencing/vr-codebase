#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;

BEGIN { unshift(@INC, './modules') }
BEGIN {
    use Test::Most;
    use_ok('Pathogens::RNASeq::GFF');
}

ok my $rna_seq_gff = Pathogens::RNASeq::GFF->new(filename => 't/data/Citrobacter_rodentium_ICC168_v1_test.gff'), 'Initialise valid GFF file';
ok $rna_seq_gff->features(), 'build features';


# Should extract feature in a valid manner
is $rna_seq_gff->features->{"continuous_feature_id"}->gene_start, 166, 'find gene start of CDS for continuous feature';
is $rna_seq_gff->features->{"continuous_feature_id"}->gene_end,   231, 'find gene end of CDS for continuous feature';
is $rna_seq_gff->features->{"continuous_feature_id"}->exon_length, 66, 'exon length for continuous feature';
is $rna_seq_gff->features->{"continuous_feature_id"}->gene_strand,  1, 'strand of continuous feature';
is $rna_seq_gff->features->{"continuous_feature_id"}->seq_id,  "FN543502", 'seq_id of continuous feature';
is $rna_seq_gff->features->{"continuous_feature_id"}->locus_tag,  "continuous_feature_locus_tag", 'continuous feature locus tag';


is $rna_seq_gff->features->{"non_cds_feature_id"}, undef, 'Dont parse non CDS features';

# discontinuous feature should be linked up
is $rna_seq_gff->features->{"discontinuous_feature_id"}->gene_start, 313, 'find gene start of CDS for discontinuous feature';
is $rna_seq_gff->features->{"discontinuous_feature_id"}->gene_end,   7420, 'find gene end of CDS for discontinuous feature';
is $rna_seq_gff->features->{"discontinuous_feature_id"}->exon_length, 3892, 'exon length for discontinuous feature';
is_deeply($rna_seq_gff->features->{"discontinuous_feature_id"}->exons, [[313, 2775],[5990, 7420]]), 'discontinous exons match';

# reverse strand
is $rna_seq_gff->features->{"reverse_strand_id"}->gene_strand,  -1, 'reverse strand';
is $rna_seq_gff->features->{"no_strand_identifier_id"}->gene_strand,  0, 'no strand identifier';

my @expected_gene_ids = ('FN543502.13','FN543502.18','ROD_00071','continuous_feature_id','discontinuous_feature_id','no_strand_identifier_id','reverse_strand_id');

is_deeply $rna_seq_gff->sorted_gene_ids, \@expected_gene_ids, 'sorting of gene ids';


# Discontinous features
ok( $rna_seq_gff = Pathogens::RNASeq::GFF->new(filename => 't/data/mm10.gff'), 'Discontinuous feature with overall CDS spanning' );
ok(my $feature = $rna_seq_gff->features()->{"10.10"}, 'get discontinous feature');
my @expected_exons = ([10206,10296],
[15002,15271],
[15608,15874],
[16453,16559]);
is_deeply($feature->exons, \@expected_exons , 'exons are as expected');
is($feature->exon_length,731 ,'overall feature length doesnt include introns');

done_testing();
