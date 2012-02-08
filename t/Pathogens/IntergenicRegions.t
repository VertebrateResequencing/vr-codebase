#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;

BEGIN { unshift(@INC, './modules') }
BEGIN {
    use Test::Most tests =>  13;
    use_ok('Pathogens::RNASeq::IntergenicRegions');
}
use Pathogens::RNASeq::GFF;
use Pathogens::RNASeq::FeaturesTabFile;
my $rna_seq_gff = Pathogens::RNASeq::GFF->new(filename => 't/data/Citrobacter_rodentium_ICC168_v1_test.gff');


# valid intergenic region generation where some features are close by
ok my $intergenic_regions = Pathogens::RNASeq::IntergenicRegions->new(
  features => $rna_seq_gff->features(),
  window_margin => 50,
  sequence_lengths => $rna_seq_gff->sequence_lengths
  ), 'initialise intergenic regions';

ok my $features = $intergenic_regions->intergenic_features, 'create features';
my @expected_keys = ('FN543502_intergenic_1_115','FN543502_intergenic_282_2727','FN543502_intergenic_5048_5096','FN543502_intergenic_7471_7640','FN543502_intergenic_8695_8707','FN543502_intergenic_9396_5346659');
my @actual_keys = sort keys(%{$features});


is_deeply \@actual_keys, \@expected_keys, 'expected keys';
is $features->{"FN543502_intergenic_7471_7640"}->gene_id, 'FN543502_intergenic_7471_7640', 'gene id';
is $features->{"FN543502_intergenic_7471_7640"}->gene_start ,7471, 'gene start';
is $features->{"FN543502_intergenic_7471_7640"}->gene_end ,7640, 'gene end';
is $features->{"FN543502_intergenic_7471_7640"}->gene_strand ,1, 'gene strand';
is $features->{"FN543502_intergenic_7471_7640"}->seq_id ,'FN543502', 'sequence id';
my @expected_exons  = ([7471,7640]);
is_deeply $features->{"FN543502_intergenic_7471_7640"}->exons, \@expected_exons, 'expected exons';


# intergenic regions where there is a tiny window margin
ok $intergenic_regions = Pathogens::RNASeq::IntergenicRegions->new(
  features => $rna_seq_gff->features(),
  window_margin => 0,
  minimum_size => 0,
  sequence_lengths => $rna_seq_gff->sequence_lengths
    ), 'initialise intergenic regions';

ok $features = $intergenic_regions->intergenic_features, 'create features';
@expected_keys = ('FN543502_intergenic_1_165','FN543502_intergenic_232_2777','FN543502_intergenic_4998_5146','FN543502_intergenic_5921_5989','FN543502_intergenic_7421_7690','FN543502_intergenic_8645_8757','FN543502_intergenic_9346_5346659');
@actual_keys = sort keys(%{$features});

is_deeply \@actual_keys, \@expected_keys, 'expected keys';

