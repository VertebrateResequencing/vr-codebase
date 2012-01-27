#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;


BEGIN { unshift(@INC, './modules') }
BEGIN {

    use Test::Most tests => 12;
    use_ok('Pathogens::RNASeq::AlignmentSlice');
}
use Pathogens::RNASeq::GFF;

my $rna_seq_gff = Pathogens::RNASeq::GFF->new(filename => 't/data/Citrobacter_rodentium_ICC168_v1_test.gff');
my $feature = $rna_seq_gff->features()->{continuous_feature_id};
$feature->exon_length(50);
$feature->gene_strand(1);
my @exons;
push @exons, [66630,66940];
$feature->exons(\@exons);

ok my $alignment_slice = Pathogens::RNASeq::AlignmentSlice->new(
  filename => 't/data/bam',
  window_margin => 10,
  total_mapped_reads => 10000,
  feature => $feature,
  _input_slice_filename => "t/data/Citrobacter_rodentium_slice"
), 'initialise alignment slice';
is $alignment_slice->_window_start, 156, 'start window';
is $alignment_slice->_window_end, 241, 'end window';
ok $alignment_slice->_slice_file_handle, 'file handle initialises okay';
ok my $rpkm_values = $alignment_slice->rpkm_values;
is $rpkm_values->{rpkm_sense}, 54000, 'rpkm sense';
is $rpkm_values->{rpkm_antisense},0, 'rpkm antisense';
is $rpkm_values->{mapped_reads_sense},27, 'mapped reads sense';
is $rpkm_values->{mapped_reads_antisense},0, 'mapped reads antisense';



# invalid filehandle
ok $alignment_slice = Pathogens::RNASeq::AlignmentSlice->new(
  filename => 't/data/bam',
  total_mapped_reads => 10000,
  window_margin => 10,
  feature => $feature,
  _input_slice_filename => "file_which_doesnt_exist"
), 'initialise invalid alignment slice';
throws_ok  {$alignment_slice->_slice_file_handle} qr/Cant view slice/ , 'invalid file should throw an error';