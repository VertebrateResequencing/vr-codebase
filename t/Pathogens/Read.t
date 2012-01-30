#!/usr/bin/env perl
use strict;
use warnings;

BEGIN { unshift(@INC, './modules') }
BEGIN {
    use Test::Most tests => 18;
    use_ok('Pathogens::RNASeq::Read');
}

ok my $alignment_slice = Pathogens::RNASeq::Read->new(
  alignment_line => 'IL25_4928:3:53:7118:13952#4	163	FN543502	66737	60	54M	=	66903	220	GGGGGGCGTTTTCCGGCGATTCTTTACTGTACATATCCAGTTGACCGTTCGGGA	BBBBBBBBBBBBBBBB@B9=B@BBBF@@@@@@@@@@@@@@@@@@?>@?@?@@B?	XT:A:U	NM:i:1	SM:i:37	AM:i:37	X0:i:1	X1:i:0	XM:i:1	XO:i:0	XG:i:0	MD:Z:0C53	RG:Z:4928_3#4',
  gene_strand => 1,
  exons => [[66337,66937],[4,5]]
  ), 'initialise';
  
is $alignment_slice->_read_position, 66737, 'read position';
is $alignment_slice->_read_length, 54, 'read length';
is $alignment_slice->_read_strand, 1, 'read strand';

ok  my %mapped_reads = %{$alignment_slice->mapped_reads}, 'build the mapped reads';
is  $mapped_reads{sense}, 1, 'identified mapped reads';
is  $mapped_reads{antisense}, 0, 'identified antisense read';


# filter out low quality reads 
my %filters = (mapping_quality => 30);
ok $alignment_slice = Pathogens::RNASeq::Read->new(
  alignment_line => 'IL25_4928:3:53:7118:13952#4	163	FN543502	66737	20	54M	=	66903	220	GGGGGGCGTTTTCCGGCGATTCTTTACTGTACATATCCAGTTGACCGTTCGGGA	BBBBBBBBBBBBBBBB@B9=B@BBBF@@@@@@@@@@@@@@@@@@?>@?@?@@B?	XT:A:U	NM:i:1	SM:i:37	AM:i:37	X0:i:1	X1:i:0	XM:i:1	XO:i:0	XG:i:0	MD:Z:0C53	RG:Z:4928_3#4',
  gene_strand => 1,
  exons => [[66337,66937],[4,5]],
  filters => \%filters
  ), 'initialise';
is $alignment_slice->_does_read_pass_filters, 0, 'filter low quality reads';
ok  %mapped_reads = %{$alignment_slice->mapped_reads}, 'build the mapped reads with filter';
is  $mapped_reads{sense},0, 'identified mapped reads with filter';
is  $mapped_reads{antisense}, 0, 'identified antisense read with filter';


# dont filter high quality reads
%filters = (mapping_quality => 30);
ok $alignment_slice = Pathogens::RNASeq::Read->new(
  alignment_line => 'IL25_4928:3:53:7118:13952#4	163	FN543502	66737	60	54M	=	66903	220	GGGGGGCGTTTTCCGGCGATTCTTTACTGTACATATCCAGTTGACCGTTCGGGA	BBBBBBBBBBBBBBBB@B9=B@BBBF@@@@@@@@@@@@@@@@@@?>@?@?@@B?	XT:A:U	NM:i:1	SM:i:37	AM:i:37	X0:i:1	X1:i:0	XM:i:1	XO:i:0	XG:i:0	MD:Z:0C53	RG:Z:4928_3#4',
  gene_strand => 1,
  exons => [[66337,66937],[4,5]],
  filters => \%filters
  ), 'initialise';
is $alignment_slice->_does_read_pass_filters, 1, 'dont filter high quality reads';
ok  %mapped_reads = %{$alignment_slice->mapped_reads}, 'build the mapped reads with high quality read and filter';
is  $mapped_reads{sense},1, 'identified mapped reads with high quality read and filter';
is  $mapped_reads{antisense}, 0, 'identified antisense read with high quality read and filter';
