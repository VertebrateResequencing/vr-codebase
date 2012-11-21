#!/usr/bin/env perl
use strict;
use warnings;

BEGIN { unshift(@INC, './modules') }
BEGIN {
    use Test::Most;
    use_ok('Pathogens::RNASeq::Read');
}

ok my $alignment_slice = Pathogens::RNASeq::Read->new(
  alignment_line => 'IL25_4928:3:53:7118:13952#4	163	FN543502	66737	60	54M	=	66903	220	GGGGGGCGTTTTCCGGCGATTCTTTACTGTACATATCCAGTTGACCGTTCGGGA	BBBBBBBBBBBBBBBB@B9=B@BBBF@@@@@@@@@@@@@@@@@@?>@?@?@@B?	XT:A:U	NM:i:1	SM:i:37	AM:i:37	X0:i:1	X1:i:0	XM:i:1	XO:i:0	XG:i:0	MD:Z:0C53	RG:Z:4928_3#4',
  gene_strand => 1,
  exons => [[66337,66937],[4,5]]
  ), 'initialise';
  
is $alignment_slice->_read_position, 66737, 'read position';
is $alignment_slice->_read_length, 54, 'read length';
is $alignment_slice->read_strand, 1, 'read strand';

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


# Passes bitwise filter
%filters = (bitwise_flag => 2);
ok my $alignment_slice2 = Pathogens::RNASeq::Read->new(
  alignment_line => 'IL25_4928:3:53:7118:13952#4	3	FN543502	66737	60	54M	=	66903	220	GGGGGGCGTTTTCCGGCGATTCTTTACTGTACATATCCAGTTGACCGTTCGGGA	BBBBBBBBBBBBBBBB@B9=B@BBBF@@@@@@@@@@@@@@@@@@?>@?@?@@B?	XT:A:U	NM:i:1	SM:i:37	AM:i:37	X0:i:1	X1:i:0	XM:i:1	XO:i:0	XG:i:0	MD:Z:0C53	RG:Z:4928_3#4',
  gene_strand => 1,
  exons => [[66337,66937],[4,5]],
  filters => \%filters
  ), 'initialise';
ok  my %mapped_reads2 = %{$alignment_slice2->mapped_reads}, 'build the mapped reads Passes bitwise filter';
is  $mapped_reads2{sense}, 1, 'identified mapped reads Passes bitwise filter';
is  $mapped_reads2{antisense}, 0, 'identified antisense read Passes bitwise filter';

# Doesnt pass bitwise filter
%filters = (bitwise_flag => 2);
ok my $alignment_slice3 = Pathogens::RNASeq::Read->new(
  alignment_line => 'IL25_4928:3:53:7118:13952#4	13	FN543502	66737	60	54M	=	66903	220	GGGGGGCGTTTTCCGGCGATTCTTTACTGTACATATCCAGTTGACCGTTCGGGA	BBBBBBBBBBBBBBBB@B9=B@BBBF@@@@@@@@@@@@@@@@@@?>@?@?@@B?	XT:A:U	NM:i:1	SM:i:37	AM:i:37	X0:i:1	X1:i:0	XM:i:1	XO:i:0	XG:i:0	MD:Z:0C53	RG:Z:4928_3#4',
  gene_strand => 1,
  exons => [[66337,66937],[4,5]],
  filters => \%filters
  ), 'initialise';
ok  my %mapped_reads3 = %{$alignment_slice3->mapped_reads}, 'build the mapped reads Doesnt pass bitwise filter';
is  $mapped_reads3{sense}, 0, 'identified mapped reads Doesnt pass bitwise filter';
is  $mapped_reads3{antisense}, 0, 'identified antisense read Doesnt pass bitwise filter';

# complex filter pass
%filters = (bitwise_flag => 5);
ok my $alignment_slice4 = Pathogens::RNASeq::Read->new(
  alignment_line => 'IL25_4928:3:53:7118:13952#4	3	FN543502	66737	60	54M	=	66903	220	GGGGGGCGTTTTCCGGCGATTCTTTACTGTACATATCCAGTTGACCGTTCGGGA	BBBBBBBBBBBBBBBB@B9=B@BBBF@@@@@@@@@@@@@@@@@@?>@?@?@@B?	XT:A:U	NM:i:1	SM:i:37	AM:i:37	X0:i:1	X1:i:0	XM:i:1	XO:i:0	XG:i:0	MD:Z:0C53	RG:Z:4928_3#4',
  gene_strand => 1,
  exons => [[66337,66937],[4,5]],
  filters => \%filters
  ), 'initialise';
ok  my %mapped_reads4 = %{$alignment_slice4->mapped_reads}, 'build the mapped reads complex filter pass';
is  $mapped_reads4{sense}, 1, 'identified mapped reads complex filter pass';
is  $mapped_reads4{antisense}, 0, 'identified antisense read complex filter pass';


# complex filter fail
%filters = (bitwise_flag => 10);
ok my $alignment_slice5 = Pathogens::RNASeq::Read->new(
  alignment_line => 'IL25_4928:3:53:7118:13952#4	21	FN543502	66737	60	54M	=	66903	220	GGGGGGCGTTTTCCGGCGATTCTTTACTGTACATATCCAGTTGACCGTTCGGGA	BBBBBBBBBBBBBBBB@B9=B@BBBF@@@@@@@@@@@@@@@@@@?>@?@?@@B?	XT:A:U	NM:i:1	SM:i:37	AM:i:37	X0:i:1	X1:i:0	XM:i:1	XO:i:0	XG:i:0	MD:Z:0C53	RG:Z:4928_3#4',
  gene_strand => 1,
  exons => [[66337,66937],[4,5]],
  filters => \%filters
  ), 'initialise';
ok  my %mapped_reads5 = %{$alignment_slice5->mapped_reads}, 'build the mapped reads complex filter fail';
is  $mapped_reads5{sense}, 0, 'identified mapped reads complex filter fail';
is  $mapped_reads5{antisense}, 0, 'identified antisense read complex filter fail';


# unmapped reads fail
ok my $alignment_slice6 = Pathogens::RNASeq::Read->new(
  alignment_line => 'IL25_4928:3:53:7118:13952#4	7	FN543502	66737	60	54M	=	66903	220	GGGGGGCGTTTTCCGGCGATTCTTTACTGTACATATCCAGTTGACCGTTCGGGA	BBBBBBBBBBBBBBBB@B9=B@BBBF@@@@@@@@@@@@@@@@@@?>@?@?@@B?	XT:A:U	NM:i:1	SM:i:37	AM:i:37	X0:i:1	X1:i:0	XM:i:1	XO:i:0	XG:i:0	MD:Z:0C53	RG:Z:4928_3#4',
  gene_strand => 1,
  exons => [[66337,66937],[4,5]]
  ), 'initialise';
ok  my %mapped_reads6 = %{$alignment_slice6->mapped_reads}, 'build the mapped reads does not pass unmapped read';
is  $mapped_reads6{sense}, 0, 'identified mapped reads does not pass unmapped read';
is  $mapped_reads6{antisense}, 0, 'identified antisense read does not pass unmapped read';
done_testing();