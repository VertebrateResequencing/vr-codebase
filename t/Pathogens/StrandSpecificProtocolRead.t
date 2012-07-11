#!/usr/bin/env perl
use strict;
use warnings;

BEGIN { unshift(@INC, './modules') }
BEGIN {
    use Test::Most tests => 9;
    use_ok('Pathogens::RNASeq::StrandSpecificProtocol::Read');
}

ok my $alignment_slice = Pathogens::RNASeq::StrandSpecificProtocol::Read->new(
  alignment_line => 'IL25_4928:3:53:13462:11492#2	147	FN543502	1918	29	1S53M	=	1689	-282	ACGCTTAATTCGCCTGGTGAAAGAATATCATCTACTCAATCCGGTGATTGTTGA	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB@BBBB?BBB??BBBBBBBBB?BB	XT:A:M	NM:i:5	SM:i:29	AM:i:29	XM:i:5	XO:i:0	XG:i:0	MD:Z:26C5G2T2C11C2	RG:Z:4928_3#2
	',
  gene_strand => 1,
  exons => [[1818,2918],[4,5]]
  ), 'initialise';

is $alignment_slice->_read_details->{flag}, 131, 'flipped flag';

# check the other bits work as per normal
is $alignment_slice->_read_position, 1918, 'read position';
is $alignment_slice->_read_length, 54, 'read length';
is $alignment_slice->read_strand, 1, 'read strand';
ok  my %mapped_reads = %{$alignment_slice->mapped_reads}, 'build the mapped reads';
is  $mapped_reads{sense}, 1, 'identified mapped reads';
is  $mapped_reads{antisense}, 0, 'identified antisense read';
