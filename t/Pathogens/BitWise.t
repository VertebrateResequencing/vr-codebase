#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;

BEGIN { unshift(@INC, './modules') }
BEGIN {
    use Test::Most tests => 11;
    use_ok('Pathogens::RNASeq::BitWise');
}

# strand specific protocol
ok my $bitwise = Pathogens::RNASeq::BitWise->new(
  filename => 't/data/rna_seq_bitwise_flags_set.bam',
  output_filename => 't/data/my_file.bam',
  protocol => 'StrandSpecificProtocol'
  ),'initialise StrandSpecificProtocol';
ok $bitwise->update_bitwise_flags(),'update bitwise flags StrandSpecificProtocol';

open(IN, '-|', 'samtools view t/data/my_file.bam | awk \'{print $2;}\'');
my $read_1 = <IN>;
my $read_2 = <IN>;
my $duplicate_read = <IN>;
chomp($read_1);
chomp($read_2);
chomp($duplicate_read);
is $read_1, 163, 'change to forward strand StrandSpecificProtocol';
is $read_2, 179, 'change to reverse strand StrandSpecificProtocol';
is $duplicate_read, 117, 'unmark duplicates StrandSpecificProtocol' ;


# Standard protocol
ok $bitwise = Pathogens::RNASeq::BitWise->new(
  filename => 't/data/rna_seq_bitwise_flags_set.bam',
  output_filename => 't/data/my_file.bam'
  ), 'initialise Standard Protocol';
ok $bitwise->update_bitwise_flags(), 'update bitwise flags Standard Protocol';

open(IN, '-|', 'samtools view t/data/my_file.bam | awk \'{print $2;}\'');
$read_1 = <IN>;
$read_2 = <IN>;
$duplicate_read = <IN>;
chomp($read_1);
chomp($read_2);
chomp($duplicate_read);
is $read_1, 179, 'change to forward strand Standard Protocol';
is $read_2, 163, 'change to reverse strand Standard Protocol';
is $duplicate_read, 117, 'unmark duplicates Standard Protocol' ;

unlink('t/data/my_file.bam');