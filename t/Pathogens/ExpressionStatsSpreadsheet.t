#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;


BEGIN { unshift(@INC, './modules') }
BEGIN {

    use Test::Most tests => 8;
    use_ok('Pathogens::RNASeq::ExpressionStatsSpreadsheet');
}
ok my $expression_results = Pathogens::RNASeq::ExpressionStatsSpreadsheet->new(
  output_filename => 'my_result_file.csv'), 'initialise';

my %result_1 = (
  seq_id  => "some_name",
  gene_id => "abc123",
  mapped_reads_sense     => 2000, 
  mapped_reads_antisense => 10,
  rpkm_sense             => 15.3245,
  rpkm_antisense         => 1.34324,
  );
  
my %result_2 = (
  seq_id  => "some_name",
  gene_id => "efg456",
  mapped_reads_sense     => 10, 
  mapped_reads_antisense => 200,
  rpkm_sense             => 0,
  rpkm_antisense         => 780.34242543543,
  );
  
ok $expression_results->add_result(\%result_1), 'add first result set';
ok $expression_results->add_result(\%result_2), 'add second result set';
ok $expression_results->build_and_close(), 'build the csv file and close';

open(IN, 'my_result_file.csv');
my $header = <IN>;
my $output_result_1 = <IN>;
my $output_result_2 = <IN>;
$header =~ s/[\r\n]//g;
$output_result_1 =~ s/[\r\n]//g;
$output_result_2 =~ s/[\r\n]//g;

is $header, '"Seq ID",GeneID,"Reads Mapping",RPKM,"Antisense Reads Mapping","Antisense RPKM"', 'header okay';
is $output_result_1, 'some_name,abc123,2000,15.3245,10,1.34324', 'result set 1';
is $output_result_2, 'some_name,efg456,10,0,200,780.34242543543', 'result set 2';

unlink('my_result_file.csv');
