#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;


BEGIN { unshift(@INC, './modules') }
BEGIN {

    use Test::Most;
    use_ok('Pathogens::RNASeq::InsertionStatsSpreadsheet');
}
ok my $expression_results = Pathogens::RNASeq::InsertionStatsSpreadsheet->new(
  output_filename => 'my_result_file.csv',
  ), 'initialise';

my %result_1 = (
  seq_id  => "some_name",
  gene_id => "abc123",
  forward_insert_sites     => 2000, 
  normalised_forward_insert_sites => 10,
  reverse_insert_sites             => 15,
  normalised_reverse_insert_sites         => 5,
  total_insert_sites => 2500,
  normalised_total_insert_sites => 300,
  );

  
my %result_2 = (
  seq_id  => "some_name",
  gene_id => "efg456",
  forward_insert_sites     => 9, 
  normalised_forward_insert_sites => 99,
  reverse_insert_sites             => 999,
  normalised_reverse_insert_sites         => 9999,
  total_insert_sites => 99999,
  normalised_total_insert_sites =>999999,
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

is $header, '"Seq ID",GeneID,"Forward insert sites","Normalised forward insert sites","Reverse insert sites","Normalised reverse insert sites","Total insert sites","Normalised total insert sites"', 'header okay';
is $output_result_1, 'some_name,abc123,2000,10,15,5,2500,300', 'result set 1';
is $output_result_2, 'some_name,efg456,9,99,999,9999,99999,999999', 'result set 2';
close(IN);
unlink('my_result_file.csv');
