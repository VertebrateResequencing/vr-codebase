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
  normalised_pos_insert_sites => 1,
  pos_insert_sites => 2,
  normalised_neg_insert_sites => 3,
  neg_insert_sites => 4,
  normalised_zero_insert_sites => 5,
  zero_insert_sites => 6,
  normalised_total_insert_sites => 7,
  total_insert_sites => 8,
  normalised_pos_insert_site_reads => 9,
  pos_insert_site_reads => 10,
  normalised_neg_insert_site_reads => 11,
  neg_insert_site_reads => 12,
  normalised_zero_insert_site_reads => 13,
  zero_insert_site_reads => 14,
  normalised_total_insert_site_reads => 15,
  total_insert_site_reads => 16,
  );

  
my %result_2 = (
  seq_id  => "some_name",
  gene_id => "efg456",
  normalised_pos_insert_sites => 91,
  pos_insert_sites => 92,
  normalised_neg_insert_sites => 93,
  neg_insert_sites => 94,
  normalised_zero_insert_sites => 95,
  zero_insert_sites => 96,
  normalised_total_insert_sites => 97,
  total_insert_sites => 98,
  normalised_pos_insert_site_reads => 99,
  pos_insert_site_reads => 910,
  normalised_neg_insert_site_reads => 911,
  neg_insert_site_reads => 912,
  normalised_zero_insert_site_reads => 913,
  zero_insert_site_reads => 914,
  normalised_total_insert_site_reads => 915,
  total_insert_site_reads => 916,
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

is $header, '"Seq ID",GeneID,"Normalised pos insert sites","Pos insert sites","Normalised neg insert sites","Neg insert sites","Normalised unknown strand insert sites","Unknown strand insert sites","Normalised total insert sites","Total insert sites","Normalised pos insertions","Pos insertions","Normalised neg insertions","Neg insertions","Normalised unknown strand insertions","Unknown strand insertions","Normalised total insertions","Total insertions"', 'header okay';
is $output_result_1, 'some_name,abc123,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16', 'result set 1';
is $output_result_2, 'some_name,efg456,91,92,93,94,95,96,97,98,99,910,911,912,913,914,915,916', 'result set 2';
close(IN);
unlink('my_result_file.csv');
