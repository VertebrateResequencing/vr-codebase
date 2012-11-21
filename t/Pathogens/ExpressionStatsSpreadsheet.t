#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;


BEGIN { unshift(@INC, './modules') }
BEGIN {

    use Test::Most;
    use_ok('Pathogens::RNASeq::ExpressionStatsSpreadsheet');
}
ok my $expression_results = Pathogens::RNASeq::ExpressionStatsSpreadsheet->new(
  output_filename => 'my_result_file.csv',
  protocol => 'StrandSpecificProtocol'
  ), 'initialise';

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

is $header, '"Seq ID",GeneID,"Antisense Reads Mapping","Antisense RPKM","Reads Mapping",RPKM', 'header okay';
is $output_result_1, 'some_name,abc123,2000,15.3245,10,1.34324', 'result set 1';
is $output_result_2, 'some_name,efg456,10,0,200,780.34242543543', 'result set 2';
close(IN);
unlink('my_result_file.csv');


#####################
## Standard protocol
#####################

ok my $expression_results_standard = Pathogens::RNASeq::ExpressionStatsSpreadsheet->new(
  output_filename => 'my_result_file_standard.csv',
  protocol => 'StandardProtocol'
  ), 'initialise';

ok $expression_results_standard->add_result(\%result_1), 'add first result set';
ok $expression_results_standard->add_result(\%result_2), 'add second result set';
ok $expression_results_standard->build_and_close(), 'build the csv file and close';

open(IN_STANDARD, 'my_result_file_standard.csv') or die "Couldnt open input file";
my $header_standard = <IN_STANDARD>;
my $output_result_1_standard = <IN_STANDARD>;
my $output_result_2_standard = <IN_STANDARD>;
$header_standard =~ s/[\r\n]//g;
$output_result_1_standard =~ s/[\r\n]//g;
$output_result_2_standard =~ s/[\r\n]//g;

is $header_standard, '"Seq ID",GeneID,"Reads Mapping",RPKM,"Antisense Reads Mapping","Antisense RPKM"', 'header okay';
is $output_result_1_standard, 'some_name,abc123,2000,15.3245,10,1.34324', 'result set 1';
is $output_result_2_standard, 'some_name,efg456,10,0,200,780.34242543543', 'result set 2';
close(IN_STANDARD);
unlink('my_result_file_standard.csv');


#####################
## Another new protocol
#####################

ok my $expression_results_tradis = Pathogens::RNASeq::ExpressionStatsSpreadsheet->new(
  output_filename => 'my_result_file_tradis.csv',
  protocol => 'TradisProtocol'
  ), 'initialise';

ok $expression_results_tradis->add_result(\%result_1), 'add first result set';
ok $expression_results_tradis->add_result(\%result_2), 'add second result set';
ok $expression_results_tradis->build_and_close(), 'build the csv file and close';

open(IN_TRADIS, 'my_result_file_tradis.csv') or die "Couldnt open input file";
my $header_tradis = <IN_TRADIS>;
my $output_result_1_tradis = <IN_TRADIS>;
my $output_result_2_tradis = <IN_TRADIS>;
$header_tradis =~ s/[\r\n]//g;
$output_result_1_tradis =~ s/[\r\n]//g;
$output_result_2_tradis =~ s/[\r\n]//g;

is $header_tradis,  '"Seq ID",GeneID,"Reads Mapping",RPKM,"Antisense Reads Mapping","Antisense RPKM"', 'header okay';
is $output_result_1_tradis, 'some_name,abc123,2000,15.3245,10,1.34324', 'result set 1';
is $output_result_2_tradis, 'some_name,efg456,10,0,200,780.34242543543', 'result set 2';
close(IN_TRADIS);
unlink('my_result_file_tradis.csv');

done_testing();
