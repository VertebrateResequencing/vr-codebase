#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;


BEGIN { unshift(@INC, './modules') }
BEGIN {

    use Test::Most;
    use_ok('Pathogens::RNASeq::Insertions');
}


ok(my $expression_results = Pathogens::RNASeq::Insertions->new(
  sequence_filename => 't/data/small_multi_sequence.bam',
  annotation_filename => 't/data/Citrobacter_rodentium_ICC168_v1_test.gff',
  output_base_filename => 'output_results'
  ), "initialise the insertion driver class");

ok($expression_results->output_spreadsheet(),'Create spreadsheet');

ok((-e "output_results.insertion.csv"), 'spredsheet file exists');
ok((-e "output_results.corrected.bam"), "Corrected bam file created");
ok((-e "output_results.corrected.bam.bai"), "Corrected bai file created");
ok((-e "output_results.corrected.bam.intergenic.FN543502.tab.gz"), "Corrected tab file for intergenic regions created");

unlink("output_results.insertion.csv");
unlink("output_results.corrected.bam");
unlink("output_results.corrected.bam.bai");
unlink("output_results.corrected.bam.intergenic.FN543502.tab.gz");
