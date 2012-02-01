#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;


BEGIN { unshift(@INC, './modules') }
BEGIN {

    use Test::Most tests => 12;
    use_ok('Pathogens::RNASeq::CoveragePlot');
}

ok my $coverage_plots_from_bam = Pathogens::RNASeq::CoveragePlot->new(
   filename => 't/data/small_multi_sequence.bam',
   output_base_filename => 't/data/coverage'
  );
ok $coverage_plots_from_bam->create_plots();

# parse output files and check they are okay
# delete temp files
