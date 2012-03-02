#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;

BEGIN { unshift(@INC, './modules') }
BEGIN {

    use Test::Most tests => 4;
    use_ok('Pathogens::RNASeq::IntergenicFeature');
}

ok my $intergenicfeature = Pathogens::RNASeq::IntergenicFeature->new(
   gene_start => 100,
   gene_end   => 200,
   seq_id     => "ABC"
  ), 'initialise intergenic feature';
is $intergenicfeature->gene_id(), "ABC_intergenic_100_200", 'build the correct gene ID';
my @expected_exons  = ([100,200]);

is_deeply $intergenicfeature->exons(), \@expected_exons, 'build correct exons';

