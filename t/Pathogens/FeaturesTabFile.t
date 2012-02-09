#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;

BEGIN { unshift(@INC, './modules') }
BEGIN {
    use Test::Most tests =>  7;
    use_ok('Pathogens::RNASeq::FeaturesTabFile');
}
use Pathogens::RNASeq::GFF;
use Pathogens::RNASeq::IntergenicRegions;

my $rna_seq_gff = Pathogens::RNASeq::GFF->new(filename => 't/data/Citrobacter_rodentium_ICC168_v1_test.gff');

my $intergenic_regions = Pathogens::RNASeq::IntergenicRegions->new(
  features => $rna_seq_gff->features(),
  window_margin => 50,
  sequence_lengths => $rna_seq_gff->sequence_lengths
  );
my $features = $intergenic_regions->intergenic_features;
my @sequence_names = keys %{$rna_seq_gff->sequence_lengths};


ok my $tab_file_results = Pathogens::RNASeq::FeaturesTabFile->new(
  output_filename => 't/data/intergenic',
  features => $intergenic_regions->intergenic_features,
  sequence_names => \@sequence_names
), 'initialise tab file';
ok $tab_file_results->create_files , 'create tab file';

# parse output files and check they are okay
ok is_input_string_found_on_given_line('FT   misc_feature    1..115',                            1, 't/data/intergenic.FN543502.tab.gz'), 'first feature coords';
ok is_input_string_found_on_given_line('FT                   /colour=2',                         2, 't/data/intergenic.FN543502.tab.gz'), 'first feature colour';
ok is_input_string_found_on_given_line('FT                   /gene="FN543502_intergenic_1_115"', 3, 't/data/intergenic.FN543502.tab.gz'), 'first feature gene name';
ok is_input_string_found_on_given_line('FT                   /seq_id="FN543502"',                4, 't/data/intergenic.FN543502.tab.gz'), 'first feature sequene id';

unlink('t/data/intergenic.FN543502.tab.gz');

sub is_input_string_found_on_given_line
{
  my($expected_string, $line_number, $filename) = @_;
  my $line_counter = 0;
  open(IN, '-|', "gzip -dc ".$filename);
  while(<IN>)
  {
    chomp;
    my $line = $_;
    $line_counter++;
    next unless($line_counter ==  $line_number);
    last if($line_counter >  $line_number);
    return 1 if($expected_string eq $line);
  }
  return 0;
}