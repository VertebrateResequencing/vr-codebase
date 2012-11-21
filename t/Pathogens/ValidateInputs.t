#!/usr/bin/env perl
use strict;
use warnings;

BEGIN { unshift(@INC, './modules') }
BEGIN {
    use Test::Most;
    use_ok('Pathogens::RNASeq::ValidateInputs');
}

# everything matches

ok my $valid_combination = Pathogens::RNASeq::ValidateInputs->new(
  sequence_filename => 't/data/1kg_lane.bam',
  annotation_filename => 't/data/1kg_lane_sq_valid.gff'
), 'initialise valid';
is $valid_combination->are_input_files_valid(), 1, 'validation passes for valid input files';


ok my $invalid_combination = Pathogens::RNASeq::ValidateInputs->new(
  sequence_filename => 't/data/1kg_lane.bam',
  annotation_filename => 't/data/1kg_lane_sq_invalid.gff'
), 'initialise valid';
is $invalid_combination->are_input_files_valid(), 0, 'validation fails for invalid input files';
is $invalid_combination->_sequence_names_match(), 0, 'sequence names dont match';
is $invalid_combination->_lengths_match(), 0, 'lengths dont match';
done_testing();


