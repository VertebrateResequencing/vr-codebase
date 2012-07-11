#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;


BEGIN { unshift(@INC, './modules') }
BEGIN {

    use Test::Most tests => 18;
    use_ok('Pathogens::RNASeq::CoveragePlot');
}

ok my $coverage_plots_from_bam = Pathogens::RNASeq::CoveragePlot->new(
   filename             => 't/data/small_multi_sequence.bam',
   output_base_filename => 't/data/coverage',
   mpileup_cmd          => "samtools mpileup"
  );
ok $coverage_plots_from_bam->create_plots();

# parse output files and check they are okay
ok is_input_string_found_on_given_line("0 0", 1,    't/data/coverage.FN543502.coverageplot.gz'), 'check main sequence coverage values first value';
ok is_input_string_found_on_given_line("1 0", 104,  't/data/coverage.FN543502.coverageplot.gz'), 'check main sequence coverage values for forward read only';
ok is_input_string_found_on_given_line("0 4", 548,  't/data/coverage.FN543502.coverageplot.gz'), 'check main sequence coverage values for reverse reads only';
ok is_input_string_found_on_given_line("7 3", 7795, 't/data/coverage.FN543502.coverageplot.gz'), 'check main sequence coverage values for both';
ok is_input_string_found_on_given_line("0 0", 8974, 't/data/coverage.FN543502.coverageplot.gz'), 'check main sequence coverage values last value';

ok is_input_string_found_on_given_line("0 0", 1,    't/data/coverage.pCROD1.coverageplot.gz'), 'check empty plasmid coverage values first value';
ok is_input_string_found_on_given_line("0 0", 59,   't/data/coverage.pCROD1.coverageplot.gz'), 'check empty plasmid coverage values last value';

ok is_input_string_found_on_given_line("0 0", 1,    't/data/coverage.pCROD2.coverageplot.gz'), 'check plasmid with 1 read coverage values first value';
ok is_input_string_found_on_given_line("0 0", 89,   't/data/coverage.pCROD2.coverageplot.gz'), 'check plasmid with 1 read coverage values before first base of read';
ok is_input_string_found_on_given_line("0 1", 90,   't/data/coverage.pCROD2.coverageplot.gz'), 'check plasmid with 1 read coverage values first base of read';
ok is_input_string_found_on_given_line("0 1", 143,  't/data/coverage.pCROD2.coverageplot.gz'), 'check plasmid with 1 read coverage values last base of read';
ok is_input_string_found_on_given_line("0 0", 144,  't/data/coverage.pCROD2.coverageplot.gz'), 'check plasmid with 1 read coverage values after last base of read';
ok is_input_string_found_on_given_line("0 0", 1000, 't/data/coverage.pCROD2.coverageplot.gz'), 'check plasmid with 1 read coverage values last value';

ok is_input_string_found_on_given_line("0 0", 1,   't/data/coverage.pCROD3.coverageplot.gz'), 'check another empty plasmid coverage values first value';
ok is_input_string_found_on_given_line("0 0", 100, 't/data/coverage.pCROD3.coverageplot.gz'), 'check another empty plasmid coverage values last value';

unlink("t/data/coverage.FN543502.coverageplot.gz");
unlink("t/data/coverage.pCROD1.coverageplot.gz");
unlink("t/data/coverage.pCROD2.coverageplot.gz");
unlink("t/data/coverage.pCROD3.coverageplot.gz");

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