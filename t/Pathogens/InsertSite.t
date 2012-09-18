#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;


BEGIN { unshift(@INC, './modules') }
BEGIN {

    use Test::Most;
    use_ok('Pathogens::RNASeq::InsertSite');
}

ok my $insert_site_plots_from_bam = Pathogens::RNASeq::InsertSite->new(
   filename             => 't/data/small_multi_sequence.bam',
   output_base_filename => 't/data/insert_site',
  );
ok $insert_site_plots_from_bam->create_plots();

# parse output files and check they are okay
ok is_input_string_found_on_given_line("0 0", 1,    't/data/insert_site.FN543502.insert_site_plot.gz'), 'check main sequence insert_site values first value';
ok is_input_string_found_on_given_line("0 2", 7899, 't/data/insert_site.FN543502.insert_site_plot.gz'), 'check main sequence insert_site value before site';
ok is_input_string_found_on_given_line("0 14", 7915,'t/data/insert_site.FN543502.insert_site_plot.gz'), 'check main sequence insert_site values for reverse reads only';
ok is_input_string_found_on_given_line("1 0", 7920, 't/data/insert_site.FN543502.insert_site_plot.gz'), 'check main sequence insert_site values after site';
ok is_input_string_found_on_given_line("0 0", 249,  't/data/insert_site.FN543502.insert_site_plot.gz'), 'various values';
ok is_input_string_found_on_given_line("1 0", 345,  't/data/insert_site.FN543502.insert_site_plot.gz'), 'various values';
ok is_input_string_found_on_given_line("3 0", 354,  't/data/insert_site.FN543502.insert_site_plot.gz'), 'various values';
ok is_input_string_found_on_given_line("1 0", 366,  't/data/insert_site.FN543502.insert_site_plot.gz'), 'various values';
ok is_input_string_found_on_given_line("0 0", 513,  't/data/insert_site.FN543502.insert_site_plot.gz'), 'various values';

ok is_input_string_found_on_given_line("0 0", 1,    't/data/insert_site.pCROD1.insert_site_plot.gz'), 'check empty plasmid insert_site values first value';
ok is_input_string_found_on_given_line("0 0", 59,   't/data/insert_site.pCROD1.insert_site_plot.gz'), 'check empty plasmid insert_site values last value';
ok is_input_string_found_on_given_line("0 0", 1,    't/data/insert_site.pCROD2.insert_site_plot.gz'), 'check plasmid with 1 read insert_site values first value';
ok is_input_string_found_on_given_line("0 1", 143,  't/data/insert_site.pCROD2.insert_site_plot.gz'), 'check plasmid with 1 read insert_site values first base of read';
ok is_input_string_found_on_given_line("0 0", 144,  't/data/insert_site.pCROD2.insert_site_plot.gz'), 'check plasmid with 1 read insert_site values after last base of read';
ok is_input_string_found_on_given_line("0 0", 1000, 't/data/insert_site.pCROD2.insert_site_plot.gz'), 'check plasmid with 1 read insert_site values last value';
ok is_input_string_found_on_given_line("0 0", 1,    't/data/insert_site.pCROD3.insert_site_plot.gz'), 'check another empty plasmid insert_site values first value';
ok is_input_string_found_on_given_line("0 0", 100,  't/data/insert_site.pCROD3.insert_site_plot.gz'), 'check another empty plasmid insert_site values last value';

unlink("t/data/insert_site.FN543502.insert_site_plot.gz");
unlink("t/data/insert_site.pCROD1.insert_site_plot.gz");
unlink("t/data/insert_site.pCROD2.insert_site_plot.gz");
unlink("t/data/insert_site.pCROD3.insert_site_plot.gz");

done_testing();

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
    print "Expected: ".$expected_string ."\t Got: ". $line."\n";
  }
  return 0;
}