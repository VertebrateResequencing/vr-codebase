#!/usr/bin/env perl
use strict;
use warnings;
use File::Temp;
use File::Spec;
use File::Compare;
use Data::Dumper;

BEGIN {
    use Test::Most;
    use_ok('Pathogens::QC::HetSNPCalculator');
}

ok my $hsc = Pathogens::QC::HetSNPCalculator->new(
    samtools         => 'samtools-1.3',
    bcftools         => 'bcftools-1.3',
    fa_ref           => 't/data/het_snp_cal_run.ref.fa',
    bam              => 't/data/het_snp_cal_run.mapped_reads.bam',
    outprefix        => 'tmp.test_het_snp_calculator',
    min_total_depth  => 4,
    min_second_depth => 2,
    max_allele_freq  => 0.9,
  ),
  'create instance of object';



my $tmp_vcf = 'tmp.get_snp_calc_test.mpileup.vcf';
unlink $tmp_vcf if -e $tmp_vcf;
$hsc->_run_mpileup($tmp_vcf);
# The vcf file could vary in many ways depending on program version etc, so
# just check the vcf is there and non-empty
is((-e $tmp_vcf and -s $tmp_vcf > 0), 1, '_run_mpileup made a vcf');
unlink $tmp_vcf;


my @got_sum_max = $hsc->_list_sum_and_max([2, 3, 1]);
my @expect = (6, 3);
is_deeply(\@got_sum_max, \@expect, 'test _list_sum_and_max');


my $vcf_line = "chrom1\t42\t.\tT\tC,<*>\t0\t.\tDP=54;ADF=22,7,0;ADR=17,8,0;AD=39,15,0;I16=22,17,7,8,1421,51911,551,20255,1928,95376,675,30457,730,16036,289,6283;\tPL\t214,0,255,255,255,255\n";
@expect = ('chrom1', 1, 1);
my @got = $hsc->_parse_vcf_line(\$vcf_line);
is_deeply(\@got, \@expect, '_parse_vcf_line het line chrom1');

$vcf_line = "chrom2\t42\t.\tT\tC,<*>\t0\t.\tDP=54;ADF=22,7,0;ADR=17,8,0;AD=39,15,0;I16=22,17,7,8,1421,51911,551,20255,1928,95376,675,30457,730,16036,289,6283;\tPL\t214,0,255,255,255,255\n";
@expect = ('chrom2', 1, 1);
@got = $hsc->_parse_vcf_line(\$vcf_line);
is_deeply(\@got, \@expect, '_parse_vcf_line het line chrom2');

$vcf_line = "chrom1\t42\t.\tT\tC,<*>\t0\t.\tDP=54;ADF=22,1,0;ADR=17,10,0;AD=39,15,0;I16=22,17,7,8,1421,51911,551,20255,1928,95376,675,30457,730,16036,289,6283;\tPL\t214,0,255,255,255,255\n";
@expect = ('chrom1', 0, 0);
@got = $hsc->_parse_vcf_line(\$vcf_line);
is_deeply(\@got, \@expect, '_parse_vcf_line not het line because ADF fail');

$vcf_line = "chrom1\t42\t.\tT\tC,<*>\t0\t.\tDP=54;ADF=22,10,0;ADR=17,1,0;AD=39,15,0;I16=22,17,7,8,1421,51911,551,20255,1928,95376,675,30457,730,16036,289,6283;\tPL\t214,0,255,255,255,255\n";
@expect = ('chrom1', 0, 0);
@got = $hsc->_parse_vcf_line(\$vcf_line);
is_deeply(\@got, \@expect, '_parse_vcf_line not het line because ADR fail');

$vcf_line = "chrom1\t42\t.\tT\tC,<*>\t0\t.\tDP=54;ADF=22,1,0;ADR=17,1,0;AD=39,15,0;I16=22,17,7,8,1421,51911,551,20255,1928,95376,675,30457,730,16036,289,6283;\tPL\t214,0,255,255,255,255\n";
@expect = ('chrom1', 0, 0);
@got = $hsc->_parse_vcf_line(\$vcf_line);
is_deeply(\@got, \@expect, '_parse_vcf_line not het line because ADF and ADR fail');

$vcf_line = "chrom1\t42\t.\tT\tC,<*>\t0\t.\tDP=54;ADF=1,22,0;ADR=1,17,0;AD=39,15,0;I16=22,17,7,8,1421,51911,551,20255,1928,95376,675,30457,730,16036,289,6283;\tPL\t214,0,255,255,255,255\n";
@expect = ('chrom1', 1, 0);
@got = $hsc->_parse_vcf_line(\$vcf_line);
is_deeply(\@got, \@expect, '_parse_vcf_line snp but not het');

$vcf_line = "chrom1\t42\t.\tT\tC,<*>\t0\t.\tDP=54;ADF=1,22,0;ADR=1,17,2,0;AD=39,15,0;I16=22,17,7,8,1421,51911,551,20255,1928,95376,675,30457,730,16036,289,6283;\tPL\t214,0,255,255,255,255\n";
throws_ok { $hsc->_parse_vcf_line(\$vcf_line) } 'Pathogens::Exception::VcfParse' , 'Throws when ADF/ADR legnths do not match';

$vcf_line = "chrom1\t42\t.\tT\tC,<*>\t0\t.\tDP=54;ADR=1,2,0;AD=39,15,0;I16=22,17,7,8,1421,51911,551,20255,1928,95376,675,30457,730,16036,289,6283;\tPL\t214,0,255,255,255,255\n";
throws_ok { $hsc->_parse_vcf_line(\$vcf_line) } 'Pathogens::Exception::VcfParse' , 'Throws when ADF missing';

$vcf_line = "chrom1\t42\t.\tT\tC,<*>\t0\t.\tDP=54;ADF=1,2,0;AD=39,15,0;I16=22,17,7,8,1421,51911,551,20255,1928,95376,675,30457,730,16036,289,6283;\tPL\t214,0,255,255,255,255\n";
throws_ok { $hsc->_parse_vcf_line(\$vcf_line) } 'Pathogens::Exception::VcfParse' , 'Throws when ADR missing';


my $vcf_in = 't/data/het_snp_cal_filter_vcf_and_count_snps.in.vcf';
my $expected_vcf = 't/data/het_snp_cal_filter_vcf_and_count_snps.expected.vcf';
my $tmp_vcf_out = 'tmp.test.HetSNPCalculator.filter_vcf_and_count_snps.out.vcf';
my $got_vcf_stats = $hsc->_filter_vcf_and_count_snps($vcf_in, $tmp_vcf_out);
my %expected_vcf_stats = (
    chrom1 => {positions => 5, hets => 1, snps => 2},
    chrom2 => {positions => 1, hets => 1, snps => 1},
);
is_deeply($got_vcf_stats, \%expected_vcf_stats, '_filter_vcf_and_count_snps');
is(compare($tmp_vcf_out, $expected_vcf), 0, 'vcf made by _filter_vcf_and_count_snps is ok');
unlink $tmp_vcf_out;


my $fai = 't/data/het_snp_cal_length_from_fai.fai';
my %expected_lengths = (contig1 => 3900, contig2 => 4620, contig3 => 3540);
my $got_lengths = $hsc->_lengths_from_fai($fai);
is_deeply(\%expected_lengths, $got_lengths, '_lengths_from_fai');


my %summary_stats = (
    length => 0,
    positions => 0,
    snps => 0,
    hets => 0
);
my $expected_tsv = 't/data/het_snp_cal_summary_stats_all_zero';
my $tmp_tsv = 'tmp.test.HetSNPCalculator.summary_report.tsv';
$hsc->_write_summary_report($tmp_tsv, \%summary_stats);
is(compare($tmp_tsv, $expected_tsv), 0, 'Summary tsv ok when all zeros');
unlink $tmp_tsv;

$summary_stats{length} = 100;
$expected_tsv = 't/data/het_snp_cal_summary_stats_positions_zero';
$tmp_tsv = 'tmp.test.HetSNPCalculator.summary_report.tsv';
$hsc->_write_summary_report($tmp_tsv, \%summary_stats);
is(compare($tmp_tsv, $expected_tsv), 0, 'Summary tsv ok when no reads mapped');
unlink $tmp_tsv;

$summary_stats{positions} = 50;
$summary_stats{snps} = 10;
$expected_tsv = 't/data/het_snp_cal_summary_stats_no_hets';
$tmp_tsv = 'tmp.test.HetSNPCalculator.summary_report.tsv';
$hsc->_write_summary_report($tmp_tsv, \%summary_stats);
is(compare($tmp_tsv, $expected_tsv), 0, 'Summary tsv ok when no het snps');
unlink $tmp_tsv;

$summary_stats{hets} = 2;
$expected_tsv = 't/data/het_snp_cal_summary_stats_has_hets';
$tmp_tsv = 'tmp.test.HetSNPCalculator.summary_report.tsv';
$hsc->_write_summary_report($tmp_tsv, \%summary_stats);
is(compare($tmp_tsv, $expected_tsv), 0, 'Summary tsv ok when has het snps');
unlink $tmp_tsv;


my %ctg_lengths = (
    ctg1 => 100,
    ctg2 => 200,
    ctg3 => 300,
);
my %snp_stats = (
    'ctg1' => {positions => 50, snps => 2, hets => 1},
    'ctg2' => {positions => 100, snps => 0, hets => 0},
);

my $tmp_tsv_per_contig = 'tmp.test.HetSNPCalculator.write_report.per_contig.tsv';
my $tmp_tsv_totals = 'tmp.test.HetSNPCalculator.write_report.totals.tsv';
$hsc->_write_reports($tmp_tsv_per_contig, $tmp_tsv_totals, \%snp_stats, \%ctg_lengths);
my $expected_per_contig = 't/data/het_snp_cal_write_report.per_contig.tsv';
my $expected_totals = 't/data/het_snp_cal_write_report.totals.tsv';
is(compare($tmp_tsv_per_contig, $expected_per_contig), 0, 'write_reports per contig file ok');
is(compare($tmp_tsv_totals, $expected_totals), 0, 'write_reports totals file ok');
unlink $tmp_tsv_per_contig;
unlink $tmp_tsv_totals;


$hsc->run();
my $tmp_bcf = $hsc->outprefix . '.bcf';
is((-e $tmp_bcf and -s $tmp_bcf > 0), 1, 'run() made a non-empty bcf');
unlink $tmp_bcf;
my $seq_breakdown = $hsc->outprefix . '_ref_seq_breakdown.tsv';
is(compare($seq_breakdown, 't/data/het_snp_cal_run_ref_seq_breakdown.tsv'), 0, 'run() ref_seq_breakdown.tsv ok');
unlink $seq_breakdown;
my $summary_file = $hsc->outprefix . '_report.txt';
is(compare($summary_file, 't/data/het_snp_cal_run_report.txt'), 0, 'run() report.txt ok');
unlink $summary_file;

done_testing();
