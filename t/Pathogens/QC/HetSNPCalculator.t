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

my $samtools_exe   = 'samtools-1.3';
my $bcftools_exe   = 'bcftools-1.3';
my $fa_ref         = 't/data/15360_1#1_het_test_staph_aureus_subsample_2000.fa';
my $reference_size = 2194961;
my $lane_path      = 't/data';
my $lane           = '15360_1#1';
my $sample_dir     = '15360_1#1_qc-sample';
my $het_report     = 'heterozygous_snps_report.txt';
my $min_rawReadDepth                 = 10;
my $min_hqNonRefBases                = 5;
my $rawReadDepth_hqNonRefBases_ratio = 0.3;
my $min_qual                         = 20;
my $hqRefReads_hqAltReads_ratio      = 0.3;

ok my $hsc = Pathogens::QC::HetSNPCalculator->new(
    samtools                         => $samtools_exe,
    bcftools                         => $bcftools_exe,
    fa_ref                           => $fa_ref,
    reference_size                   => $reference_size,
    lane_path                        => $lane_path,
    lane                             => $lane,
    sample_dir                       => $sample_dir,
    het_report                       => $het_report,
    min_rawReadDepth                 => 10,
    min_hqNonRefBases                => 5,
    rawReadDepth_hqNonRefBases_ratio => 0.3,
    min_qual                         => 20,
    hqRefReads_hqAltReads_ratio      => 0.3,
  ),
  'create instance of object';

#is( $hsc->samtools,          $samtools_exe,      'Samtools exe' );
#is( $hsc->bcftools,          $bcftools_exe,      'Bcftools exe' );
#is( $hsc->fa_ref,            $fa_ref,            'Fasta ref' );
#is( $hsc->reference_size,    $reference_size,    'Reference size' );
#is( $hsc->lane_path,         $lane_path,         'Lane path' );
#is( $hsc->lane,              $lane,              'Lane' );
#is( $hsc->sample_dir,        $sample_dir,        'Sample dir' );
#is( $hsc->het_report,        $het_report,        'Het report file name' );
#is( $hsc->min_rawReadDepth,  $min_rawReadDepth,  'DP threshold' );
#is( $hsc->min_hqNonRefBases, $min_hqNonRefBases, 'DV threshold' );
#is(
#    $hsc->rawReadDepth_hqNonRefBases_ratio,
#    $rawReadDepth_hqNonRefBases_ratio,
#    'DP/DV threshold'
#);
#is( $hsc->min_qual, $min_qual, 'Min QUAL threshold' );
#is( $hsc->hqRefReads_hqAltReads_ratio,
#    $hqRefReads_hqAltReads_ratio, 'DP4 Ref/Alt threshold' );
#
#is( $hsc->full_path, 't/data/15360_1#1_qc-sample', 'Full path' );
#is(
#    $hsc->mpileup_command,
#q(samtools-1.1.30 mpileup --skip-indels -d 500 -t INFO/AD,INFO/ADF,INFO/ADR -C50 -u -f t/data/15360_1#1_het_test_staph_aureus_subsample_2000.fa t/data/15360_1#1_qc-sample/15360_1#1.bam > t/data/15360_1#1_qc-sample/15360_1#1_temp_vcf.vcf),
#    'mpileup command'
#);
#is(
#    $hsc->total_genome_covered_command,
#q(bcftools-1.2 query -f "%CHROM\n" -i "DP > 0" t/data/15360_1#1_qc-sample/15360_1#1_temp_vcf.vcf > t/data/15360_1#1_qc-sample/15360_1#1_total_genome_covered.csv),
#    'total number of snps command'
#);
#is(
#    $hsc->snp_call_command,
#q(bcftools-1.2 call -vm -O z t/data/15360_1#1_qc-sample/15360_1#1_temp_vcf.vcf > t/data/15360_1#1_qc-sample/15360_1#1_snp_called.vcf.gz),
#    'snp call command'
#);
#is(
#   $hsc->all_snps_command,
#   q(bcftools-1.2 query -f "%CHROM %POS\n" -i "MIN(DP) >= 10 & MIN(DV) >= 5 & MIN(DV/DP)>= 0.3 & QUAL >= 20 & (GT='0/0' | GT='1/1' | GT='0/1' | GT='1/2')" t/data/15360_1#1_qc-sample/15360_1#1_snp_called.vcf.gz > t/data/15360_1#1_qc-sample/15360_1#1_all_snps_list.csv),
#   'all snps command'
#);
#is(
#    $hsc->bcf_query_command,
#q(bcftools-1.2 query -f "%CHROM %POS\n" -i "MIN(DP) >= 10 & MIN(DV) >= 5 & MIN(DV/DP)>= 0.3 & QUAL >= 20 & (GT='0/1' | GT='1/2') & ((DP4[0]+DP4[1])/(DP4[2]+DP4[3]) > 0.3)" t/data/15360_1#1_qc-sample/15360_1#1_snp_called.vcf.gz > t/data/15360_1#1_qc-sample/15360_1#1_filtered_snp_called_list.csv),
#    'bcf query filter command'
#);
#
#
##throws_ok { $hsc->total_genome_covered } qr/Backtrace:/,
##  'Total number of SNPs file doesnt exist yet';
#
#is( $hsc->total_genome_covered, 16893, 'Total number of SNPs (Het and Hom)' );
#
##Expected files
#my $temp_master_file =
#  File::Spec->catfile( 't/data', '15360_1#1_temp_vcf_master.vcf' );
#my $snp_called_master_file =
#  File::Spec->catfile( 't/data', '15360_1#1_snp_called_master.vcf.gz' );
#my $filtered_master_file = File::Spec->catfile( 't/data',
#    '15360_1#1_filtered_snp_called_list_master.csv' );
#
#is( compare( $hsc->temp_vcf, $temp_master_file ), 0, 'Temp vcf file' );
#is( compare( $hsc->snp_called_vcf, $snp_called_master_file ),
#    0, 'Snp called vcf file' );
#is( compare( $hsc->filtered_snp_called_csv, $filtered_master_file ),
#    0, 'Filtered csv file' );
#
#open( my $fh2, '<', $hsc->het_report_path )
#  or die "Couldn't open the het report file for reading";
#my @lines = <$fh2>;
#close($fh2);
#my @values = split( /\t/, $lines[1] );
#
#$values[2] =~ s/\n//;
#is( $values[0], '1', 'No. Het SNPs' );
#is( $values[1], '4.55588960350548e-05',
#    '% Het SNPs (Total Genome)' );
#is( $values[2], '0.00591961167347422',
#    '% Het SNPs (Genome Covered)' );
#chomp($values[3]);
#is( $values[3], '0.819672131147541',
#    '% Het SNPs (Total No. of SNPs)' );
#
#if ( -e $hsc->het_report_path ) {
#    unlink( $hsc->het_report_path );
#}
#
#
#
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


done_testing();
