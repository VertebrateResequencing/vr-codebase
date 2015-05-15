#!/usr/bin/env perl
use strict;
use warnings;
use File::Temp;
use File::Spec;
use File::Compare;
use Data::Dumper;

BEGIN {
    use Test::Most ;
    use_ok('Pathogens::QC::HetSNPCalculator');
}

my $samtools_exe = 'samtools-1.1.30';
my $bcftools_exe = 'bcftools-1.2';
my $fa_ref = 't/data/15360_1#1_het_test_staph_aureus_subsample_2000.fa';
my $reference_size = 2194961;
my $lane_path = 't/data';
my $lane = '15360_1#1';
my $sample_dir = '15360_1#1_qc-sample';
my $het_report = 'heterozygosity_report.txt';
my $min_rawReadDepth = 10;
my $min_hqNonRefBases = 5;
my $rawReadDepth_hqNonRefBases_ratio = 0.3;
my $min_qual = 20;
my $hqRefReads_hqAltReads_ratio = 0.3;

ok my $hsc = Pathogens::QC::HetSNPCalculator->new(
						  samtools => $samtools_exe,
						  bcftools => $bcftools_exe,
						  fa_ref => $fa_ref,
						  reference_size => $reference_size,
						  lane_path => $lane_path,
						  lane => $lane,
						  sample_dir => $sample_dir,
						  het_report => $het_report,
						  min_rawReadDepth => 10,
						  min_hqNonRefBases => 5,
						  rawReadDepth_hqNonRefBases_ratio => 0.3,
						  min_qual => 20,
						  hqRefReads_hqAltReads_ratio => 0.3,
						 ), 'create instance of object';

is ( $hsc->samtools, $samtools_exe, 'Samtools exe' );
is ( $hsc->bcftools, $bcftools_exe, 'Bcftools exe' );
is ( $hsc->fa_ref, $fa_ref, 'Fasta ref' );
is ( $hsc->reference_size, $reference_size, 'Reference size' );
is ( $hsc->lane_path, $lane_path, 'Lane path' );
is ( $hsc->lane, $lane, 'Lane' );
is ( $hsc->sample_dir, $sample_dir, 'Sample dir' );
is ( $hsc->het_report, $het_report, 'Het report file name' );
is ( $hsc->min_rawReadDepth, $min_rawReadDepth, 'DP threshold' );
is ( $hsc->min_hqNonRefBases, $min_hqNonRefBases, 'DV threshold' );
is ( $hsc->rawReadDepth_hqNonRefBases_ratio, $rawReadDepth_hqNonRefBases_ratio, 'DP/DV threshold' );
is ( $hsc->min_qual, $min_qual, 'Min QUAL threshold' );
is ( $hsc->hqRefReads_hqAltReads_ratio, $hqRefReads_hqAltReads_ratio, 'DP4 Ref/Alt threshold' );

is ( $hsc->full_path, 't/data/15360_1#1_qc-sample', 'Full path');
is ( $hsc->mpileup_command, q(samtools-1.1.30 mpileup -d 500 -t INFO/DPR,DV -C50 -ugf t/data/15360_1#1_het_test_staph_aureus_subsample_2000.fa t/data/15360_1#1_qc-sample/15360_1#1.bam | bgzip > t/data/15360_1#1_qc-sample/15360_1#1_temp_vcf.vcf.gz), 'mpileup command');
is ( $hsc->total_number_of_snps_command, q(bcftools-1.2 query -f "%CHROM\n" -i "DP > 0" t/data/15360_1#1_qc-sample/15360_1#1_temp_vcf.vcf.gz > t/data/15360_1#1_qc-sample/15360_1#1_total_number_of_snps.csv), 'total number of snps command');
is ( $hsc->snp_call_command, q(bcftools-1.2 call -vm -O z t/data/15360_1#1_qc-sample/15360_1#1_temp_vcf.vcf.gz > t/data/15360_1#1_qc-sample/15360_1#1_snp_called.vcf.gz), 'snp call command');
is ( $hsc->bcf_query_command, q(bcftools-1.2 query -f "%CHROM %POS\n" -i "MIN(DP) >= 10 & MIN(DV) >= 5 & MIN(DV/DP)>= 0.3 & QUAL >= 20 & (GT='1/0' | GT='0/1' | GT='1/2') & ((DP4[0]+DP4[1])/(DP4[2]+DP4[3]) > 0.3)" t/data/15360_1#1_qc-sample/15360_1#1_snp_called.vcf.gz > t/data/15360_1#1_qc-sample/15360_1#1_filtered_snp_called_list.csv), 'bcf query filter command');

my $csv_file_with_empty_lines = 't/data/15360_1#1_total_number_of_snps_empty_lines.csv';
open (my $fh, '<', $csv_file_with_empty_lines) or die "$csv_file_with_empty_lines: $!";
is ( $hsc->_count_file_rows($fh), 4, 'total number of SNPs');

is ( $hsc->_calculate_percentage(1,400), 0.25, 'percentage calculation');

throws_ok { $hsc->_calculate_percentage(1,0) } qr/The value used as total/ , 'Illegal division';

throws_ok { $hsc->get_total_number_of_snps } qr/Backtrace:/, 'Total number of SNPs file doesnt exist yet';

ok ( $hsc->get_number_of_het_snps, 'Calculating number of heterozygous SNPs');
is ( $hsc->number_of_het_snps, '1', 'Number of heterozigous SNPs');

ok ( $hsc->get_total_number_of_snps, 'Calculating total number of SNPs');

#Expected files
my $temp_master_file = File::Spec->catfile('t/data','15360_1#1_temp_vcf_master.vcf.gz');
my $snp_called_master_file = File::Spec->catfile('t/data','15360_1#1_snp_called_master.vcf.gz');
my $filtered_master_file = File::Spec->catfile('t/data','15360_1#1_filtered_snp_called_list_master.csv');
my $total_number_of_snps_master_file = File::Spec->catfile('t/data','15360_1#1_total_number_of_snps_master.csv');

#Actual output
my $temp_vcf_file = File::Spec->catfile('t/data/15360_1#1_qc-sample','15360_1#1_temp_vcf.vcf.gz');
my $snp_called_vcf_file = File::Spec->catfile('t/data/15360_1#1_qc-sample','15360_1#1_snp_called.vcf.gz');
my $filtered_csv_file = File::Spec->catfile('t/data/15360_1#1_qc-sample','15360_1#1_filtered_snp_called_list.csv');
my $total_number_of_snps_csv_file = File::Spec->catfile('t/data/15360_1#1_qc-sample','15360_1#1_total_number_of_snps.csv');

is(compare($temp_vcf_file,$temp_master_file),0, 'Temp vcf file');
is(compare($snp_called_vcf_file,$snp_called_master_file),0, 'Snp called vcf file');
is(compare($filtered_csv_file,$filtered_master_file),0, 'Filtered csv file');
is(compare($total_number_of_snps_csv_file,$total_number_of_snps_master_file),0, 'Total number of snps (Het and Hom)');

is ( $hsc->total_number_of_snps, '16893', 'Total number of SNPs (Het and Hom)');

ok ($hsc->get_percentages_of_het_snps, 'Calculating percentages');
is ( $hsc->het_snps_genome_percentage, '4.55588960350548e-05', 'Percentage of heterozygous SNPs for the whole genome');
is ( $hsc->het_snps_total_snps_percentage, '0.00591961167347422', 'Percentage of heterozygous SNPs for all the SNPs called (Het and Hom)');


ok ( $hsc->write_het_report, 'Writing heterozygosity report' );
is ( $hsc->het_report_path, 't/data/heterozygosity_report.txt', 'Het report path' );

my $het_report_file = $hsc->het_report_path;
open (my $fh2, '<', $het_report_file) or die "$het_report_file: $!";
my @lines = <$fh2>;
close($fh2);
is (scalar @lines, 1, 'Only one line in the het report');

chomp($lines[0]);
is ($lines[0], "1\t4.55588960350548e-05\t0.00591961167347422", 'Content of het report file');

my $het_report_master_file = File::Spec->catfile('t/data','15360_1#1_heterozygosity_report_master.txt');
is(compare($het_report_file,$het_report_master_file),0, 'Heterozygosity report file');

if (-e $het_report_file) {
  unlink($het_report_file);
}

ok ( $hsc->remove_temp_vcfs_and_csvs, 'Removing files created');


done_testing();
