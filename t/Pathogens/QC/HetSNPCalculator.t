#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;

BEGIN {
    use Test::Most ;
    use_ok('Pathogens::QC::HetSNPCalculator');
}

my $samtools_exe = 'samtools-1.1.30';
my $bcftools_exe = 'bcftools-1.2';
my $fa_ref = 'path/to_refs/Neisseria_meningitidis_serogroup_C_FAM18_v1.fa';
my $reference_size = 2194961;
my $lane_path = '/lane/path';
my $lane = '1234_4#2';
my $sample_dir = 'qc-sample';
my $het_report = 'heterozygosity_report.txt';
my $bcft_min_dp = 10;
my $bcft_min_dv = 5;
my $bcft_dp_dv_ratio = 0.3;
my $bcft_min_qual = 20;
my $bcft_dp4_ref_allele_ratio = 0.3;

ok my $hsc = Pathogens::QC::HetSNPCalculator->new(
						  samtools => $samtools_exe,
						  bcftools => $bcftools_exe,
						  fa_ref => $fa_ref,
						  reference_size => $reference_size,
						  lane_path => $lane_path,
						  lane => $lane,
						  sample_dir => $sample_dir,
						  heterozygosity_report_file_name => $het_report,
						  bcft_min_dp => 10,
						  bcft_min_dv => 5,
						  bcft_dp_dv_ratio => 0.3,
						  bcft_min_qual => 20,
						  bcft_dp4_ref_allele_ratio => 0.3,
						 ), 'create instance of object';

is ( $hsc->samtools, $samtools_exe, 'Samtools exe' );
is ( $hsc->bcftools, $bcftools_exe, 'Bcftools exe' );
is ( $hsc->fa_ref, $fa_ref, 'Fasta ref' );
is ( $hsc->reference_size, $reference_size, 'Reference size' );
is ( $hsc->lane_path, $lane_path, 'Lane path' );
is ( $hsc->lane, $lane, 'Lane' );
is ( $hsc->sample_dir, $sample_dir, 'Sample dir' );
is ( $hsc->heterozygosity_report_file_name, $het_report, 'Het report file name' );
is ( $hsc->bcft_min_dp, $bcft_min_dp, 'DP threshold' );
is ( $hsc->bcft_min_dv, $bcft_min_dv, 'DV threshold' );
is ( $hsc->bcft_dp_dv_ratio, $bcft_dp_dv_ratio, 'DP/DV threshold' );
is ( $hsc->bcft_min_qual, $bcft_min_qual, 'Min QUAL threshold' );
is ( $hsc->bcft_dp4_ref_allele_ratio, $bcft_dp4_ref_allele_ratio, 'DP4 Ref/Alt threshold' );

is ( $hsc->full_path, '/lane/path/qc-sample', 'Full path');
is ( $hsc->mpileup_command, q(samtools-1.1.30 mpileup -d 500 -t INFO/DPR,DV -C50 -ugf path/to_refs/Neisseria_meningitidis_serogroup_C_FAM18_v1.fa /lane/path/qc-sample/1234_4#2.bam | bgzip > /lane/path/qc-sample/1234_4#2_temp_vcf.vcf.gz), 'mpileup command');
is ( $hsc->total_number_of_snps_command, q(bcftools-1.2 call -m -f GQ,GP /lane/path/qc-sample/1234_4#2_temp_vcf.vcf.gz | egrep -v "^#|DP=0" | wc -l), 'total number of snps command');
is ( $hsc->snp_call_command, q(bcftools-1.2 call -vm -O z /lane/path/qc-sample/1234_4#2_temp_vcf.vcf.gz > /lane/path/qc-sample/1234_4#2_snp_called.vcf.gz), 'snp call command');
is ( $hsc->bcf_query_filter_command, q(bcftools-1.2 filter -i "(DP4[0]+DP4[1])/(DP4[2]+DP4[3]) > 0.3" /lane/path/qc-sample/1234_4#2_snp_called.vcf.gz | bcftools-1.2 filter -i "MIN(DP) >= 10 & MIN(DV) >= 5 & MIN(DV/DP)>= 0.3 & QUAL >= 20 & (GT='1/0' | GT='0/1' | GT='1/2')" - > /lane/path/qc-sample/1234_4#2_filtered_snp_called.vcf), 'bcf query filter command');
print Dumper($hsc);

done_testing();
