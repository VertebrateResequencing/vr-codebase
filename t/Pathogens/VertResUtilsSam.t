#!/usr/bin/env perl
use strict;
use warnings;

BEGIN {
    use Test::Most ;
    use_ok('VertRes::Utils::Sam');
}

ok my $sam_util = VertRes::Utils::Sam->new(java_memory => 10), 'made sam_util';
isa_ok $sam_util, 'VertRes::Base';

# data files
my $simple_bam = 't/data/simple.bam';
my $fastq_base = 'fastq_base';
my $fastq_name = 't/data/'.$fastq_base;
my $cover_bam  = 't/data/small_slice.bam';
my $not_bam    = 't/data/io_test.txt'; # not a bam file.


# test sam2fastq
VertRes::Utils::Sam->new(verbose => 0, quiet => 1, java_memory => 800)->bam2fastq($simple_bam, $fastq_base);
# check files produced
is -s $fastq_name.'_1.fastq',           149478,'fastq_base.fastq_1.fastq produced';
is -s $fastq_name.'_2.fastq',           135478,'fastq_base.fastq_2.fastq produced';
is -s $fastq_name.'_1.fastq.fastqcheck', 12214,'fastq_base.fastq_1.fastq.fastqcheck produced';
is -s $fastq_name.'_2.fastq.fastqcheck', 10198,'fastq_base.fastq_2.fastq.fastqcheck produced';
# check intermediate files removed when restarting
rename($fastq_name.'_1.fastq',$fastq_name.'_1.fastq.tmp.fastq');
rename($fastq_name.'_2.fastq',$fastq_name.'_2.fastq.tmp.fastq');
unlink($simple_bam.'.bc');
unlink($fastq_name.'_1.fastq.fastqcheck',$fastq_name.'_2.fastq.fastqcheck');
VertRes::Utils::Sam->new(verbose => 0, quiet => 1, java_memory => 800)->bam2fastq($simple_bam, $fastq_base);
is -s $fastq_name.'_1.fastq',          149478,'fastq_base.fastq_1.fastq produced after restart';
is -s $fastq_name.'_2.fastq',          135478,'fastq_base.fastq_2.fastq produced after restart';
is -e $fastq_name.'_1.fastq.tmp.fastq', undef,'fastq_base.fastq_1.fastq.tmp.fastq removed after restart';
is -e $fastq_name.'_2.fastq.tmp.fastq', undef,'fastq_base.fastq_2.fastq.tmp.fastq removed after restart';
# clean up
unlink($simple_bam.'.bc');
unlink($fastq_name.'_1.fastq',$fastq_name.'_2.fastq');
unlink($fastq_name.'_1.fastq.fastqcheck',$fastq_name.'_2.fastq.fastqcheck');


# test coverage and coverage_depth

open(my $copy_stderr, ">&STDERR"); open(STDERR, '>/dev/null'); # Redirect STDERR

# test coverage()
ok my @cover_bins = $sam_util->coverage($cover_bam), 'coverage default ran ok';
is $cover_bins[0], '248', 'coverage expected result ok';
ok @cover_bins = $sam_util->coverage($cover_bam,1,50,100,150,200), 'coverage for bins ran ok'; 
is join(',',@cover_bins), '248,220,194,165,141', 'coverage for bins gives expected result';
ok @cover_bins = $sam_util->coverage($cover_bam,1,150,100,50,200), 'coverage for unsorted bins ran ok'; 
is join(',',@cover_bins), '248,165,194,220,141', 'coverage for unsorted bins gives expected result';
throws_ok {$sam_util->coverage('a_file_that_does_not_exist')} qr/VertRes::Utils::Sam::coverage/,'coverage throws for nonexistent file';
throws_ok {$sam_util->coverage($not_bam)} qr/VertRes::Utils::Sam::coverage/,'coverage throws for bad file';
# test coverage_depth()
ok my ($cover_bases,$depth_mean,$depth_sd) = $sam_util->coverage_depth($cover_bam,300), 'coverage_depth runs ok';
is sprintf("%d,%.2f,%.2f",$cover_bases,$depth_mean,$depth_sd),'248,165.30,117.74','coverage_depth gives expected result';
throws_ok {$sam_util->coverage_depth($cover_bam,'foobar')} qr/VertRes::Utils::Sam::coverage_depth/, 'coverage_depth throws for garbage input';
throws_ok {$sam_util->coverage_depth($cover_bam,0)} qr/VertRes::Utils::Sam::coverage_depth/, 'coverage_depth throws for ref size = zero';
throws_ok {$sam_util->coverage_depth($cover_bam,200)} qr/VertRes::Utils::Sam::coverage_depth/, 'coverage_depth throws for ref size too small';
throws_ok {$sam_util->coverage_depth('a_file_that_does_not_exist',300)} qr/VertRes::Utils::Sam::coverage_depth/, 'coverage_depth throws for nonexistent file';
throws_ok {$sam_util->coverage_depth($not_bam,300)} qr/VertRes::Utils::Sam::coverage_depth/, 'coverage_depth throws for bad file';

close(STDERR); open(STDERR, ">&", $copy_stderr); # Restore STDERR

done_testing();
exit;
