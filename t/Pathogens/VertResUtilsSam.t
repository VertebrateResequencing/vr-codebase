#!/usr/bin/env perl
use strict;
use warnings;

BEGIN {
    use Test::Most ;
    use_ok('VertRes::Utils::Sam');
}

ok my $sam_util = VertRes::Utils::Sam->new(java_memory => 10), 'made sam_util';
isa_ok $sam_util, 'VertRes::Base';

# test coverage and coverage_depth
my $bam     = 't/data/small_slice.bam';
my $not_bam = 't/data/io_test.txt'; # not a bam file.

open(my $copy_stderr, ">&STDERR"); open(STDERR, '>/dev/null'); # Redirect STDERR
# test coverage()
ok my @cover_bins = $sam_util->coverage($bam), 'coverage default ran ok';
is $cover_bins[0], '248', 'coverage expected result ok';
ok @cover_bins = $sam_util->coverage($bam,1,50,100,150,200), 'coverage for bins ran ok'; 
is join(',',@cover_bins), '248,204,158,84,0', 'coverage for bins gives expected result';
ok @cover_bins = $sam_util->coverage($bam,1,150,100,50,200), 'coverage for unsorted bins ran ok'; 
is join(',',@cover_bins), '248,84,158,204,0', 'coverage for unsorted bins gives expected result';
throws_ok {$sam_util->coverage('a_file_that_does_not_exist')} qr/VertRes::Utils::Sam::coverage/,'coverage throws for nonexistent file';
throws_ok {$sam_util->coverage($not_bam)} qr/VertRes::Utils::Sam::coverage/,'coverage throws for bad file';
# test coverage_depth()
ok my ($cover_bases,$depth_mean,$depth_sd) = $sam_util->coverage_depth($bam,300), 'coverage_depth runs ok';
is sprintf("%d,%.2f,%.2f",$cover_bases,$depth_mean,$depth_sd),'248,93.95,64.72','coverage_depth gives expected result';
throws_ok {$sam_util->coverage_depth($bam,'foobar')} qr/VertRes::Utils::Sam::coverage_depth/, 'coverage_depth throws for garbage input';
throws_ok {$sam_util->coverage_depth($bam,0)} qr/VertRes::Utils::Sam::coverage_depth/, 'coverage_depth throws for ref size = zero';
throws_ok {$sam_util->coverage_depth($bam,200)} qr/VertRes::Utils::Sam::coverage_depth/, 'coverage_depth throws for ref size too small';
throws_ok {$sam_util->coverage_depth('a_file_that_does_not_exist',300)} qr/VertRes::Utils::Sam::coverage_depth/, 'coverage_depth throws for nonexistent file';
throws_ok {$sam_util->coverage_depth($not_bam,300)} qr/VertRes::Utils::Sam::coverage_depth/, 'coverage_depth throws for bad file';
close(STDERR); open(STDERR, ">&", $copy_stderr); # Restore STDERR

done_testing();
exit;
