#!/usr/bin/env perl
use strict;
use warnings;

BEGIN {
    use Test::Most ;
    use_ok('Pathogens::Parser::GenomeCoverage');
}

ok my $gc = Pathogens::Parser::GenomeCoverage->new( bamcheck => 't/data/small_slice.bam.bc' ), 'create instance of object';


# test errors
ok my $gc_err = Pathogens::Parser::GenomeCoverage->new( bamcheck => 'a_file_that_does_not_exist' ), 'create instance for nonfile';
throws_ok {$gc_err->coverage()} qr/GenomeCoverage::coverage/,'coverage throws for nonexistent file';
throws_ok {$gc_err->coverage_depth()} qr/GenomeCoverage::coverage/,'coverage_depth throws for nonexistent file';
throws_ok {$gc_err->ref_size('xxx')} qr/GenomeCoverage::ref_size/,'throws for garbage input';


# test coverage()
ok my @cover_bins = $gc->coverage(), 'coverage default runs';
is $cover_bins[0], '248', 'coverage default gives expected result';
ok @cover_bins = $gc->coverage(1,50,100,150,200), 'coverage for bins runs'; 
is join(',',@cover_bins), '248,221,194,165,145', 'coverage for bins gives expected result';
ok @cover_bins = $gc->coverage(1,150,100,50,200), 'coverage for unsorted bins runs'; 
is join(',',@cover_bins), '248,165,194,221,145', 'coverage for unsorted bins gives expected result';


# test coverage_depth()
throws_ok {$gc->coverage_depth} qr/Reference size must be set for mean coverage depth calculation./,'coverage_depth throws if ref_size not set.'; 
ok $gc->ref_size(200), 'set ref_size too low.';
throws_ok {$gc->coverage_depth} qr/Total bases found exceeds size of reference sequence./,'coverage_depth throws if ref_size too low';
ok $gc->ref_size(300), 'set ref_size to correct value.';
ok my ($cover_bases,$depth_mean,$depth_sd) = $gc->coverage_depth(), 'coverage_depth runs';
is sprintf("%d,%.2f,%.2f",$cover_bases,$depth_mean,$depth_sd),'248,167.00,118.05','coverage_depth gives expected result';
is $gc->_bam_coverage, 248, 'Bases covered check from samtools gives expected result';

done_testing();
exit;
