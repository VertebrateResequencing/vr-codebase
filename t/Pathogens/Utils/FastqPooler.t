#!/usr/bin/env perl
use strict;
use warnings;
use File::Compare;
use File::Path;

BEGIN {
    use Test::Most;
    use_ok('Pathogens::Utils::FastqPooler');
}


my $unpaired_path1 = 't/data/FastqPooler_unpaired1/1234';
my $unpaired_path2 = 't/data/FastqPooler_unpaired2/1235';
my $paired_path1 = 't/data/FastqPooler_paired1/1234';
my $paired_path2 = 't/data/FastqPooler_paired2/1235';

my @expected_fwd = ($unpaired_path1 . "_1.fastq.gz");
my @expected_rev;
my ($got_fwd, $got_rev) = Pathogens::Utils::FastqPooler::_get_files_to_cat([$unpaired_path1]);
is_deeply($got_fwd, \@expected_fwd, '_get_files_to_cat unpaired 1 fwd');
is_deeply($got_rev, \@expected_rev, '_get_files_to_cat unpaired 1 rev');

push(@expected_fwd, $unpaired_path2 . "_1.fastq.gz");
($got_fwd, $got_rev) = Pathogens::Utils::FastqPooler::_get_files_to_cat([$unpaired_path1, $unpaired_path2]);
is_deeply($got_fwd, \@expected_fwd, '_get_files_to_cat unpaired 1 and 2 fwd');
is_deeply($got_rev, \@expected_rev, '_get_files_to_cat unpaired 1 and 2 rev');

@expected_fwd = ($paired_path1 . "_1.fastq.gz");
@expected_rev = ($paired_path1 . "_2.fastq.gz");
($got_fwd, $got_rev) = Pathogens::Utils::FastqPooler::_get_files_to_cat([$paired_path1]);
is_deeply($got_fwd, \@expected_fwd, '_get_files_to_cat paired 1 fwd');
is_deeply($got_rev, \@expected_rev, '_get_files_to_cat paired 1 rev');




my $infile1 = 't/data/FastqPooler_cat_files.in1.gz';
my $infile2 = 't/data/FastqPooler_cat_files.in2.gz';
my $expected1 = 't/data/FastqPooler_cat_files.expected1';
my $expected1_2 = 't/data/FastqPooler_cat_files.expected1_2';
my $tmp_out = 'tmp.test.FastqPooler.cat_files.gz';
Pathogens::Utils::FastqPooler::_cat_files([$infile1], $tmp_out);
is( compare( $tmp_out, $expected1 ), 0, '_cat_files 1' );
unlink($tmp_out);

Pathogens::Utils::FastqPooler::_cat_files([$infile1, $infile2], $tmp_out);
is( compare( $tmp_out, $expected1_2 ), 0, '_cat_files 1 and 2' );
unlink($tmp_out);


my $tmp_out = 'tmp.test.FastqPooler';
if (-d $tmp_out) {
    File::Path::remove_tree($tmp_out);
}

ok my $pooler = Pathogens::Utils::FastqPooler->new(
    output_directory => $tmp_out,
    lane_paths => [$unpaired_path1]
), 'Create instance of object';
$pooler->run();
is( compare( "$tmp_out/forward.fastq", "t/data/FastqPooler_unpaired1_expected.fastq" ), 0, 'Unpaired 1 fastq' );
File::Path::remove_tree($tmp_out);


ok $pooler = Pathogens::Utils::FastqPooler->new(
    output_directory => $tmp_out,
    lane_paths => [$unpaired_path1, $unpaired_path2]
), 'Create instance of object';
$pooler->run();
is( compare( "$tmp_out/forward.fastq", "t/data/FastqPooler_unpaired1_and_2_expected.fastq" ), 0, 'Unpaired 1 and 2 fastq' );
File::Path::remove_tree($tmp_out);


ok $pooler = Pathogens::Utils::FastqPooler->new(
    output_directory => $tmp_out,
    lane_paths => [$paired_path1]
), 'Create instance of object';
$pooler->run();
is( compare( "$tmp_out/forward.fastq", "t/data/FastqPooler_paired1_expected.forward.fastq" ), 0, 'Paired 1 forward fastq' );
is( compare( "$tmp_out/reverse.fastq", "t/data/FastqPooler_paired1_expected.reverse.fastq" ), 0, 'Paired 1 reverse fastq' );
File::Path::remove_tree($tmp_out);


ok $pooler = Pathogens::Utils::FastqPooler->new(
    output_directory => $tmp_out,
    lane_paths => [$paired_path1, $paired_path2]
), 'Create instance of object';
$pooler->run();
is( compare( "$tmp_out/forward.fastq", "t/data/FastqPooler_paired1_and_2_expected.forward.fastq" ), 0, 'Paired 1 and 2 forward fastq' );
is( compare( "$tmp_out/reverse.fastq", "t/data/FastqPooler_paired1_and_2_expected.reverse.fastq" ), 0, 'Paired 1 and 2 reverse fastq' );
File::Path::remove_tree($tmp_out);

done_testing();
