#!/usr/bin/env perl
use strict;
use warnings;
use Test::MockObject;
use Test::Exception;

# does it compile?
BEGIN { unshift(@INC, './modules') }
BEGIN {
    use Test::Most;
    use_ok('Pathogens::Import::CompressAndValidate');
    my $valid_fastq_wrapper = Test::MockObject->new();
    $valid_fastq_wrapper->fake_module( 'Pathogens::Import::ValidateFastqConversion', test => sub{1} );
    $valid_fastq_wrapper->fake_new( 'Pathogens::Import::ValidateFastqConversion' );
    $valid_fastq_wrapper->mock('is_total_reads_valid', sub{ 1 });
}

# Pass in one fastq file - Do we get expected output files?
my $fastq_1 = 't/data/SRR001629.99999_1.fastq';
my $fastq_2 = 't/data/SRR001629.99999_2.fastq';
my $nonfile = 't/data/does_not_exist.fastq';
ok my $validator = Pathogens::Import::CompressAndValidate->new(
    irods_filename  => 'a_file_that_does_not_exist',
    fastq_filenames => [$fastq_1]
    ), 'Object instance created.';
is $validator->is_compressed_and_validated(), 1, 'Compress and validate completes for single fastq.';

# Check files produced.
is -e $fastq_1.'.gz',                    1, 'Compressed fastq file exists';
is -s $fastq_1.'.gz',                  329, 'Compressed fastq file correct size';
is -e $fastq_1.'.md5',                   1, 'Uncompressed fastq md5 file exists';
is -e $fastq_1.'.gz.md5',                1, 'Compressed fastq md5 file exists';
is -e $fastq_1.'.gz.fastqcheck',         1, 'Final fastqcheck file exists';
is -s $fastq_1.'.gz.fastqcheck',     28700, 'Final fastqcheck file correct size';
is -e $fastq_1.'.gz.fastqcheck.tmp', undef, 'Fastqcheck.tmp file removed';

# Check compressed reads
open(my $fh_1, "gunzip -c $fastq_1.gz | wc -l |");
is <$fh_1>, "12\n", 'Expected number of reads from file.';
close $fh_1;

# Check md5 files.
open my $fh_2, $fastq_1.'.md5';
is <$fh_2>, "1a5665ea34f204387f5b2513d7e8e004  $fastq_1\n",    'Expected checksum from fastq.md5 file.';
close $fh_2;
open my $fh_3, $fastq_1.'.gz.md5';
is <$fh_3>, "ec3547080cb487060469fafc368c9ce7  $fastq_1.gz\n", 'Expected checksum from compressed fastq.md5 file.';
close $fh_3;

# Unlink temp files
unlink($fastq_1.'.gz',$fastq_1.'.md5',$fastq_1.'.gz.md5',$fastq_1.'.gz.fastqcheck',$fastq_1.'.gz.fastqcheck.tmp');


# Pass in two fastq files - do we get expected output files?
ok my $validator_pair = Pathogens::Import::CompressAndValidate->new(
    irods_filename  => 'a_file_that_does_not_exist',
    fastq_filenames => [$fastq_1,$fastq_2]
    ), 'Object instance created.';
is $validator_pair->is_compressed_and_validated(), 1, 'Compress and validate completes for paired fastq.';

# Check files produced.
is -e $fastq_1.'.gz',                    1, 'Compressed fastq file exists';
is -s $fastq_1.'.gz',                  329, 'Compressed fastq file correct size';
is -e $fastq_1.'.md5',                   1, 'Uncompressed fastq md5 file exists';
is -e $fastq_1.'.gz.md5',                1, 'Compressed fastq md5 file exists';
is -e $fastq_1.'.gz.fastqcheck',         1, 'Final fastqcheck file exists';
is -s $fastq_1.'.gz.fastqcheck',     28700, 'Final fastqcheck file correct size';
is -e $fastq_1.'.gz.fastqcheck.tmp', undef, 'Fastqcheck.tmp file removed';

# Check compressed reads
open(my $fh_4, "gunzip -c $fastq_1.gz | wc -l |");
is <$fh_4>, "12\n", 'Expected number of reads from file.';
close $fh_4;

# Check md5 files.
open my $fh_5, $fastq_1.'.md5';
is <$fh_5>, "1a5665ea34f204387f5b2513d7e8e004  $fastq_1\n",    'Expected checksum from fastq.md5 file.';
close $fh_5;
open my $fh_6, $fastq_1.'.gz.md5';
is <$fh_6>, "ec3547080cb487060469fafc368c9ce7  $fastq_1.gz\n", 'Expected checksum from compressed fastq.md5 file.';
close $fh_6;

# Check files produced.
is -e $fastq_2.'.gz',                    1, 'Compressed fastq file exists';
is -s $fastq_2.'.gz',                  325, 'Compressed fastq file correct size';
is -e $fastq_2.'.md5',                   1, 'Uncompressed fastq md5 file exists';
is -e $fastq_2.'.gz.md5',                1, 'Compressed fastq md5 file exists';
is -e $fastq_2.'.gz.fastqcheck',         1, 'Final fastqcheck file exists';
is -s $fastq_2.'.gz.fastqcheck',     16000, 'Final fastqcheck file correct size';
is -e $fastq_2.'.gz.fastqcheck.tmp', undef, 'Fastqcheck.tmp file removed';

# Check compressed reads
open(my $fh_7, "gunzip -c $fastq_2.gz | wc -l |");
is <$fh_7>, "12\n", 'Expected number of reads from file.';
close $fh_7;

# Check md5 files.
open my $fh_8, $fastq_2.'.md5';
is <$fh_8>, "b8cba4bdb437cb71ef6059cfa82aed24  $fastq_2\n",    'Expected checksum from fastq.md5 file.';
close $fh_8;
open my $fh_9, $fastq_2.'.gz.md5';
is <$fh_9>, "fd2c496f1832d4b8bc7510e6bee2fdb2  $fastq_2.gz\n", 'Expected checksum from compressed fastq.md5 file.';
close $fh_9;

# Unlink temp files for forward reads only
unlink($fastq_1.'.gz',$fastq_1.'.md5',$fastq_1.'.gz.md5',$fastq_1.'.gz.fastqcheck',$fastq_1.'.gz.fastqcheck.tmp');


# Re-run with output from previous files in place 
my $validator_rerun = Pathogens::Import::CompressAndValidate->new(
    irods_filename  => 'a_file_that_does_not_exist',
    fastq_filenames => [$fastq_1,$fastq_2]
    );
is $validator_rerun->is_compressed_and_validated(), 1, 'Compress and validate completes when output files already present.';


# Test garbage in error
my $validator_input_error = Pathogens::Import::CompressAndValidate->new(
    irods_filename  => 'a_file_that_does_not_exist',
    fastq_filenames => [$fastq_1,$nonfile]
    );
open(my $copy_stderr, ">&STDERR"); open(STDERR, '>/dev/null'); # Redirect STDERR
throws_ok { $validator_input_error->is_compressed_and_validated() } qr/Pathogens::Import::CompressAndValidate::_compress_and_checksum/, 'Compress and validate throws when input files not present.';
close(STDERR); open(STDERR, ">&", $copy_stderr); # Restore STDERR

# iRODS validation fail
my $invalid_fastq_wrapper = Test::MockObject->new();
$invalid_fastq_wrapper->fake_module( 'Pathogens::Import::ValidateFastqConversion', test => sub{1} );
$invalid_fastq_wrapper->fake_new( 'Pathogens::Import::ValidateFastqConversion' );
$invalid_fastq_wrapper->mock('is_total_reads_valid', sub{ 0 }); # Error!!

my $validator_error = Pathogens::Import::CompressAndValidate->new(
    irods_filename  => 'a_file_that_does_not_exist',
    fastq_filenames => [$fastq_1,$fastq_2]
    );
is $validator_error->is_compressed_and_validated(), 0, 'Compress and validate returns zero on error.';


# Unlink temp files
unlink($fastq_1.'.gz',$fastq_1.'.md5',$fastq_1.'.gz.md5',$fastq_1.'.gz.fastqcheck',$fastq_1.'.gz.fastqcheck.tmp');
unlink($fastq_2.'.gz',$fastq_2.'.md5',$fastq_2.'.gz.md5',$fastq_2.'.gz.fastqcheck',$fastq_2.'.gz.fastqcheck.tmp');
unlink($nonfile.'.gz');

done_testing();
