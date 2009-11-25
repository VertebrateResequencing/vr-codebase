#!/usr/bin/perl -w
use strict;
use warnings;
use File::Spec;

BEGIN {
    use Test::Most tests => 10;
    
    use_ok('VertRes::Wrapper::fastqcheck');
    use_ok('VertRes::Utils::FileSystem');
}

my $fw = VertRes::Wrapper::fastqcheck->new(quiet => 1);
isa_ok $fw, 'VertRes::Wrapper::WrapperI';
is $fw->quiet, 1, 'quiet set via new';

# setup files
my $io = VertRes::IO->new();
my $fsu = VertRes::Utils::FileSystem->new();
my $fastq_file1 = File::Spec->catfile('t', 'data', '2822_6_1_1000.fastq');
ok -s $fastq_file1, 'first input fastq file ready to test on';
my $fastq_file2 = File::Spec->catfile('t', 'data', 'SRR001629_1.fastq.gz');
ok -s $fastq_file2, 'second input fastq file ready to test on';
my $temp_dir = $fsu->tempdir();
my $out_file = File::Spec->catfile($temp_dir, 'fastqcheck.file');

# test on an uncompressed fastq with no spaces in the ids
$fw->run($fastq_file1, $out_file);
cmp_ok $fw->run_status, '>=', 1, 'ran ok on an easy fastq';
$io->file($out_file);
is $io->num_lines(), 65, 'fastqcheck output had correct number of lines';

# test on a compressed fastq with spaces in the ids
unlink($out_file);
$fw->run($fastq_file2, $out_file);
cmp_ok $fw->run_status, '>=', 1, 'ran ok on a hard fastq';
$io->file($out_file);
is $io->num_lines($out_file), 243, 'fastqcheck output had correct number of lines';

exit;
