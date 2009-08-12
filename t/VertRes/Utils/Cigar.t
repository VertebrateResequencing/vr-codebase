#!/usr/bin/perl -w
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 15;
    
    use_ok('VertRes::Utils::Cigar');
    use_ok('VertRes::IO');
}

my $cigar_util = VertRes::Utils::Cigar->new();
isa_ok $cigar_util, 'VertRes::Base';

# setup our input files
my $io = VertRes::IO->new();
my $fq1_file = $io->catfile('t', 'data', '2822_6_1_1000.fastq');
ok -s $fq1_file, 'first fastq file ready to test with';
my $fq2_file = $io->catfile('t', 'data', '2822_6_2_1000.fastq');
ok -s $fq2_file, 'first fastq file ready to test with';
my $cigar1_file = $io->catfile('t', 'data', 'ssaha1.cigar');
ok -s $cigar1_file, 'first cigar file ready to test with';
my $cigar2_file = $io->catfile('t', 'data', 'ssaha2.cigar');
ok -s $cigar2_file, 'second cigar file ready to test with';
my $temp_dir = $io->tempdir();
my $sam_out = $io->catfile($temp_dir, 'cigar_to_sam.sam');

my $hash = $cigar_util->cigar_by_reads($cigar1_file);
is keys %{$hash}, 881, 'cigar_by_reads returned the correct number of read keys for first cigar';
$hash = $cigar_util->cigar_by_reads($cigar2_file);
is keys %{$hash}, 759, 'cigar_by_reads returned the correct number of read keys for second cigar';

is $cigar_util->cigar_to_sam([$fq1_file, $fq2_file], [$cigar1_file, $cigar2_file], 2000, undef, $sam_out), 1596, 'cigar_to_sam claimed to write the correct number of entries';
$io->file($sam_out);
is $io->num_lines, 1596, 'cigar_to_sam generated sam with correct number of lines';

is $cigar_util->cigar_to_sam([$fq1_file], [$cigar1_file], 2000, undef, $sam_out), 855, 'cigar_to_sam claimed to write the correct number of entries - first only';
$io->file($sam_out);
is $io->num_lines, 855, 'cigar_to_sam generated sam with correct number of lines - first only';

is $cigar_util->cigar_to_sam([$fq2_file], [$cigar2_file], 2000, undef, $sam_out), 741, 'cigar_to_sam claimed to write the correct number of entries - second only';
$io->file($sam_out);
is $io->num_lines, 741, 'cigar_to_sam generated sam with correct number of lines - second only';

exit;
