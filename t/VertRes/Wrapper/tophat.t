#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use File::Spec;

BEGIN {
    use Test::Most;
    
    use_ok('VertRes::Wrapper::tophat');
    use_ok('VertRes::Utils::Mappers::tophat');
    use VertRes::Wrapper::fastqcheck;
    use VertRes::Utils::FileSystem;
}


my $tophat = VertRes::Wrapper::tophat->new(quiet => 1);
isa_ok $tophat, 'VertRes::Wrapper::WrapperI';
is $tophat->quiet, 1, 'quiet set via new';

is $tophat->exe, 'tophat', 'exe ok';
like $tophat->version, qr/\d+\.\d+/, 'version ok';


# prepare our test files; copy them to a temp dir where we will do the tests
my $fsu = VertRes::Utils::FileSystem->new();
my $temp_dir = $fsu->tempdir;
my $read1_basename = '2822_6_1_1000.fastq';
my $read2_basename = '2822_6_2_1000.fastq';
my $ref_basename = 'S_suis_P17.fa';
my $read1 = File::Spec->catfile($temp_dir, $read1_basename);
my $read2 = File::Spec->catfile($temp_dir, $read2_basename);
my $ref = File::Spec->catfile($temp_dir, $ref_basename);
copy(File::Spec->catfile('t', 'data', $read1_basename), $read1);
copy(File::Spec->catfile('t', 'data', $read2_basename), $read2);
copy(File::Spec->catfile('t', 'data', $ref_basename), $ref);
ok -s $read1, 'test file 1 ready to use';
ok -s $read2, 'test file 2 ready to use';
ok -s $ref, 'test file 3 ready to use';

# the files we expect to be created
my $sai1 = File::Spec->catfile($temp_dir, '2822_6_1_1000.sai');
my $sai2 = File::Spec->catfile($temp_dir, '2822_6_2_1000.sai');
my $mapping = File::Spec->catfile($temp_dir, 'mapping.sam');
my @ref_index_files;
foreach my $suffix (qw(1.ebwt 2.ebwt 3.ebwt 4.ebwt rev.1.ebwt rev.2.ebwt)) {
    push(@ref_index_files, File::Spec->catfile($temp_dir, 'S_suis_P17.fa.'.$suffix));
}

my $fw = VertRes::Wrapper::fastqcheck->new(quiet => 1);
my $out_file1 = File::Spec->catfile($temp_dir, '2822_6_1_1000.fastq.fastqcheck');
my $out_file2 = File::Spec->catfile($temp_dir, '2822_6_2_1000.fastq.fastqcheck');
$fw->run($read1, $out_file1);
$fw->run($read2, $out_file2);


# run the whole mapping
$tophat->do_mapping(ref => $ref,
                 read1 => $read1,
                 read2 => $read2,
                 output => $mapping,
                 insert_size => 500);
is $tophat->run_status, 2, 'status after mapping is ok';
my ($lines, $mapped) = check_sam($mapping);
is $lines, 524, 'output sam not truncated';
cmp_ok $mapped, '>=', 522, 'mapped enough reads';
foreach my $file ($mapping, @ref_index_files) {
    ok -s $file, "output file exists $file";
    unlink($file);
}

# individual commands
$tophat->setup_reference($ref);
is $tophat->run_status, 2, 'status after index is ok';
foreach my $file (@ref_index_files) {
    ok -s $file, "output file exists $file";
}

done_testing;
exit;

sub check_sam {
    my $sam = shift;
    open(my $sfh, $sam) || return (0, 0);
    my $lines = 0;
    my $mapped = 0;
    while (<$sfh>) {
        next unless $_ =~ /\S/;
        $lines++;
        my @a = split;
        $mapped++ if $a[3];
    }
    return ($lines, $mapped);
}

