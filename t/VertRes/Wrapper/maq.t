#!/usr/bin/perl -w
use strict;
use warnings;
use File::Copy;

BEGIN {
    use Test::Most tests => 36;
    
    use_ok('VertRes::Wrapper::maq');
    use_ok('VertRes::IO');
}

my $maq = VertRes::Wrapper::maq->new(quiet => 1);
isa_ok $maq, 'VertRes::Wrapper::WrapperI';
is $maq->quiet, 1, 'quiet set via new';

is $maq->exe, 'maq', 'exe ok';
like $maq->version, qr/\d+\.\d+/, 'version ok';

# prepare our test files; copy them to a temp dir where we will do the tests
my $io = VertRes::IO->new();
my $temp_dir = $io->tempdir;
my $read1_basename = '2822_6_1_1000.fastq';
my $read2_basename = '2822_6_2_1000.fastq';
my $bad_read1_basename = 'bad_1.fastq.gz';
my $bad_read2_basename = 'bad_2.fastq.gz';
my $ref_basename = 'S_suis_P17.fa';
my $read1 = $io->catfile($temp_dir, $read1_basename);
my $read2 = $io->catfile($temp_dir, $read2_basename);
my $bad_read1 = $io->catfile($temp_dir, $bad_read1_basename);
my $bad_read2 = $io->catfile($temp_dir, $bad_read2_basename);
my $ref = $io->catfile($temp_dir, $ref_basename);
copy($io->catfile('t', 'data', $read1_basename), $read1);
copy($io->catfile('t', 'data', $read2_basename), $read2);
copy($io->catfile('t', 'data', $bad_read1_basename), $bad_read1);
copy($io->catfile('t', 'data', $bad_read2_basename), $bad_read2);
copy($io->catfile('t', 'data', $ref_basename), $ref);
ok -s $read1, 'test file 1 ready to use';
ok -s $read2, 'test file 2 ready to use';
ok -s $bad_read1, 'bad test file 1 ready to use';
ok -s $bad_read2, 'bad test file 2 ready to use';
ok -s $ref, 'test file 3 ready to use';

# the files we expect to be created
my $bfa  = $ref.'.bfa';
my $bfq1 = $io->catfile($temp_dir, '2822_6_1_1000.fastq.bfq');
my $bfq2 = $io->catfile($temp_dir, '2822_6_2_1000.fastq.bfq');
my $map_file = $io->catfile($temp_dir, 'mapping.map');
my $unmapped = $map_file.'.unmapped';
my $mapstat = $map_file.'.mapstat';
my $sam_file = $io->catfile($temp_dir, 'mapping.sam');

# run the whole mapping
$maq->do_mapping(ref => $ref,
                 read1 => $read1,
                 read2 => $read2,
                 output => $sam_file,
                 a => 200);
is $maq->run_status, 1, 'status after mapping is ok';
my ($lines, $mapped) = check_sam($sam_file);
is $lines, 2000, 'output sam not truncated';
cmp_ok $mapped, '>=', 1670, 'mapped enough reads';
foreach my $file ($bfa, $unmapped, $sam_file, $map_file) {
    ok -s $file, 'do_mapping output file exists';
    unlink($file);
}

# individual commands
$maq->fasta2bfa($ref, $bfa);
is $maq->run_status, 1, 'status after fasta2bfa is ok';
ok -s $bfa, 'bfa file exists';

$maq->fastq2bfq($read1, $bfq1);
is $maq->run_status, 1, 'status after fasta2bfq 1 is ok';
ok -s $bfq1, 'bfq 1 file created';
$maq->fastq2bfq($read2, $bfq2);
is $maq->run_status, 1, 'status after fasta2bfq 2 is ok';
ok -s $bfq2, 'bfq 2 file created';

$maq->map($map_file, $bfa, [$bfq1, $bfq2], a => 1000, e => 120);
is $maq->run_status, 1, 'status after map is ok';
ok -s $map_file, 'mapping file created';

$maq->mapstat($map_file, $mapstat);
is $maq->run_status, 1, 'status after mapstat is ok';
my $mapped2 = 0;
ok open(my $msfh, $mapstat), 'mapstat file could be opened';
while (<$msfh>) {
    if (/Total number of reads: (\d+)/) {
        $mapped2 = $1;
        last;
    }
}
is $mapped, $mapped2, 'mapped the same number of reads both times';

# try mapping when one or both fastq contain all bad reads
unlink($map_file);
warnings_like {$maq->do_mapping(ref => $ref,
                                read1 => $bad_read1,
                                read2 => $bad_read2,
                                output => $sam_file,
                                a => 200);}
             [{carped => qr/Set clip point for .+bad_1.fastq.gz to 49\s/},
              {carped => qr/Set clip point for .+bad_2.fastq.gz to 1\s/},
              {carped => qr/fastq file .+bad_2.fastq.gz will be skipped/}],
             'do_mapping issues appropriate warnings when one of the fastq is bad';
is $maq->run_status, 1, 'status after mapping is ok';
($lines, $mapped) = check_sam($sam_file);
is $lines, 200, 'output sam not truncated';
is $mapped, 0, 'did not map any reads, since using the wrong reference';

unlink($map_file);
warnings_like {$maq->do_mapping(ref => $ref,
                                read1 => $bad_read2,
                                read2 => $bad_read2,
                                output => $sam_file,
                                a => 200);}
             [{carped => qr/Set clip point for .+bad_2.fastq.gz to 1\s/},
              {carped => qr/Set clip point for .+bad_2.fastq.gz to 1\s/},
              {carped => qr/fastq file .+bad_2.fastq.gz will be skipped/},
              {carped => qr/fastq file .+bad_2.fastq.gz will be skipped/},
              {carped => qr/none of the input fastqs had long enough sequences/}],
             'do_mapping issues appropriate warnings when both of fastqs are bad';
is $maq->run_status, 1, 'status after mapping is ok';
($lines, $mapped) = check_sam($sam_file);
is $lines, 200, 'output sam not truncated';

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
