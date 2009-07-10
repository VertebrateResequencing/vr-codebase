#!/usr/bin/perl -w
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 15;
    
    use_ok('VertRes::Wrapper::samtools');
    use_ok('VertRes::IO');
}

my $st = VertRes::Wrapper::samtools->new(quiet => 1);
isa_ok $st, 'VertRes::Wrapper::WrapperI';
is $st->quiet, 1, 'quiet set via new';

# setup files
my $io = VertRes::IO->new();
my $sam_input_file = $io->catfile('t', 'data', 'simple.sam');
ok -s $sam_input_file, 'input sam file ready to test on';
my $fai_file = $io->catfile('t', 'data', 'S_suis_P17.dna.fai');
ok -s $fai_file, 'fai file ready to test with';
my $temp_dir = $io->tempdir();
my $bam_out_file = $io->catfile($temp_dir, 'out.bam');

# test the fancy multi-step methods first
$st->sam_to_fixed_sorted_bam($sam_input_file, $bam_out_file, $fai_file);
cmp_ok $st->run_status, '>=', 1, 'sam_to_fixed_sorted_bam ran ok';

$st->run_method('open');
my $checkfh = $st->view($bam_out_file);
ok $checkfh, 'got a filehandle to check the bam output';
my $count = 0;
my @expected = qw(4927 5086 9816 9816 10121);
while (<$checkfh>) {
    $count++ if /^IL/;
    if ($count <= 5) {
        my @a = split;
        my $e = shift @expected;
        is $a[3], $e, 'bam was correctly sorted';
    }
}
is $count, 2000, 'sam -> bam -> sam didn\'t lose any reads';

# still lots more tests to do...
TODO: {
    local $TODO = "lots more tests to write...";
    ok 0;
}

exit;
