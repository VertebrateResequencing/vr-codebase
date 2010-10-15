#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;
use File::Copy;

BEGIN {
    use Test::Most tests => 37;
    
    use_ok('VertRes::Wrapper::samtools');
    use_ok('VertRes::Utils::FileSystem');
}

my $st = VertRes::Wrapper::samtools->new(quiet => 1);
isa_ok $st, 'VertRes::Wrapper::WrapperI';
is $st->quiet, 1, 'quiet set via new';

# setup files
my $fsu = VertRes::Utils::FileSystem->new();
my $sam_input_file = File::Spec->catfile('t', 'data', 'simple.sam');
ok -s $sam_input_file, 'input sam file ready to test on';
my $ref_file = File::Spec->catfile('t', 'data', 'S_suis_P17.fa');
ok -s $ref_file, 'ref file ready to test with';
my $fai_file = File::Spec->catfile('t', 'data', 'S_suis_P17.fa.fai');
ok -s $fai_file, 'fai file ready to test with';
my $temp_dir = $fsu->tempdir();
my $bam_out_file = File::Spec->catfile($temp_dir, 'out.bam');
my $ref_file_copy = File::Spec->catfile($temp_dir, 'ref.fa');
copy($ref_file, $ref_file_copy);
ok -s $ref_file_copy, 'copy of ref file in tmp dir ready to test with';

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

# merging an RG-tagged bam with itself (simulating merging splits at the lane
# level) should generate a merged bam with a single RG tag (unlike picard-tools
# merger, which uniquifies the RG tags)
my $bam_input_file = File::Spec->catfile('t', 'data', 'rgtagged.bam');
ok -s $bam_input_file, 'input rgtagged bam file ready to test on';
$st->merge_and_check($bam_out_file, [$bam_input_file, $bam_input_file]);
cmp_ok $st->run_status, '>=', 1, 'merge_and_check ran ok';
$checkfh = $st->view($bam_out_file, undef, h => 1);
ok $checkfh, 'got a filehandle to check the bam output';
$count = 0;
my $hcount = 0;
my $rg = '';
my $saw_rgs = 0;
my %rgs;
@expected = qw(2 2 3 3 6);
while (<$checkfh>) {
    if (/^@/) {
        $hcount++;
        
        if (/^\@RG\tID:(\S+)/) {
            $rg = $1;
            $saw_rgs++;
        }
        
        next;
    }
    $count++;
    if ($count <= 5) {
        my @a = split;
        my $e = shift @expected;
        is $a[3], $e, "bam entry $count was as expected";
    }
    
    /\tRG:Z:(\S+)/;
    $rgs{$1}++;
}
is $rg, 'SRR003435', 'merge had the correct RG id in the header';
is $saw_rgs, 1, 'merge had just one RG header line';
is $hcount, 115, 'merge had the correct number of header lines';
is $count, 3770, 'merge had the correct number of entries';
is keys %rgs, 1, 'merge gave all records the same RG tag';
my ($record_rg) = keys %rgs;
is $record_rg, 'SRR003435', 'merge gave all records the correct RG tag';
my ($record_rg_count) = values %rgs;
is $record_rg_count, 3770, 'merge gave all records an RG tag';

# faidx
$st->faidx($ref_file_copy);
TODO: {
    local $TODO = "even though faidx creates a fai file, the standard check_output_files() doesn't think it exists?!";
    cmp_ok $st->run_status, '>=', 1, 'faidx ran ok';
}
ok -s $ref_file_copy.'.fai', 'faidx created a .fai file';
is_deeply [$st->faidx($ref_file, 'Streptococcus_suis:1-5', 'Streptococcus_suis:2007487-2007491')], ['atgaa', 'aaaat'], 'faidx with regions worked';
is_deeply [$st->faidx($ref_file, 'Streptococcus_suis:1-65')], ['atgaaccaagaacaacttttttggcaacgatttattgaattggcaaaggtaaattttaagccatc'], 'faidx with a region longer than 60bp worked';
is_deeply [$st->faidx($ref_file, 'Streptococcus_suis:2007492')], [], 'faidx with bad region returns undef';

# still lots more tests to do...
TODO: {
    local $TODO = "lots more tests to write...";
    ok 0;
}

exit;
