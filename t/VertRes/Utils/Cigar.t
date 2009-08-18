#!/usr/bin/perl -w
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 33;
    
    use_ok('VertRes::Utils::Cigar');
    use_ok('VertRes::IO');
    use_ok('VertRes::Parser::sam');
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

is $cigar_util->cigar_to_sam([$fq1_file, $fq2_file], [$cigar1_file, $cigar2_file], 2000, undef, $sam_out), 2000, 'cigar_to_sam claimed to write the correct number of entries';
my ($mapped, $unmapped) = sam_lines_check($sam_out);
is_deeply [$mapped, $unmapped], [1635, 365], 'cigar_to_sam generated sam with correct number of mapped and unmapped reads';

is $cigar_util->cigar_to_sam([$fq1_file], [$cigar1_file], 2000, undef, $sam_out), 1000, 'cigar_to_sam claimed to write the correct number of entries - first only';
($mapped, $unmapped) = sam_lines_check($sam_out);
is_deeply [$mapped, $unmapped], [878, 122], 'cigar_to_sam generated sam with correct number of mapped and unmapped reads - first only';

is $cigar_util->cigar_to_sam([$fq2_file], [$cigar2_file], 2000, undef, $sam_out), 1000, 'cigar_to_sam claimed to write the correct number of entries - second only';
($mapped, $unmapped) = sam_lines_check($sam_out);
is_deeply [$mapped, $unmapped], [757, 243], 'cigar_to_sam generated sam with correct number of mapped and unmapped reads - second only';


# above tests are with non-454 fastqs. These fastqs are real G1K 454 reads and
# caused problems during the mapping pipeline:
$fq1_file = $io->catfile('t', 'data', 'SRR001629_1.fastq');
ok -s $fq1_file, 'SRR001629_1 fastq file ready to test with';
$fq2_file = $io->catfile('t', 'data', 'SRR001629_2.fastq');
ok -s $fq2_file, 'SRR001629_2 fastq file ready to test with';
$cigar1_file = $io->catfile('t', 'data', 'SRR001629_1.cigar');
ok -s $cigar1_file, 'SRR001629_1 cigar file ready to test with';
$cigar2_file = $io->catfile('t', 'data', 'SRR001629_2.cigar');
ok -s $cigar2_file, 'SRR001629_2 cigar file ready to test with';
is $cigar_util->cigar_to_sam([$fq1_file, $fq2_file], [$cigar1_file, $cigar2_file], 2000, undef, $sam_out), 2000, 'cigar_to_sam worked with unfiltered SRR001629 reads';
($mapped, $unmapped) = sam_lines_check($sam_out);
is_deeply [$mapped, $unmapped], [1736, 264], 'cigar_to_sam generated sam with correct number of mapped and unmapped reads - unfiltered SRR001629';

$fq1_file = $io->catfile('t', 'data', 'SRR001629_1.filtered.fastq');
ok -s $fq1_file, 'SRR001629_1 filtered fastq file ready to test with';
$fq2_file = $io->catfile('t', 'data', 'SRR001629_2.filtered.fastq');
ok -s $fq2_file, 'SRR001629_2 filtered fastq file ready to test with';
is $cigar_util->cigar_to_sam([$fq1_file, $fq2_file], [$cigar1_file, $cigar2_file], 2000, undef, $sam_out), 1770, 'cigar_to_sam worked with filtered SRR001629 reads';
($mapped, $unmapped) = sam_lines_check($sam_out);
is_deeply [$mapped, $unmapped], [1736, 34], 'cigar_to_sam generated sam with correct number of mapped and unmapped reads - filtered SRR001629';

# for this read pair it was outputting read _1 in the sam instead of read _2
$fq1_file = $io->catfile('t', 'data', 'SRR001629.99999_1.fastq');
ok -s $fq1_file, 'SRR001629.99999_1 filtered fastq file ready to test with';
$fq2_file = $io->catfile('t', 'data', 'SRR001629.99999_2.fastq');
ok -s $fq2_file, 'SRR001629.99999_2 filtered fastq file ready to test with';
$cigar1_file = $io->catfile('t', 'data', 'SRR001629.99999_1.cigar');
ok -s $cigar1_file, 'SRR001629.99999_1.cigar cigar file ready to test with';
$cigar2_file = $io->catfile('t', 'data', 'SRR001629.99999_2.cigar');
ok -s $cigar2_file, 'SRR001629.99999_1.cigar cigar file ready to test with';
is $cigar_util->cigar_to_sam([$fq1_file, $fq2_file], [$cigar1_file, $cigar2_file], 2000, undef, $sam_out), 6, 'cigar_to_sam worked with mini .99999 cigars';
$io->file($sam_out);
my $fh = $io->fh;
<$fh>;
<$fh>;
<$fh>;
my $line = <$fh>;
# the 137 flag is important: in the old mapping pipline code, this came out as
# the invalid flag 139 which means 'proper pair' yet also 'mate unmapped'.
# 137 doesn't claim to be properly paired... scratch that, now that we futher
# updated to make repetative reads mapping q 0, we have our mate after all
# and get flag 131.
is $line, "SRR001629.99999\t131\t6\t123545210\t254\t62M\t=\t123547570\t-2397\tTTTCCCAGTCAAGGGAAGCAGTAAGGGATTGTGCTATCCGGCCCAGGTACTATGCTTTTCCC\t".'>998889;@@;;;;>;;9>><<<<<<B<<<>@>>>>@@@99999999>>@@>>@@9988888'."\n", 'the sam line output for the .99999 cigar was correct';

# could do with complete tests to check the flags output in all scenarios...
TODO: {
    local $TODO = "Need to test flag output of cigar_to_sam in more scnearios";
    ok 0, 'cigar_to_sam flag tests';
}

exit;

sub sam_lines_check {
    my $sam = shift;
    my $parser = VertRes::Parser::sam->new(file => $sam);
    my $rh = $parser->result_holder;
    
    my ($mapped, $unmapped) = (0, 0);
    while ($parser->next_result) {
        my $flag = $rh->{FLAG};
        if ($parser->is_mapped($flag)) {
            $mapped++;
        }
        else {
            $unmapped++;
        }
    }
    return ($mapped, $unmapped);
}
