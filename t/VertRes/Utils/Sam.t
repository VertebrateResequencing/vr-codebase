#!/usr/bin/perl -w
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 29;
    
    use_ok('VertRes::Utils::Sam');
    use_ok('VertRes::IO');
    use_ok('VertRes::Wrapper::samtools');
}

my $sam_util = VertRes::Utils::Sam->new();
isa_ok $sam_util, 'VertRes::Base';

# setup our input files
my $io = VertRes::IO->new();
my $bam1_file = $io->catfile('t', 'data', 'nsort1.bam');
ok -s $bam1_file, 'first bam file ready to test with';
my $bam2_file = $io->catfile('t', 'data', 'nsort2.bam');
ok -s $bam2_file, 'second bam file ready to test with';
my $headless_sam = $io->catfile('t', 'data', 'simple.sam');
ok -s $headless_sam, 'headerless sam file ready to test with';
my $ref = $io->catfile('t', 'data', 'S_suis_P17.dna');
ok -s $ref, 'ref file ready to test with';
my $dict = $io->catfile('t', 'data', 'S_suis_P17.dict');
ok -s $dict, 'dict file ready to test with';

# bams_are_similar
is $sam_util->bams_are_similar($bam1_file, $bam2_file), 1, 'bams are similar';

# add_sam_header
my $temp_dir = $io->tempdir();
my $temp_sam = $io->catfile($temp_dir, 'test.sam');
system("cp $headless_sam $temp_sam");
ok $sam_util->add_sam_header($temp_sam,
                             sample_name => 'NA00000',
                             run_name => '7563',
                             library => 'alib',
                             platform => 'SLX',
                             centre => 'Sanger',
                             insert_size => 2000,
                             lane => 'SRR00000',
                             ref_fa => $ref,
                             ref_dict => $dict,
                             ref_name => 'SsP17'), 'add_sam_header manual args test';
my @expected = ("\@HD\tVN:1.0\tSO:coordinate",
                "\@SQ\tSN:Streptococcus_suis\tLN:2007491\tAS:SsP17\tM5:c52b2f0394e12013e2573e9a38f51031\tUR:file:t/data/S_suis_P17.dna",
                "\@RG\tID:SRR00000\tPU:7563\tLB:alib\tSM:NA00000\tPI:2000\tCN:Sanger\tPL:SLX");
is_deeply [get_sam_header($temp_sam)], \@expected, 'generated the correct header';
my $wc = `wc -l $temp_sam`;
my ($lines) = $wc =~ /^(\d+)/;
is $lines, 2003, 'correct number of lines in sam after adding header';
# adding a header to a file that already has one replaces it:
ok $sam_util->add_sam_header($temp_sam,
                             sample_name => 'NA00001',
                             run_name => '7563',
                             library => 'alib',
                             platform => 'SLX',
                             centre => 'Sanger',
                             insert_size => 2000,
                             lane => 'SRR00001',
                             ref_fa => $ref,
                             ref_dict => $dict,
                             ref_name => 'SsP17'), 'add_sam_header manual args test';
@expected = ("\@HD\tVN:1.0\tSO:coordinate",
             "\@SQ\tSN:Streptococcus_suis\tLN:2007491\tAS:SsP17\tM5:c52b2f0394e12013e2573e9a38f51031\tUR:file:t/data/S_suis_P17.dna",
             "\@RG\tID:SRR00001\tPU:7563\tLB:alib\tSM:NA00001\tPI:2000\tCN:Sanger\tPL:SLX");
my @header_lines = get_sam_header($temp_sam);
is_deeply \@header_lines, \@expected, 'generated the correct header for the second time';
# we also add RG tags to the body
my @records = get_sam_body($temp_sam);
is @header_lines + @records, 2003, 'correct number of lines in sam after adding header a second time';
is substr($records[0], -14), "\tRG:Z:SRR00001", 'correct RG tag added to record';

TODO: {
    local $TODO = "Difficult to test add_sam_header with non-manual args since must fake things...";
    ok 0, 'add_sam_header auto args test';
    # ok $sam_util->add_sam_header($temp_sam,
    #                         sequence_index => 'sequence.index',
    #                         lane_path => '/path/to/lane'), 'add_sam_header auto args test';
}

# sam_to_fixed_sorted_bam (basically a shortcut to VertRes::Wrapper::samtools::sam_to_fixed_sorted_bam - no need to test thoroughly here)
my $sorted_bam = $io->catfile($temp_dir, 'sorted.bam');
ok $sam_util->sam_to_fixed_sorted_bam($temp_sam, $sorted_bam, $ref, quiet => 1), 'sam_to_fixed_sorted_bam test';
@records = get_bam_body($sorted_bam);
is @records, 2000, 'sorted bam had the correct number of records';
like $records[0], qr/\tNM:i:0\tMD:Z:54/, 'record with no prior NM/MD tags now has the correct NM and MD flags';
unlike $records[2], qr/\tNM:i:\d.+\tNM:i:\d/, 'record with existing NM tag didn\'t duplicate the NM flags';

# rmdup (just a shortcut to VertRes::Wrapper::samtools::rmdup - no need to test thoroughly here)
my $rmdup_bam = $io->catfile($temp_dir, 'rmdup.bam');
ok $sam_util->rmdup($sorted_bam, $rmdup_bam, single_ended => 0, quiet => 1), 'rmdup test, paired';
ok $sam_util->rmdup($sorted_bam, $rmdup_bam, single_ended => 1, quiet => 1), 'rmdup test, single ended';

# merge (pretty much just a shortcut to VertRes::Wrapper::picard::merge_and_check - no need to test thoroughly here)
my $merge_bam = $io->catfile($temp_dir, 'merge.bam');
ok $sam_util->merge($merge_bam, $rmdup_bam), 'merge on single bam test';
ok -l $merge_bam, 'merge on single bam created a symlink';
ok $sam_util->merge($merge_bam, $rmdup_bam, $sorted_bam), 'merge on multiple bam test';
ok -f $merge_bam, 'merge on multiple bams created a new file';
ok -s $merge_bam > -s $rmdup_bam, 'merge on multiple bams created a new file bigger than one of the originals';

exit;

sub get_sam_header {
    my $sam_file = shift;
    open(my $samfh, $sam_file) || die "Could not open sam file '$sam_file'\n";
    my @header_lines;
    while (<$samfh>) {
        chomp;
        last unless /^@/;
        push(@header_lines, $_);
    }
    return @header_lines;
}

sub get_sam_body {
    my $sam_file = shift;
    open(my $samfh, $sam_file) || die "Could not open sam file '$sam_file'\n";
    my @records;
    while (<$samfh>) {
        chomp;
        next if /^@/;
        push(@records, $_);
    }
    return @records;
}

sub get_bam_body {
    my $bam_file = shift;
    my $samtools = VertRes::Wrapper::samtools->new(quiet => 1);
    $samtools->run_method('open');
    my $bamfh = $samtools->view($bam_file);
    my @records;
    while (<$bamfh>) {
        chomp;
        next if /^@/;
        push(@records, $_);
    }
    return @records;
}
