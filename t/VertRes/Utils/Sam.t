#!/usr/bin/perl -w
use strict;
use warnings;
use File::Spec;

BEGIN {
    use Test::Most tests => 88;
    
    use_ok('VertRes::Utils::Sam');
    use_ok('VertRes::Wrapper::samtools');
    use_ok('VertRes::Utils::FileSystem');
}

my $sam_util = VertRes::Utils::Sam->new(java_memory => 1000);
isa_ok $sam_util, 'VertRes::Base';

# setup our input files
my $bam1_file = File::Spec->catfile('t', 'data', 'nsort1.bam');
ok -s $bam1_file, 'first bam file ready to test with';
my $bam2_file = File::Spec->catfile('t', 'data', 'nsort2.bam');
ok -s $bam2_file, 'second bam file ready to test with';
my $headed_bam = File::Spec->catfile('t', 'data', 'headed2.bam');
ok -s $headed_bam, 'headed bam file ready to test with';
my $headless_sam = File::Spec->catfile('t', 'data', 'simple.sam');
ok -s $headless_sam, 'headerless sam file ready to test with';
my $ref = File::Spec->catfile('t', 'data', 'S_suis_P17.fa');
ok -s $ref, 'ref file ready to test with';
my $dict = File::Spec->catfile('t', 'data', 'S_suis_P17.dict');
ok -s $dict, 'dict file ready to test with';

# bams_are_similar
is $sam_util->bams_are_similar($bam1_file, $bam2_file), 1, 'bams are similar';

# add_sam_header
my $fsu = VertRes::Utils::FileSystem->new();
my $temp_dir = $fsu->tempdir();
my $temp_sam = File::Spec->catfile($temp_dir, 'test.sam');
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
                             ref_name => 'SsP17',
                             project => 'SRP000001'), 'add_sam_header manual args test';
my @expected = ("\@HD\tVN:1.0\tSO:coordinate",
                "\@SQ\tSN:Streptococcus_suis\tLN:2007491\tAS:SsP17\tM5:c52b2f0394e12013e2573e9a38f51031\tUR:file:t/data/S_suis_P17.fa",
                "\@RG\tID:SRR00000\tLB:alib\tSM:NA00000\tPU:7563\tPI:2000\tCN:Sanger\tPL:ILLUMINA\tDS:SRP000001");
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
                             ref_name => 'SsP17',
                             project => 'SRP000001',
                             program => 'bwa',
                             program_version => '0.4.9'), 'add_sam_header manual args test';
@expected = ("\@HD\tVN:1.0\tSO:coordinate",
             "\@SQ\tSN:Streptococcus_suis\tLN:2007491\tAS:SsP17\tM5:c52b2f0394e12013e2573e9a38f51031\tUR:file:t/data/S_suis_P17.fa",
             "\@RG\tID:SRR00001\tLB:alib\tSM:NA00001\tPU:7563\tPI:2000\tCN:Sanger\tPL:ILLUMINA\tDS:SRP000001",
             "\@PG\tID:bwa\tVN:0.4.9");
my @header_lines = get_sam_header($temp_sam);
is_deeply \@header_lines, \@expected, 'generated the correct header for the second time';
# we also add RG tags to the body
my @records = get_sam_body($temp_sam);
is @header_lines + @records, 2004, 'correct number of lines in sam after adding header a second time';
my $found_rgs = 0;
foreach my $record (@records) {
    $found_rgs++ if substr($record, -14) eq "\tRG:Z:SRR00001";
}
is $found_rgs, @records, 'correct RG tag added to all records';

TODO: {
    local $TODO = "Difficult to test add_sam_header with non-manual args since must fake things...";
    ok 0, 'add_sam_header auto args test';
    # ok $sam_util->add_sam_header($temp_sam,
    #                         sequence_index => 'sequence.index',
    #                         lane_path => '/path/to/lane'), 'add_sam_header auto args test';
}

# sam_to_fixed_sorted_bam (basically a shortcut to VertRes::Wrapper::samtools::sam_to_fixed_sorted_bam - no need to test thoroughly here)
my $sorted_bam = File::Spec->catfile($temp_dir, 'sorted.bam');
ok $sam_util->sam_to_fixed_sorted_bam($temp_sam, $sorted_bam, $ref, quiet => 1), 'sam_to_fixed_sorted_bam test';
@records = get_bam_body($sorted_bam);
is @records, 2000, 'sorted bam had the correct number of records';
like $records[0], qr/\tNM:i:0\tMD:Z:54/, 'record with no prior NM/MD tags now has the correct NM and MD flags';
unlike $records[2], qr/\tNM:i:\d.+\tNM:i:\d/, 'record with existing NM tag didn\'t duplicate the NM flags';
# (a samtools fillmd bug meant that sometimes the RG tag would get removed)
$found_rgs = 0;
foreach my $record (@records) {
    $found_rgs++ if $record =~ /\tRG:Z:SRR00001/;
}
is $found_rgs, @records, 'correct RG tag still present on all records';

# rewrite_bam_header
ok $sam_util->rewrite_bam_header($sorted_bam, invalid => { sample_name => 'NA00002', library => 'blib', centre => 'NCBI' }), 'rewrite_bam_header ran ok with an invalid readgroup';
@header_lines = get_bam_header($sorted_bam);
is $header_lines[2], "\@RG\tID:SRR00001\tLB:alib\tSM:NA00001\tPU:7563\tPI:2000\tCN:Sanger\tPL:ILLUMINA\tDS:SRP000001", 'rewrite_bam_header didn\'t change the header when readgroup not in the bam';
my %new_header_args = (SRR00001 => { sample_name => 'NA00002', library => 'blib', centre => 'NCBI', project => 'SRP000001' });
is $sam_util->check_bam_header($sorted_bam, %new_header_args), 1, 'before rewriting header, check returns true';
ok $sam_util->rewrite_bam_header($sorted_bam, %new_header_args), 'rewrite_bam_header ran ok';
@header_lines = get_bam_header($sorted_bam);
is $header_lines[2], "\@RG\tID:SRR00001\tLB:blib\tSM:NA00002\tPU:7563\tPI:2000\tCN:NCBI\tPL:ILLUMINA\tDS:SRP000001", 'rewrite_bam_header actually changed the header';
@records = get_bam_body($sorted_bam);
is @records, 2000, 'rewrite_bam_header didn\'t change the number of records';
is $sam_util->check_bam_header($sorted_bam, %new_header_args), 0, 'after rewriting header, check returns false';

# rmdup (just a shortcut to VertRes::Wrapper::samtools::rmdup - no need to test thoroughly here)
my $rmdup_bam = File::Spec->catfile($temp_dir, 'rmdup.bam');
ok $sam_util->rmdup($sorted_bam, $rmdup_bam, single_ended => 0, quiet => 1), 'rmdup test, paired';
my @rmdup_pe = get_bam_body($rmdup_bam);
unlink($rmdup_bam);
ok $sam_util->rmdup($sorted_bam, $rmdup_bam, single_ended => 1, quiet => 1), 'rmdup test, single ended';
my @rmdup_se = get_bam_body($rmdup_bam);
cmp_ok scalar(@rmdup_pe), '<', scalar(@rmdup_se), 'rmdup pe gave fewer reads than se mode';

# markdup (just a shortcut to VertRes::Wrapper::picard::markdup - no need to test thoroughly here)
ok $sam_util->markdup($sorted_bam, $rmdup_bam), 'markdup test';
@rmdup_pe = get_bam_body($rmdup_bam);
is @rmdup_pe, 2000, 'markdup gave the correct number of reads';

# merge (pretty much just a shortcut to VertRes::Wrapper::picard::merge_and_check - no need to test thoroughly here)
my $merge_bam = File::Spec->catfile($temp_dir, 'merge.bam');
ok $sam_util->merge($merge_bam, $rmdup_bam), 'merge on single bam test';
is readlink($merge_bam), 'rmdup.bam', 'merge on single bam created the correct relative symlink';
ok $sam_util->merge($merge_bam, $rmdup_bam, $sorted_bam), 'merge on multiple bam test';
ok -f $merge_bam, 'merge on multiple bams created a new file';
ok -s $merge_bam > -s $rmdup_bam, 'merge on multiple bams created a new file bigger than one of the originals';

# calculate_flag
is $sam_util->calculate_flag(), 0, 'calculate_flag no args test';
is $sam_util->calculate_flag(paired_tech => 1), 1, 'calculate_flag paired_tech test';
warning_like {is $sam_util->calculate_flag(self_unmapped => 1, paired_map => 1), 4, 'calculate_flag self_unmapped test';} {carped => qr/was set/}, 'calculate_flag extra paired_map warning test';
warning_like {is $sam_util->calculate_flag(mate_unmapped => 1), 9, 'calculate_flag mate_unmapped test';} {carped => qr/forcing paired_tech on/}, 'calculate_flag missing paired_tech warning test';
warning_like {is $sam_util->calculate_flag('1st_in_pair' => 1, '2nd_in_pair' => 1), 0, 'calculate_flag both test';} {carped => qr/forcing both off/}, 'calculate_flag both warning test';
is $sam_util->calculate_flag(self_reverse => 1), 16, 'calculate_flag self_reverse test';
is $sam_util->calculate_flag(paired_tech => 1, mate_reverse => 1), 33, 'calculate_flag mate_reverse test';
is $sam_util->calculate_flag(not_primary => 1), 256, 'calculate_flag not_primary test';
is $sam_util->calculate_flag(failed_qc => 1), 512, 'calculate_flag failed_qc test';
is $sam_util->calculate_flag(duplicate => 1), 1024, 'calculate_flag duplicate test';
is $sam_util->calculate_flag(paired_tech => 1, paired_map => 1, mate_reverse => 1, '1st_in_pair' => 1), 99, 'calculate_flag mapped pair first test';
is $sam_util->calculate_flag(paired_tech => 1, paired_map => 1, self_reverse => 1, '2nd_in_pair' => 1), 147, 'calculate_flag mapped pair second test';
is $sam_util->calculate_flag(paired_tech => 1, self_unmapped => 1, mate_reverse => 1, '1st_in_pair' => 1), 101, 'calculate_flag pair self unmapped test';
is $sam_util->calculate_flag(paired_tech => 1, mate_unmapped => 1, self_reverse => 1, '2nd_in_pair' => 1), 153, 'calculate_flag pair mate unmapped args test';
is $sam_util->calculate_flag(paired_tech => 1, self_unmapped => 1, '1st_in_pair' => 1), 69, 'calculate_flag pair both unmapped first test';
is $sam_util->calculate_flag(paired_tech => 1, self_unmapped => 1, '2nd_in_pair' => 1), 133, 'calculate_flag pair both unmapped second test';
is $sam_util->calculate_flag(paired_tech => 1, self_unmapped => 1, mate_unmapped => 1, '1st_in_pair' => 1), 77, 'calculate_flag both unmapped v2 first test';
is $sam_util->calculate_flag(paired_tech => 1, self_unmapped => 1, mate_unmapped => 1, '2nd_in_pair' => 1), 141, 'calculate_flag both unmapped v2 second test';

# bam_statistics and bas
# stats independently verified with samtools flagstat,
# Math::NumberCruncher and Statistics::Robust::Scale:
# perl -MVertRes::Utils::FastQ -MMath::NumberCruncher -MVertRes::Parser::sam -Mstrict -we 'my $fqu = VertRes::Utils::FastQ->new(); my $pars = VertRes::Parser::sam->new(file => "t/data/simple.sam"); my $rh = $pars->result_holder; my %stats; while ($pars->next_result) { my $flag = $rh->{FLAG}; if ($pars->is_mapped($flag)) { my $seq = $rh->{SEQ}; $stats{mapped_bases} += length($seq); $stats{mapped_reads}++ if $pars->is_sequencing_paired($flag); foreach my $qual ($fqu->qual_to_ints($rh->{QUAL})) { push(@{$stats{qs}}, $qual); } if ($rh->{MAPQ} > 0) { push(@{$stats{isizes}}, $rh->{ISIZE}) if $rh->{ISIZE} > 0; } } } while (my ($stat, $val) = each %stats) { print "$stat => $val\n"; } print "avg isize: ", Math::NumberCruncher::Mean($stats{isizes}), "\n"; print "sd isize: ", Math::NumberCruncher::StandardDeviation($stats{isizes}), "\n"; print "med isize: ", Math::NumberCruncher::Median($stats{isizes}), "\n"; use Statistics::Robust::Scale "MAD"; print "mad isize: ", MAD($stats{isizes}), "\n"; print "avg qual: ", Math::NumberCruncher::Mean($stats{qs}), "\n";'
# , and percent_mismatch with:
# perl -e '$bases = 0; $matches = 0; open($fh, "samtools pileup -sf t/data/S_suis_P17.fa sorted.bam |"); while (<$fh>) { @s = split; $bases += $s[3]; $matches += $s[4] =~ tr/.,/.,/; } print "$bases / $matches\n"; $p = 100 - ((100/$bases) * $matches); print "percent mismatches: $p\n";'
is_deeply {$sam_util->bam_statistics($sorted_bam)}, {SRR00001 => {total_bases => 115000,
                                                                  mapped_bases => 58583,
                                                                  total_reads => 2000,
                                                                  mapped_reads => 1084,
                                                                  mapped_reads_paired_in_seq => 1084,
                                                                  mapped_reads_properly_paired => 1070,
                                                                  percent_mismatch => '2.05',
                                                                  avg_qual => '23.32',
                                                                  avg_isize => 286,
                                                                  sd_isize => '74.10',
                                                                  median_isize => 275,
                                                                  mad => 48}}, 'bam_statistics test';
my $given_bas = File::Spec->catfile($temp_dir, 'test.bas');
ok $sam_util->bas($sorted_bam, '20100208', $given_bas), 'bas() ran ok';
my $expected_bas = File::Spec->catfile('t', 'data', 'example.bas');
ok open(my $ebfh, $expected_bas), 'opened expected .bas';
@expected = <$ebfh>;
ok open(my $tbfh, $given_bas), 'opened result .bas';
my @given = <$tbfh>;
close($tbfh);
is_deeply \@given, \@expected, 'bas output was as expected';

# rewrite_bas_meta
ok $sam_util->rewrite_bas_meta($given_bas, invalid => { sample_name => 'NA00003', library => 'clib', centre => 'FOO' }), 'rewrite_bas_meta ran ok with an invalid readgroup';
open($tbfh, $given_bas);
@given = <$tbfh>;
close($tbfh);
is $given[1], $expected[1], 'rewrite_bas_meta didn\'t change anything when readgroup not in the bas';
ok $sam_util->rewrite_bas_meta($given_bas, SRR00001 => { sample_name => 'NA00003', library => 'clib', platform => 'LS454', project => 'SRP000002' }), 'rewrite_bas_meta ran ok';
open($tbfh, $given_bas);
@given = <$tbfh>;
close($tbfh);
is_deeply \@given, [$expected[0], join("\t", qw(NA00003.LS454.bwa.SRP000002.20100208 d4f062c6eed8d0205517979801305bfd SRP000002 NA00003 LS454 clib SRR00001 115000 58583 2000 1084 1084 1070 2.05 23.32 286 74.10 275 48))."\n"], 'rewrite_bas_meta actually changed the bas';

# stats method
ok $sam_util->stats('20100208', $sorted_bam), 'stats test';
foreach my $file (File::Spec->catfile($temp_dir, 'sorted.bam.flagstat'), File::Spec->catfile($temp_dir, 'sorted.bam.bas')) {
    ok -s $file, 'stats output file exists';
}

# split_bam_by_sequence method
my @splits = $sam_util->split_bam_by_sequence($headed_bam,
                                              output_dir => $temp_dir,
                                              ignore => '^N[TC]_\d+',
                                              make_unmapped => 0,
                                              non_chr => 0);
my @expected_splits;
foreach (1..22, 'X', 'MT') {
    push(@expected_splits, File::Spec->catfile($temp_dir, 'chrom'.$_.'.headed2.bam'));
}
@expected_splits = sort @expected_splits;
@splits = sort @splits;
is_deeply \@splits, \@expected_splits, 'split_bam_by_sequence claimed to make the appropriate files';
my $actually_created = 0;
my $created_lines = 0;
foreach my $split (@expected_splits) {
    $actually_created += -s $split ? 1 : 0;
    $created_lines += get_bam_body($split);
}
is $actually_created, 24, 'split_bam_by_sequence actually created all the split bams';
is $created_lines, 30, 'split_bam_by_sequence created split bams with appropriate numbers of entries';
@splits = $sam_util->split_bam_by_sequence($headed_bam,
                                           output_dir => $temp_dir,
                                           ignore => '^N[TC]_\d+',
                                           make_unmapped => 1,
                                           non_chr => 0);
push(@expected_splits, File::Spec->catfile($temp_dir, 'unmapped.headed2.bam'));
$actually_created = 0;
$created_lines = 0;
foreach my $split (@expected_splits) {
    $actually_created += -s $split ? 1 : 0;
    $created_lines += get_bam_body($split);
}
is $actually_created, 25, 'split_bam_by_sequence with make_unmapped gives an unmapped bam as well';
is $created_lines, 1998, 'and the unmapped gives us more entries';
@splits = $sam_util->split_bam_by_sequence($headed_bam,
                                           output_dir => $temp_dir);
push(@expected_splits, File::Spec->catfile($temp_dir, 'nonchrom.headed2.bam'));
$actually_created = 0;
$created_lines = 0;
foreach my $split (@expected_splits) {
    $actually_created += -s $split ? 1 : 0;
    $created_lines += get_bam_body($split);
    unlink($split);
}
is $actually_created, 26, 'split_bam_by_sequence defaults gives a nonchrom bam as well';
is $created_lines, 2000, 'and the nonchrom gives us all entries';
@splits = $sam_util->split_bam_by_sequence($headed_bam,
                                           output_dir => $temp_dir,
                                           pretend => 1);
$actually_created = 0;
foreach my $split (@expected_splits) {
    $actually_created += -s $split ? 1 : 0;
}
is @splits, 26, 'split_bam_by_sequence can pretend to make all splits';
is $actually_created, 0, 'split_bam_by_sequence pretend doesn\'t actually make any splits';

# add_unmapped
my $om_sam_orig = File::Spec->catfile('t', 'data', 'only_mapped.sam');
my $om_sam = File::Spec->catfile($temp_dir, 'only_mapped.sam');
system("cp $om_sam_orig $om_sam");
@records = get_sam_body($om_sam);
is @records, 6, 'only_mapped.sam starts with 6 records';
my $om_1_fq = File::Spec->catfile('t', 'data', 'only_mapped_1.fastq');
ok -s $om_1_fq, 'only mapped fq1 ready to test with';
my $om_2_fq = File::Spec->catfile('t', 'data', 'only_mapped_2.fastq');
ok -s $om_2_fq, 'only mapped fq2 ready to test with';
ok $sam_util->add_unmapped($om_sam, $om_1_fq, $om_2_fq), 'add_unmapped returned true';
@records = get_sam_body($om_sam);
is @records, 10, 'after adding back unmapped we have 10 records';

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

sub get_bam_header {
    my $bam_file = shift;
    my $samtools = VertRes::Wrapper::samtools->new(quiet => 1);
    $samtools->run_method('open');
    my $bamfh = $samtools->view($bam_file, undef, H => 1);
    my @header_lines;
    while (<$bamfh>) {
        chomp;
        push(@header_lines, $_);
    }
    return @header_lines;
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
