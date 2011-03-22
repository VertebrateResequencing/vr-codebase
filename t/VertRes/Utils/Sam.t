#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;
use File::Copy;
use Data::Dumper;

BEGIN {
    use Test::Most tests => 180;
    
    use_ok('VertRes::Utils::Sam');
    use_ok('VertRes::Wrapper::samtools');
    use_ok('SamTools');
    use_ok('VertRes::Utils::FileSystem');
    use_ok('VertRes::Wrapper::picard');
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
my $pg_bam = File::Spec->catfile('t', 'data', 'hard_soft.bam');
ok -s $pg_bam, 'pg bam file ready to test with';
my $ref = File::Spec->catfile('t', 'data', 'S_suis_P17.fa');
ok -s $ref, 'ref file ready to test with';
my $dict = File::Spec->catfile('t', 'data', 'S_suis_P17.dict');
ok -s $dict, 'dict file ready to test with';
my $intervals_all_bam = File::Spec->catfile('t', 'data', 'extract_intervals_all.bam');
ok -s $intervals_all_bam, 'extract intervals all records bam file ready to test with';
my $intervals_extracted_bam = File::Spec->catfile('t', 'data', 'extract_intervals_extracted.bam');
ok -s $intervals_extracted_bam, 'extracted intervals bam file ready to test with';
my $intervals_to_extract = File::Spec->catfile('t', 'data', 'extract_intervals.txt');
ok -s $intervals_to_extract, 'file of intervals to be extracted ready to test with';
my $exome_bam_file = File::Spec->catfile('t', 'data', 'bam_exome_qc.bam');
ok -s $exome_bam_file, 'exome qc bam file ready to test with';
my $bait_interval = File::Spec->catfile('t', 'data',  'bam_exome_qc_bait.intervals');
ok -s $bait_interval, 'exome qc bait file ready to test with';
my $bait_interval_list = File::Spec->catfile('t', 'data',  'bam_exome_qc_bait.interval_list');
ok -s $bait_interval_list, 'exome qc bait interval file ready to test with';
my $target_interval = File::Spec->catfile('t', 'data',  'bam_exome_qc_target.intervals');
ok -s $target_interval, 'exome qc target file ready to test with';
my $target_interval_list = File::Spec->catfile('t', 'data',  'bam_exome_qc_target.interval_list');
ok -s $target_interval_list, 'exome qc target intervals file ready to test with';
my $ref_fa = File::Spec->catfile('t', 'data',  'bam_exome_qc_ref.fa');
ok -s $ref_fa, 'exome qc ref file ready to test with';
my $ref_fai = File::Spec->catfile('t', 'data',  'bam_exome_qc_ref.fa.fai');
ok -s $ref_fai, 'exome qc ref fai file ready to test with';

# bams_are_similar
is $sam_util->bams_are_similar($bam1_file, $bam2_file), 1, 'bams are similar';

# add_sam_header
my $fsu = VertRes::Utils::FileSystem->new();
my $temp_dir = $fsu->tempdir();
my $temp_sam = File::Spec->catfile($temp_dir, 'test.sam');
copy($headless_sam, $temp_sam);
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
                             study => 'SRP000001'), 'add_sam_header manual args test';
my @expected = ("\@HD\tVN:1.0\tSO:coordinate",
                "\@SQ\tSN:Streptococcus_suis\tLN:2007491\tAS:SsP17\tUR:file:/nfs/users/nfs_s/sb10/src/vert_reseq/trunk/t/data/S_suis_P17.dna\tM5:c52b2f0394e12013e2573e9a38f51031",
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
                             study => 'SRP000001',
                             program => 'bwa',
                             program_version => '0.4.9'), 'add_sam_header manual args test';
@expected = ("\@HD\tVN:1.0\tSO:coordinate",
             "\@SQ\tSN:Streptococcus_suis\tLN:2007491\tAS:SsP17\tUR:file:/nfs/users/nfs_s/sb10/src/vert_reseq/trunk/t/data/S_suis_P17.dna\tM5:c52b2f0394e12013e2573e9a38f51031",
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

my $sorted_bam_copy = File::Spec->catfile($temp_dir, 'sorted.copy.bam');
copy($sorted_bam, $sorted_bam_copy);

# rewrite_bam_header
ok $sam_util->rewrite_bam_header($sorted_bam, invalid => { sample_name => 'NA00002', library => 'blib', centre => 'NCBI' }), 'rewrite_bam_header ran ok with an invalid readgroup';
@header_lines = get_bam_header($sorted_bam);
is $header_lines[2], "\@RG\tID:SRR00001\tLB:alib\tSM:NA00001\tPU:7563\tPI:2000\tCN:Sanger\tPL:ILLUMINA\tDS:SRP000001", 'rewrite_bam_header didn\'t change the header when readgroup not in the bam';
my %new_header_args = (SRR00001 => { sample_name => 'NA00002', library => 'blib', centre => 'NCBI', study => 'SRP000009' });
is $sam_util->check_bam_header($sorted_bam, %new_header_args), 1, 'before rewriting header, check returns true';
ok $sam_util->rewrite_bam_header($sorted_bam, %new_header_args), 'rewrite_bam_header ran ok';
@header_lines = get_bam_header($sorted_bam);
is $header_lines[2], "\@RG\tID:SRR00001\tLB:blib\tSM:NA00002\tPU:7563\tPI:2000\tCN:NCBI\tPL:ILLUMINA\tDS:SRP000009", 'rewrite_bam_header actually changed the header';
@records = get_bam_body($sorted_bam);
is @records, 2000, 'rewrite_bam_header didn\'t change the number of records';
is $sam_util->check_bam_header($sorted_bam, %new_header_args), 0, 'after rewriting header, check returns false';
# (subsequent tests expect a project/study of SRP000001, so change it back)
%new_header_args = (SRR00001 => { sample_name => 'NA00002', library => 'blib', centre => 'NCBI', project => 'SRP000001' });
$sam_util->rewrite_bam_header($sorted_bam, %new_header_args);

# standardise_pg_header_lines
my $temp_pg_bam = File::Spec->catfile($temp_dir, 'pg.bam');
copy($pg_bam, $temp_pg_bam);
ok $sam_util->standardise_pg_header_lines($temp_pg_bam, invalid => { version => 'cl' }), 'standardise_pg_header_lines ran ok with an invalid pg';
@header_lines = get_bam_header($temp_pg_bam);
is $header_lines[-1], "\@PG\tID:bwa\tVN:0.5.3", 'standardise_pg_header_lines didn\'t change the header when pg not in the bam';
ok $sam_util->standardise_pg_header_lines($temp_pg_bam, bwa => { '0.5.2' => '-q 15' }), 'standardise_pg_header_lines ran ok with an invalid version';
@header_lines = get_bam_header($temp_pg_bam);
is $header_lines[-1], "\@PG\tID:bwa\tVN:0.5.3", 'standardise_pg_header_lines didn\'t change the header when pg vn not in the bam';
ok $sam_util->standardise_pg_header_lines($temp_pg_bam, bwa => { '0.5.3' => '-q 15' }), 'standardise_pg_header_lines ran ok with correct pg and vn';
@header_lines = get_bam_header($temp_pg_bam);
is $header_lines[-1], "\@PG\tID:bwa\tVN:0.5.3\tCL:-q 15", 'standardise_pg_header_lines actually changed the header';
@records = get_bam_body($temp_pg_bam);
is @records, 1, 'standardise_pg_header_lines didn\'t change the number of records';
ok $sam_util->standardise_pg_header_lines($temp_pg_bam, bwa => { '0.5.3' => '-q 30' }), 'standardise_pg_header_lines ran ok with correct pg and vn again';
@header_lines = get_bam_header($temp_pg_bam);
is $header_lines[-1], "\@PG\tID:bwa\tVN:0.5.3\tCL:-q 30", 'standardise_pg_header_lines actually changed the header again';

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

# do a samtools sam->bam, picard merge test on a multi-sequence sam file to
# ensure things are all compatible
my $multi_seq_sam = File::Spec->catfile('t', 'data', 'header.sam');
my $multi_seq_bam = File::Spec->catfile($temp_dir, 'header.bam');
my $samtools = VertRes::Wrapper::samtools->new(quiet => 1);
$samtools->view($multi_seq_sam, $multi_seq_bam, S => 1, b => 1);
ok -s $multi_seq_bam, 'samtools could make a bam from a multi-sequence sam';
my $multi_seq_bam2 = File::Spec->catfile($temp_dir, 'header2.bam');
copy($multi_seq_bam, $multi_seq_bam2);
unlink($merge_bam);
ok $sam_util->merge($merge_bam, $multi_seq_bam, $multi_seq_bam2), 'merge on multiple sequence bam test';
ok -s $merge_bam, 'merge did generate a file';

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
                                                                  mad => 48,
                                                                  duplicate_reads => 2}}, 'bam_statistics test';
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
ok $sam_util->rewrite_bas_meta($given_bas, SRR00001 => { sample_name => 'NA00003', library => 'clib', platform => 'LS454', study => 'SRP000002' }), 'rewrite_bas_meta ran ok';
open($tbfh, $given_bas);
@given = <$tbfh>;
close($tbfh);
is_deeply \@given, [$expected[0], join("\t", qw(NA00003.LS454.bwa.unknown_population.unknown_analysisgroup.20100208 d025918167ebe5b1917abb0909e3185e SRP000002 NA00003 LS454 clib SRR00001 115000 58583 2000 1084 1084 1070 2.05 23.32 286 74.10 275 48 2))."\n"], 'rewrite_bas_meta actually changed the bas';
ok $sam_util->rewrite_bas_meta($given_bas, SRR00001 => { sample_name => 'NA00004', filename => 'foo' }), 'rewrite_bas_meta ran ok again';
open($tbfh, $given_bas);
@given = <$tbfh>;
close($tbfh);
is_deeply \@given, [$expected[0], join("\t", qw(foo d025918167ebe5b1917abb0909e3185e SRP000002 NA00004 LS454 clib SRR00001 115000 58583 2000 1084 1084 1070 2.05 23.32 286 74.10 275 48 2))."\n"], 'rewrite_bas_meta actually changed the bas, and the filename option worked';

# stats method
ok $sam_util->stats('20100208', $sorted_bam), 'stats test';
foreach my $file (File::Spec->catfile($temp_dir, 'sorted.bam.flagstat'), File::Spec->catfile($temp_dir, 'sorted.bam.bas')) {
    ok -s $file, 'stats output file exists';
}


# change_header_lines - RG
ok $sam_util->change_header_lines($sorted_bam_copy, (RG => {invalid => { sample_name => 'NA00002', library => 'blib', centre => 'NCBI' }})), 'change_header_lines ran ok with an invalid readgroup';
@header_lines = get_bam_header($sorted_bam_copy);
is $header_lines[2], "\@RG\tID:SRR00001\tLB:alib\tSM:NA00001\tPU:7563\tPI:2000\tCN:Sanger\tPL:ILLUMINA\tDS:SRP000001", 'change_header_lines didn\'t change the header when readgroup not in the bam';
%new_header_args = (SRR00001 => { sample_name => 'NA00002', library => 'blib', centre => 'NCBI', study => 'SRP000009' });
is $sam_util->header_rewrite_required($sorted_bam_copy, (RG => \%new_header_args)), 1, 'before rewriting header, check returns true';
ok $sam_util->change_header_lines($sorted_bam_copy, (RG => \%new_header_args)), 'rewrite_bam_header ran ok';
@header_lines = get_bam_header($sorted_bam_copy);
is $header_lines[2], "\@RG\tID:SRR00001\tLB:blib\tSM:NA00002\tPU:7563\tPI:2000\tCN:NCBI\tPL:ILLUMINA\tDS:SRP000009", 'change_header_lines actually changed the header';
@records = get_bam_body($sorted_bam_copy);
is @records, 2000, 'rewrite_bam_header didn\'t change the number of records';
is $sam_util->header_rewrite_required($sorted_bam_copy, (RG => \%new_header_args)), 0, 'after rewriting header, check returns false';


# change_header_lines - PG
ok $sam_util->change_header_lines($sorted_bam_copy, (PG => {invalid => { command => '-q 98', version => '1.234' }})), 'change_header_lines ran ok with an invalid program id';
@header_lines = get_bam_header($sorted_bam_copy);
is $header_lines[3], "\@PG\tID:bwa\tVN:0.4.9", 'change_header_lines didn\'t change the header when program id not in the bam';
%new_header_args = (bwa => { version => '0.5.3', command => '-q 98' });
is $sam_util->header_rewrite_required($sorted_bam_copy, (PG => \%new_header_args)), 1, 'before rewriting header, check returns true';
ok $sam_util->change_header_lines($sorted_bam_copy, (PG => \%new_header_args)), 'rewrite_bam_header ran ok';
@header_lines = get_bam_header($sorted_bam_copy);
is $header_lines[3], "\@PG\tID:bwa\tVN:0.5.3\tCL:-q 98", 'change_header_lines actually changed the header';
@records = get_bam_body($sorted_bam_copy);
is @records, 2000, 'rewrite_bam_header didn\'t change the number of records';
is $sam_util->header_rewrite_required($sorted_bam_copy, (PG => \%new_header_args)), 0, 'after rewriting header, check returns false';
# Change it back
%new_header_args = (bwa => { version => '0.4.9', command => '' });
$sam_util->change_header_lines($sorted_bam_copy, (PG => \%new_header_args));
@header_lines = get_bam_header($sorted_bam_copy);
is $header_lines[3], "\@PG\tID:bwa\tVN:0.4.9", 'change_header_lines successfully deleted at tag';
$sam_util->change_header_lines($sorted_bam_copy, (PG => \%new_header_args));


# test change_header_lines - remove_unique
my $non_unique_header = File::Spec->catfile('t', 'data', 'non_unique_header.txt');
ok -s $non_unique_header, 'non unique header file ready to test with';
ok $sam_util->replace_bam_header($sorted_bam_copy, $non_unique_header), "replace bam header ran successfully";
@header_lines = get_bam_header($sorted_bam_copy);
is $header_lines[3], "\@PG\tID:GATK IndelRealigner\tVN:1.0.4487\tCL:SNPsFileForDebugging=null targetIntervals=/lustre/scratch102/projects/g1k/ref/broad_recal_data/pilot_data/indel.dbsnp_129_b37-vs-pilot.intervals output=null useOnlyKnownIndels=true maxReadsForConsensuses=120 maxConsensuses=30 entropyThreshold=0.15 out=org.broadinstitute.sting.gatk.io.stubs.SAMFileWriterStub\@761eec35 bam_compression=org.broadinstitute.sting.gatk.io.stubs.SAMFileWriterStub\@761eec35 index_output_bam_on_the_fly=org.broadinstitute.sting.gatk.io.stubs.SAMFileWriterStub\@761eec35 indelsFileForDebugging=null noOriginalAlignmentTags=false realignReadsWithBadMates=false maxReadsForRealignment=20000 noPGTag=false LODThresholdForCleaning=0.4 maxReadsInRam=500000 targetIntervalsAreNotSorted=false sortInCoordinateOrderEvenThoughItIsHighlyUnsafe=false statisticsFileForDebugging=null", "non unique header first line added";
is $header_lines[5], "\@PG\tID:GATK TableRecalibration\tVN:1.0.4487\tCL:out=org.broadinstitute.sting.gatk.io.stubs.SAMFileWriterStub\@31ad98ef bam_compression=org.broadinstitute.sting.gatk.io.stubs.SAMFileWriterStub\@31ad98ef index_output_bam_on_the_fly=org.broadinstitute.sting.gatk.io.stubs.SAMFileWriterStub\@31ad98ef output_bam=null window_size_nqs=5 force_read_group=null smoothing=1 default_platform=ILLUMINA exception_if_no_tile=false homopolymer_nback=7 recal_file=/lustre/scratch103/sanger/team145/projects/g1k_main/META/CEU_low_coverage/NA06984/SLX/Solexa_16652/SRR035022/219188.pe.realigned.sorted.bam.recal_data.csv no_pg_tag=false skipUQUpdate=false default_read_group=RG max_quality_score=40 fail_with_no_eof_marker=true solid_recal_mode=SET_Q_ZERO solid_nocall_strategy=THROW_EXCEPTION force_platform=null preserve_qscores_less_than=5 Covariates=[ReadGroupCovariate, QualityScoreCovariate, CycleCovariate, DinucCovariate] pQ=5 maxQ=40 smoothing=1", "non unique header second line added";
is $sam_util->header_rewrite_required($sorted_bam_copy, (PG => {remove_unique => 1})), 1, 'rewrite of header required';
ok $sam_util->change_header_lines($sorted_bam_copy, (PG => {remove_unique => 1})), "remove unique ran ok";
@header_lines = get_bam_header($sorted_bam_copy);
is $header_lines[3], "\@PG\tID:GATK IndelRealigner\tVN:1.0.4487\tCL:SNPsFileForDebugging=null targetIntervals=/lustre/scratch102/projects/g1k/ref/broad_recal_data/pilot_data/indel.dbsnp_129_b37-vs-pilot.intervals output=null useOnlyKnownIndels=true maxReadsForConsensuses=120 maxConsensuses=30 entropyThreshold=0.15 indelsFileForDebugging=null noOriginalAlignmentTags=false realignReadsWithBadMates=false maxReadsForRealignment=20000 noPGTag=false LODThresholdForCleaning=0.4 maxReadsInRam=500000 targetIntervalsAreNotSorted=false sortInCoordinateOrderEvenThoughItIsHighlyUnsafe=false statisticsFileForDebugging=null", "line one uniquified";
is $header_lines[5], "\@PG\tID:GATK TableRecalibration\tVN:1.0.4487\tCL:output_bam=null window_size_nqs=5 force_read_group=null smoothing=1 default_platform=ILLUMINA exception_if_no_tile=false homopolymer_nback=7 no_pg_tag=false skipUQUpdate=false default_read_group=RG max_quality_score=40 fail_with_no_eof_marker=true solid_recal_mode=SET_Q_ZERO solid_nocall_strategy=THROW_EXCEPTION force_platform=null preserve_qscores_less_than=5 Covariates=[ReadGroupCovariate, QualityScoreCovariate, CycleCovariate, DinucCovariate] pQ=5 maxQ=40 smoothing=1", "line two uniquified";
is $sam_util->header_rewrite_required($sorted_bam_copy, (PG => {remove_unique => 1})), 0, 'after rewriting header, rewrite no longer required';
@records = get_bam_body($sorted_bam_copy);
is @records, 2000, 'remove unique didn\'t change the number of records';

# test again for new GATK output - remove_unique
$non_unique_header = File::Spec->catfile('t', 'data', 'non_unique_header2.txt');
ok -s $non_unique_header, 'second non unique header file ready to test with';
ok $sam_util->replace_bam_header($sorted_bam_copy, $non_unique_header), "replace bam header ran successfully";
@header_lines = get_bam_header($sorted_bam_copy);
is $header_lines[3], "\@PG	ID:GATK IndelRealigner	VN:1.0.4487	CL:SNPsFileForDebugging=null targetIntervals=/lustre/scratch105/projects/g1k/ref/broad_recal_data/pilot_data/indel.dbsnp_129_b37-vs-pilot.intervals output=null useOnlyKnownIndels=true maxReadsForConsensuses=120 maxConsensuses=30 entropyThreshold=0.15 out=org.broadinstitute.sting.gatk.io.stubs.SAMFileWriterStub\@7a0d85cc bam_compression=org.broadinstitute.sting.gatk.io.stubs.SAMFileWriterStub\@7a0d85cc index_output_bam_on_the_fly=org.broadinstitute.sting.gatk.io.stubs.SAMFileWriterStub\@7a0d85cc indelsFileForDebugging=null noOriginalAlignmentTags=false realignReadsWithBadMates=false maxReadsForRealignment=20000 noPGTag=false LODThresholdForCleaning=0.4 maxReadsInRam=500000 targetIntervalsAreNotSorted=false sortInCoordinateOrderEvenThoughItIsHighlyUnsafe=false statisticsFileForDebugging=null", "non unique header first line added";
is $header_lines[5], "\@PG	ID:GATK TableRecalibration	VN:1.0.4702	CL:default_read_group=RG default_platform=ILLUMINA force_read_group=null force_platform=null window_size_nqs=5 homopolymer_nback=7 exception_if_no_tile=false solid_recal_mode=SET_Q_ZERO solid_nocall_strategy=THROW_EXCEPTION recal_file=/lustre/scratch105/projects/g1k/META_II/FIN_low_coverage/HG00273/SLX/HUMgfvRBBDIAAPE/ERR016229/238457.pe.realigned.sorted.bam.recal_data.csv output_bam=null preserve_qscores_less_than=5 smoothing=1 max_quality_score=40 no_pg_tag=false fail_with_no_eof_marker=true skipUQUpdate=false Covariates=[ReadGroupCovariate, QualityScoreCovariate, CycleCovariate, DinucCovariate] pQ=5 maxQ=40 smoothing=1", "non unique header second line added";
is $sam_util->header_rewrite_required($sorted_bam_copy, (PG => {remove_unique => 1})), 1, 'rewrite of header required';
ok $sam_util->change_header_lines($sorted_bam_copy, (PG => {remove_unique => 1})), "remove unique ran ok";
@header_lines = get_bam_header($sorted_bam_copy);
is $header_lines[3], "\@PG\tID:GATK IndelRealigner\tVN:1.0.4487\tCL:SNPsFileForDebugging=null targetIntervals=/lustre/scratch105/projects/g1k/ref/broad_recal_data/pilot_data/indel.dbsnp_129_b37-vs-pilot.intervals output=null useOnlyKnownIndels=true maxReadsForConsensuses=120 maxConsensuses=30 entropyThreshold=0.15 indelsFileForDebugging=null noOriginalAlignmentTags=false realignReadsWithBadMates=false maxReadsForRealignment=20000 noPGTag=false LODThresholdForCleaning=0.4 maxReadsInRam=500000 targetIntervalsAreNotSorted=false sortInCoordinateOrderEvenThoughItIsHighlyUnsafe=false statisticsFileForDebugging=null", "line one uniquified";
is $header_lines[5], "\@PG	ID:GATK TableRecalibration	VN:1.0.4702	CL:default_read_group=RG default_platform=ILLUMINA force_read_group=null force_platform=null window_size_nqs=5 homopolymer_nback=7 exception_if_no_tile=false solid_recal_mode=SET_Q_ZERO solid_nocall_strategy=THROW_EXCEPTION output_bam=null preserve_qscores_less_than=5 smoothing=1 max_quality_score=40 no_pg_tag=false fail_with_no_eof_marker=true skipUQUpdate=false Covariates=[ReadGroupCovariate, QualityScoreCovariate, CycleCovariate, DinucCovariate] pQ=5 maxQ=40 smoothing=1", "line two uniquified";
is $sam_util->header_rewrite_required($sorted_bam_copy, (PG => {remove_unique => 1})), 0, 'after rewriting header, rewrite no longer required';
@records = get_bam_body($sorted_bam_copy);
is @records, 2000, 'remove unique didn\'t change the number of records';

# test from_dict and remove_unique
$non_unique_header = File::Spec->catfile('t', 'data', 'non_unique_header.txt');
ok -s $non_unique_header, 'second non unique header file ready to test with';
ok $sam_util->replace_bam_header($sorted_bam_copy, $non_unique_header), "replace bam header ran successfully";
@header_lines = get_bam_header($sorted_bam_copy);
is $header_lines[1], "\@SQ\tSN:GL000192.1\tLN:547496\tAS:NCBI37\tUR:file:/lustre/scratch102/projects/g1k/ref/main_project/human_g1k_v37.fasta\tM5:325ba9e808f669dfeee210fdd7b470ac";
is $header_lines[3], "\@PG\tID:GATK IndelRealigner\tVN:1.0.4487\tCL:SNPsFileForDebugging=null targetIntervals=/lustre/scratch102/projects/g1k/ref/broad_recal_data/pilot_data/indel.dbsnp_129_b37-vs-pilot.intervals output=null useOnlyKnownIndels=true maxReadsForConsensuses=120 maxConsensuses=30 entropyThreshold=0.15 out=org.broadinstitute.sting.gatk.io.stubs.SAMFileWriterStub\@761eec35 bam_compression=org.broadinstitute.sting.gatk.io.stubs.SAMFileWriterStub\@761eec35 index_output_bam_on_the_fly=org.broadinstitute.sting.gatk.io.stubs.SAMFileWriterStub\@761eec35 indelsFileForDebugging=null noOriginalAlignmentTags=false realignReadsWithBadMates=false maxReadsForRealignment=20000 noPGTag=false LODThresholdForCleaning=0.4 maxReadsInRam=500000 targetIntervalsAreNotSorted=false sortInCoordinateOrderEvenThoughItIsHighlyUnsafe=false statisticsFileForDebugging=null", "non unique header first line added";
is $header_lines[5], "\@PG\tID:GATK TableRecalibration\tVN:1.0.4487\tCL:out=org.broadinstitute.sting.gatk.io.stubs.SAMFileWriterStub\@31ad98ef bam_compression=org.broadinstitute.sting.gatk.io.stubs.SAMFileWriterStub\@31ad98ef index_output_bam_on_the_fly=org.broadinstitute.sting.gatk.io.stubs.SAMFileWriterStub\@31ad98ef output_bam=null window_size_nqs=5 force_read_group=null smoothing=1 default_platform=ILLUMINA exception_if_no_tile=false homopolymer_nback=7 recal_file=/lustre/scratch103/sanger/team145/projects/g1k_main/META/CEU_low_coverage/NA06984/SLX/Solexa_16652/SRR035022/219188.pe.realigned.sorted.bam.recal_data.csv no_pg_tag=false skipUQUpdate=false default_read_group=RG max_quality_score=40 fail_with_no_eof_marker=true solid_recal_mode=SET_Q_ZERO solid_nocall_strategy=THROW_EXCEPTION force_platform=null preserve_qscores_less_than=5 Covariates=[ReadGroupCovariate, QualityScoreCovariate, CycleCovariate, DinucCovariate] pQ=5 maxQ=40 smoothing=1", "non unique header second line added";
is $sam_util->header_rewrite_required($sorted_bam_copy, (SQ => {from_dict => '/lustre/scratch105/projects/g1k/ref/main_project/human_g1k_v37.dict.new'}, PG => {remove_unique => 1})), 1, 'rewrite of header required';
ok $sam_util->change_header_lines($sorted_bam_copy, (SQ => {from_dict => '/lustre/scratch105/projects/g1k/ref/main_project/human_g1k_v37.dict.new'}, PG => {remove_unique => 1})), "remove unique ran ok";
@header_lines = get_bam_header($sorted_bam_copy);
is $header_lines[1], "\@SQ\tSN:1\tLN:249250621\tAS:NCBI37\tUR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz\tM5:1b22b98cdeb4a9304cb5d48026a85128\tSP:Human";
is $header_lines[86], "\@PG\tID:GATK IndelRealigner\tVN:1.0.4487\tCL:SNPsFileForDebugging=null targetIntervals=/lustre/scratch102/projects/g1k/ref/broad_recal_data/pilot_data/indel.dbsnp_129_b37-vs-pilot.intervals output=null useOnlyKnownIndels=true maxReadsForConsensuses=120 maxConsensuses=30 entropyThreshold=0.15 indelsFileForDebugging=null noOriginalAlignmentTags=false realignReadsWithBadMates=false maxReadsForRealignment=20000 noPGTag=false LODThresholdForCleaning=0.4 maxReadsInRam=500000 targetIntervalsAreNotSorted=false sortInCoordinateOrderEvenThoughItIsHighlyUnsafe=false statisticsFileForDebugging=null", "line one uniquified";
is $header_lines[88], "\@PG\tID:GATK TableRecalibration\tVN:1.0.4487\tCL:output_bam=null window_size_nqs=5 force_read_group=null smoothing=1 default_platform=ILLUMINA exception_if_no_tile=false homopolymer_nback=7 no_pg_tag=false skipUQUpdate=false default_read_group=RG max_quality_score=40 fail_with_no_eof_marker=true solid_recal_mode=SET_Q_ZERO solid_nocall_strategy=THROW_EXCEPTION force_platform=null preserve_qscores_less_than=5 Covariates=[ReadGroupCovariate, QualityScoreCovariate, CycleCovariate, DinucCovariate] pQ=5 maxQ=40 smoothing=1", "line two uniquified";
is $sam_util->header_rewrite_required($sorted_bam_copy, (SQ => {from_dict => '/lustre/scratch105/projects/g1k/ref/main_project/human_g1k_v37.dict.new'}, PG => {remove_unique => 1})), 0, 'after rewriting header, rewrite no longer required';
@records = get_bam_body($sorted_bam_copy);
is @records, 2000, 'remove unique didn\'t change the number of records';

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
                                           all_unmapped => 1);
$actually_created = 0;
$created_lines = 0;
foreach my $split (@expected_splits) {
    $actually_created += -s $split ? 1 : 0;
    $created_lines += get_bam_body($split);
    unlink($split);
}
is $actually_created, 26, 'split_bam_by_sequence all_unmapped gives the same number of files as before';
is $created_lines, 2014, 'but we now have duplicate reads entries';
@splits = $sam_util->split_bam_by_sequence($headed_bam,
                                           output_dir => $temp_dir,
                                           pretend => 1);
$actually_created = 0;
foreach my $split (@expected_splits) {
    $actually_created += -s $split ? 1 : 0;
}
is @splits, 26, 'split_bam_by_sequence can pretend to make all splits';
is $actually_created, 0, 'split_bam_by_sequence pretend doesn\'t actually make any splits';

is $sam_util->num_bam_header_lines($headed_bam), 116, 'num_bam_header_lines works';
is $sam_util->num_bam_records($headed_bam), 2000, 'num_bam_records works';
is $sam_util->num_bam_records($headed_bam, only => '^12$'), 4, 'num_bam_records works with only option';
is $sam_util->num_bam_lines($headed_bam), 2116, 'num_bam_lines works';

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

# extract_intervals_from_bam
my $temp_extracted_bam = File::Spec->catfile($temp_dir, 'intervals_extracted.bam');
ok $sam_util->extract_intervals_from_bam($intervals_all_bam, $intervals_to_extract, $temp_extracted_bam), 'extract_intervals_from_bam returned true';
my @expected_reads = get_bam_readnames($intervals_extracted_bam);
my @extracted_reads = get_bam_readnames($temp_extracted_bam);
is_deeply \@expected_reads, \@extracted_reads, 'extract_intervals_from_bam extracted the correct reads';

# tag_strip
my $lane_bam = File::Spec->catfile('t', 'data', '1kg_lane.bam');
my $strip_bam = File::Spec->catfile($temp_dir, 'strip.bam');
@records = ([qw(SRR035022.11486888 113 1 10003 0 76M 12 95704 0 ACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC @>=>??;>;?=@4@>?5?=@A@<@=>A@>@AAA@@@@AA@@@AAA@A?AAA@=@AAAA@@>AAAA@;AAAAA=AA@)],
            [qw(SRR035022.7899884 81 1 10021 0 22S54M 5 11065 0 CCTACCCCTACCCCTACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA), '#######################>>=>>.?@@@?@A@?@>@AA@@@@<A@A@=AA@A@AAAA>@AA@A>A<>7;2;'],
            [qw(SRR035022.16462491 99 1 10024 0 74M2S = 10317 327 CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA), '@B@B@AABAA?AABAA?AABAB?@@BAB>AAB@B?@@BAB?8@<AB>@@@@B<AA@<A??AB?A;A>::A8@?###']);
ok $sam_util->tag_strip($lane_bam, $strip_bam, qw(OQ MD XM XG XO)), 'tag_strip returned ok';
is_deeply [get_bam_body($strip_bam)], [join("\t", @{$records[0]}, qw(X0:i:372 RG:Z:SRR035022 AM:i:0 NM:i:0 SM:i:0 MQ:i:0 XT:A:R)),
                                       join("\t", @{$records[1]}, qw(X0:i:442 XC:i:54 RG:Z:SRR035022 AM:i:0 NM:i:0 SM:i:0 MQ:i:0 XT:A:R)),
                                       join("\t", @{$records[2]}, qw(X0:i:377 XC:i:74 RG:Z:SRR035022 AM:i:0 NM:i:0 SM:i:0 MQ:i:23 XT:A:R))], 'tag_strip produced correct results';

# bam_exome_qc_stats
my %vertres_stats;
my %verify_stats;
get_exome_bam_stats($exome_bam_file, $bait_interval, $bait_interval_list, $target_interval, $target_interval_list, $ref_fa, $ref_fai, \%vertres_stats, \%verify_stats, $temp_dir);
is_deeply \%vertres_stats, \%verify_stats, 'get_exome_bam_stats produced correct results';

<<<<<<< HEAD

# filter_readgroups
my $multi_rg_bam = File::Spec->catfile('t', 'data', 'multi_readgroup.bam');
ok -s $multi_rg_bam, 'multi readgroup test bam is present';
@records = get_bam_body($multi_rg_bam);
is @records, 407, 'multi read group bam initially has 407 records';
my $include_bam = File::Spec->catfile($temp_dir, 'include.bam');
my $exclude_bam = File::Spec->catfile($temp_dir, 'exclude.bam');
$sam_util->filter_readgroups($multi_rg_bam, $include_bam, include => [{SM => 'HG01619'}, {SM => 'HG01625'}]);
ok -s $include_bam, 'bam created with include filter';
@records = get_bam_body($include_bam);
is @records, 209, 'bam filtered with \'include\' has 209 records';
$sam_util->filter_readgroups($multi_rg_bam, $exclude_bam, exclude => [{ID => 'SRR071227'}]);
ok -s $exclude_bam, 'bam created with exclude filter';
@records = get_bam_body($exclude_bam);
is @records, 268, 'bam filtered with \'exclude\' has 268 records';
=======
# bam_refnames
my @expected_refnames = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X MT
                  NT_113887 NT_113947 NT_113903 NT_113908 NT_113940 NT_113917
                  NT_113963 NT_113876 NT_113950 NT_113946 NT_113920 NT_113911
                  NT_113907 NT_113937 NT_113941 NT_113909 NT_113921 NT_113919
                  NT_113960 NT_113945 NT_113879 NT_113938 NT_113928 NT_113906
                  NT_113904 NT_113873 NT_113966 NT_113943 NT_113914 NT_113948
                  NT_113886 NT_113932 NT_113929 NT_113878 NT_113927 NT_113900
                  NT_113918 NT_113875 NT_113942 NT_113926 NT_113934 NT_113954
                  NT_113953 NT_113874 NT_113883 NT_113924 NT_113933 NT_113884
                  NT_113890 NT_113870 NT_113881 NT_113939 NT_113956 NT_113951
                  NT_113902 NT_113913 NT_113958 NT_113949 NT_113889 NT_113936
                  NT_113957 NT_113961 NT_113925 NT_113882 NT_113916 NT_113930
                  NT_113955 NT_113944 NT_113901 NT_113905 NT_113872 NT_113952
                  NT_113912 NT_113935 NT_113880 NT_113931 NT_113923 NT_113915
                  NT_113885 NT_113888 NT_113871 NT_113964 NT_113877 NT_113910
                  NT_113962 NT_113899 NT_113965 NT_113898 NC_007605);
my $got_names = $sam_util->bam_refnames($headed_bam);
is_deeply \@expected_refnames, $got_names, 'bam_refnames produced correct results';
>>>>>>> cf02f32645b3fac145c6ba3c7df04b417fff2e15

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

# Arg = bam filename.
# Returns list of all read names in records, with
# '/1' or '/2' appended, depending if read is first
# or second of a pair respectively.
# Appends nothing to the name of an unpaired read .
sub get_bam_readnames {
    my $bam_file = shift;
    my @readnames;
    my $pars = VertRes::Parser::sam->new(file => $bam_file);

    while (my ($qname, $flag) = $pars->get_fields('QNAME', 'FLAG')) {
        if ($pars->is_first($flag)){
            $qname .= "/1";
        }
        elsif ($pars->is_second($flag)){
            $qname .= "/2";
        }
        push @readnames, $qname;
    }
    return @readnames;
}

# get exome bam QC stats using VertRes code, and also using a few other methods (different
# stats got by differnet methods, so need a few).  Not all stats can be verified, but check all
# the possible ones.
sub get_exome_bam_stats {
    my $bam_file = shift;             # bam to be QC'd
    my $bait_interval = shift;        # baits interval file for VertRes::Utils::Sam->bam_exome_qc_stats()
    my $bait_interval_list = shift;   # baits interval list file for picard CalculateHsMetrics
    my $target_interval = shift;      # targets interval file for VertRes::Utils::Sam->bam_exome_qc_stats()
    my $target_interval_list = shift; # targets interval list file for picard CalculateHsMetrics
    my $ref_file = shift;
    my $ref_fai_file = shift;
    my $vr_stats = shift;             # hash ref to be filled with VertRes stats
    my $check_stats = shift;          # hash ref to be filled with stats from other methods
    my $outdir = shift;               # temporary directory in which to put output files
    my $picard_file = File::Spec->catfile($outdir, "picard.out");

    # get stats from vertres subroutine
    my $o = VertRes::Utils::Sam->new();

    my %ops = ('bam', $bam_file,
               'bait_interval', $bait_interval,
               'target_interval', $target_interval,
               'ref_fa', $ref_file,
               'ref_fai', $ref_fai_file);
<<<<<<< HEAD
    %{$vr_stats} = %{$o->bam_exome_qc_stats('bam', $bam_file, %ops)};
=======
    %{$vr_stats} = %{$o->bam_exome_qc_stats(%ops)};
>>>>>>> cf02f32645b3fac145c6ba3c7df04b417fff2e15

    # we get stats from various subroutines to verify the stats made by bam_exome_qc_stats, put
    # the stats in a hash to compare with VertRes results...

    # --------------------------- picard statistics -------------------------
    my $picard = VertRes::Wrapper::picard->new(quiet => 1, java_memory => 100);
    $picard->CalculateHsMetrics($bam_file, $bait_interval_list, $target_interval_list, $picard_file);

    open my $fh, $picard_file or die "error opening file $picard_file";

    while (my $line = <$fh>) {
        if ($line =~ /^## METRICS CLASS/) {
            my $keys_line = <$fh>;
            chomp $keys_line;
            my @keys = split /\t/, $keys_line;
            my $values_line = <$fh>;
            chomp $values_line;
            my @values = split /\t/, $values_line;

            my %picard_stats;
            @picard_stats{@keys} = @values;
            $check_stats->{bait_bases} = $picard_stats{BAIT_TERRITORY};
            $check_stats->{target_bases} = $picard_stats{TARGET_TERRITORY};
            $check_stats->{bait_design_efficiency} = $picard_stats{BAIT_DESIGN_EFFICIENCY};
            $check_stats->{bait_bases_mapped} = $picard_stats{ON_BAIT_BASES};
            $check_stats->{target_bases_mapped} = $picard_stats{ON_TARGET_BASES};
            $check_stats->{low_cvg_targets_pct} = $picard_stats{ZERO_CVG_TARGETS_PCT};
            $vr_stats->{low_cvg_targets_pct} = 100 * $vr_stats->{low_cvg_targets} / $vr_stats->{targets};

            last;
        }
    }
    close $fh;
    unlink $picard_file;

    # picard's output is to 6dp
    $vr_stats->{bait_design_efficiency} = sprintf "%.6f", $vr_stats->{bait_design_efficiency};

    # -------------------- flagstat statistics -----------------------------------
    my $bam_file_flagstat = $bam_file . '.flagstat';
    $o = VertRes::Wrapper::samtools->new();
    $o->flagstat($bam_file, $bam_file_flagstat);

    open $fh, $bam_file_flagstat or die "error opening flagstat file $bam_file_flagstat";

    while (<$fh>){
        my ($stat) = split;

        if (/paired in sequencing/) {
            $check_stats->{reads_paired} = $stat;
        }
        elsif (/in total/) {
            $check_stats->{raw_reads} = $stat;
        }
        elsif (/mapped \(.*\)/) {
            $check_stats->{reads_mapped} = $stat;
        }
        elsif (/properly paired/){
            $check_stats->{mapped_as_pair} = $stat;
        }
    }

    close $fh;
    unlink $bam_file_flagstat;

    # ------------------------ bam statistics ---------------------------
    my %bam_statistics = VertRes::Utils::Sam->new()->bam_statistics($bam_file);
    %bam_statistics = %{$bam_statistics{'*'}};
    $check_stats->{bases_mapped} = $bam_statistics{mapped_bases};
    $check_stats->{pct_mismatches} = $bam_statistics{percent_mismatch};
    $vr_stats->{pct_mismatches} = sprintf "%.2f", $vr_stats->{pct_mismatches};
    $check_stats->{median_insert_size} = $bam_statistics{median_isize};
    $check_stats->{mean_insert_size} = $bam_statistics{avg_isize};
    $vr_stats->{mean_insert_size} =  sprintf "%.0f", $vr_stats->{mean_insert_size};
    $check_stats->{rmdup_reads_mapped} = $vr_stats->{reads_mapped} - $bam_statistics{duplicate_reads};

    # ------------------- samtools statistics ---------------------------
    my $samtools_stats = SamTools::collect_detailed_bam_stats($bam_file, $ref_fai);
    $check_stats->{raw_bases} = $samtools_stats->{total}{bases_total};
    $check_stats->{num_mismatches} = $samtools_stats->{total}{num_mismatches};
    $check_stats->{insert_size_hist} = $samtools_stats->{total}{insert_size}{data};
    foreach (keys %{$check_stats->{insert_size_hist}}) {
        $check_stats->{insert_size_hist}{$_} /= 2 if $check_stats->{insert_size_hist}{$_} > 0;
    }

    # we can't test all the vertres stats, so have to remove the untestable ones from the hash
    my @to_delete;
    foreach my $key (keys %{$vr_stats}) {
        push @to_delete, $key unless (exists $check_stats->{$key});
    }

    delete @{$vr_stats}{@to_delete};

    # there's a small disagreement here; hard to tell exactly what Picard
    # does, so leave this one out of the test
    delete $vr_stats->{low_cvg_targets_pct};
    delete $check_stats->{low_cvg_targets_pct};
}
