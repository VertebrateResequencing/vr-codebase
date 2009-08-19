#!/usr/bin/perl -w
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 48;
    
    use_ok('VertRes::Parser::sam');
}

my $ps = VertRes::Parser::sam->new();
isa_ok $ps, 'VertRes::Parser::ParserI';
isa_ok $ps, 'VertRes::IO';
isa_ok $ps, 'VertRes::Base';

ok my $rh = $ps->result_holder(), 'result_holder returned something';
is ref($rh), 'HASH', 'result_holder returns a hash ref';
is keys %{$rh}, 0, 'the result_holder starts off empty';

ok ! $ps->next_result, 'next_result returns false when we have no file set';

my $s_file = $ps->catfile('t', 'data', 'simple.sam');
ok -e $s_file, 'file we will test with exists';
ok $ps->file($s_file), 'file set into parser';

ok $ps->next_result, 'next_result now works';
my %col_to_name = (0  => 'QNAME',
                   1  => 'FLAG',
                   2  => 'RNAME',
                   3  => 'POS',
                   4  => 'MAPQ',
                   5  => 'CIGAR',
                   6  => 'MRNM',
                   7  => 'MPOS',
                   8  => 'ISIZE',
                   9  => 'SEQ',
                   10 => 'QUAL');
is_deeply $rh, {QNAME => 'IL3_2822:6:1:3:1989', FLAG => '69', RNAME => '*', POS => '0', MAPQ => '0', CIGAR => '*', MRNM => '*', MPOS => '0', ISIZE => '0', SEQ => 'AAAAAAAAAAAAAAAAAAAAAAAAAAANNAAAAAAAANAANAAAAAAAAAAAAAAAAAAAA', QUAL => '44444477774774774447477747($$(447474\'$(($(7777477744477774477'}, 'result_holder contains correct info for first line';
ok $ps->next_result, 'next_result worked again';
is $rh->{SEQ}, 'NANANNAAANNANNANAAAAAAAAANAAANNNNNNNNNNNNNNNNNNANNNNNN', 'result_holder contains correct info for second line';

# check the last line as well
while ($ps->next_result) {
    next;
}
is $rh->{QNAME}, 'IL3_2822:6:1:46:1716', 'result_holder contains correct qname for last line';
is $rh->{NM}, 1, 'result_holder contains correct NM for last line';
is $rh->{MD}, '47A6', 'result_holder contains correct MD for last line';

# test all the flag-related methods
my $flag = $rh->{FLAG};
is $flag, 163, 'last flag was correct';
is $ps->is_sequencing_paired($flag), 1, 'is_sequencing_paired test';
is $ps->is_mapped_paired($flag), 1, 'is_mapped_paired test';
is $ps->is_mapped($flag), 1, 'is_mapped test';
is $ps->is_mate_mapped($flag), 1, 'is_mate_mapped test';
is $ps->is_reverse_strand($flag), 0, 'is_reverse_strand test';
is $ps->is_mate_reverse_strand($flag), 1, 'is_mate_reverse_strand test';
is $ps->is_first($flag), 0, 'is_first test';
is $ps->is_second($flag), 1, 'is_second test';
is $ps->is_primary($flag), 1, 'is_primary test';
is $ps->passes_qc($flag), 1, 'passes_qc test';
is $ps->is_duplicate($flag), 0, 'is_duplicate test';

# test the header methods on this headerless sam, then on a headed sam
is $ps->sam_version, undef, 'no header info on headless sam';
my $headed_sam = $ps->catfile('t', 'data', 'NA11918SLX.headed.sam');
ok -e $headed_sam, 'headed sam file ready for testing';
ok $ps->file($headed_sam), 'headed sam set into parser';
ok $ps->next_result, 'next_result worked on a headed sam';
is $rh->{QNAME}, 'SRR003436.685', 'QNAME fine for first record';
is $ps->sam_version, '1.0', 'sam_version worked on headed sam';
is $ps->group_order, 'none', 'group_order correct';
is $ps->sort_order, 'coordinate', 'sort_order correct';
is $ps->program, 'bwa', 'program correct';
is $ps->program_version, '0.4.9', 'program_version correct';
is $ps->command_line, 'bwa -foo bar', 'command_line correct';
is $ps->sequence_info(1, 'LN'), '247249719', 'sequence_info 1st line specific correct';
is $ps->sequence_info('NT_113947', 'LN'), '4262', 'sequence_info last line specific correct';
is $ps->sequence_info(16, 'UR'), 'file:/nfs/sf8/G1K/ref/human_b36_female.fa', 'sequence_info other specific correct';
my %all_rgs = $ps->readgroup_info();
is_deeply [sort keys %all_rgs], ['SRR003435', 'SRR003436', 'SRR003447'], 'readgroup_info all had all the readgroups';
my $rg_info = $all_rgs{SRR003435};
my $expected_info = {PL => 'ILLUMINA',
                     PU => 'BI.PE.080728_SL-XBE_0001_FC3044GAAXX.080731_SL-XBE_0007_FC3044GAAXX.1',
                     LB => 'Solexa-6391',
                     PI => 500,
                     SM => 'NA11918',
                     CN => 'BI'};
is_deeply $rg_info, $expected_info, 'readgroup info from all hash contains the correct info';
is_deeply {$ps->readgroup_info('SRR003435')}, $expected_info, 'readgroup_info(id) gave correct info';
is $ps->readgroup_info('SRR003435', 'LB'), 'Solexa-6391', 'readgroup_info specific was correct';
while ($ps->next_result) {
    next;
}
is $rh->{QNAME}, 'SRR003447.1000', 'QNAME fine for last record';

exit;
