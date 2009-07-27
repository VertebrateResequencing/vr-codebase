#!/usr/bin/perl -w
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 29;
    
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

exit;
