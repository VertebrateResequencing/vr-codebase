#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;

BEGIN {
    use Test::Most tests => 21;
    
    use_ok('VertRes::Parser::bam');
}

my $pb = VertRes::Parser::bam->new();
isa_ok $pb, 'VertRes::Parser::ParserI';
isa_ok $pb, 'VertRes::IO';
isa_ok $pb, 'VertRes::Base';

ok my $rh = $pb->result_holder(), 'result_holder returned something';
is ref($rh), 'HASH', 'result_holder returns a hash ref';
is keys %{$rh}, 0, 'the result_holder starts off empty';

ok ! $pb->next_result, 'next_result returns false when we have no file set';

my $b_file = File::Spec->catfile('t', 'data', 'simple2.bam');
ok -e $b_file, 'bam file we will test with exists';
ok $pb->file($b_file), 'file set into parser';

$pb->get_fields('QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'MRNM', 'MPOS', 'ISIZE', 'SEQ', 'QUAL');
ok $pb->next_result, 'next_result now works';
is_deeply $rh, {QNAME => 'IL3_2822:6:1:3:1989', FLAG => '69', RNAME => '*', POS => '0', MAPQ => '0', CIGAR => '*', MRNM => '*', MPOS => '0', ISIZE => '0', SEQ => 'AAAAAAAAAAAAAAAAAAAAAAAAAAANNAAAAAAAANAANAAAAAAAAAAAAAAAAAAAA', QUAL => '44444477774774774447477747($$(447474\'$(($(7777477744477774477'}, 'result_holder contains correct info for first line';
$pb->get_fields('SEQ');
ok $pb->next_result, 'next_result worked again';
is $rh->{SEQ}, 'NANANNAAANNANNANAAAAAAAAANAAANNNNNNNNNNNNNNNNNNANNNNNN', 'result_holder contains correct info for second line';
is $rh->{QNAME}, undef, 'result_holder contains undef for an unrequested field';

# check the last line as well
$pb->get_fields('QNAME', 'NM', 'MD');
while ($pb->next_result) {
    next;
}
is $rh->{QNAME}, 'IL3_2822:6:1:46:1716', 'result_holder contains correct qname for last line';
is $rh->{NM}, 1, 'result_holder contains correct NM for last line';
is $rh->{MD}, '47A6', 'result_holder contains correct MD for last line';

# more tests in a different file
$b_file = File::Spec->catfile('t', 'data', 'simple.bam');
ok -e $b_file, 'second bam file we will test with exists';
$pb = VertRes::Parser::bam->new(file => $b_file);
$pb->get_fields('QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'MRNM', 'MPOS', 'ISIZE', 'SEQ', 'QUAL', 'XT', 'NM', 'SM', 'AM', 'X0', 'X1', 'XM', 'XO', 'XG', 'MD', 'SEQ_LENGTH', 'MAPPED_SEQ_LENGTH');
$pb->next_result;
$rh = $pb->result_holder();
is_deeply $rh, {QNAME => 'IL3_2822:6:1:40:1780', FLAG => 163, RNAME => 'Streptococcus_suis', POS => 4927, MAPQ => 37, CIGAR => '54M', MRNM => 'Streptococcus_suis', MPOS => 5086, ISIZE => 213, SEQ => 'CTAGAAGACGGAAAATCTGCCCGCACAGTTGAGTTCACAGATGAAGAACAAAAA', QUAL => '?>??;:@<?????6)5><?8=;??=;;?2:?>>2<>:>:?:3@97?1764091=', XT => 'U', NM => 0, SM => 37, AM => 37, X0 => 1, X1 => 0, XM => 0, XO => 0, XG => 0, MD => 54, SEQ_LENGTH => 54, MAPPED_SEQ_LENGTH => 54}, 'next_result test on first line';
$pb->next_result;
is $rh->{ISIZE}, -213, 'get_fields second line mate gets a negative isize';
#my $last_fields;
#my @fields;
#while ($pb->next_result) {
#    @fields = $pb->get_fields('QNAME', 'FLAG', 'RG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'MRNM', 'MPOS', 'ISIZE', 'SEQ', 'SEQ_LENGTH', 'QUAL', 'XT');
#    $last_fields = [@fields];
#}
#is_deeply $last_fields, ['IL3_2822:6:1:46:1362', 133, '*', '*', 0, 0, '*', '*', 0, 0, 'CGTTGTAGCTGAGGCTGACATTGAATTTTCACATCTCGAAAGCTTCATCAGTCT', 54, '??,<-(\'4+9/&<A>8((2)4<?835?>.9+(.\'%39@<8<3\'6,.)-==./.(', '*'], 'get_fields test on last line';
#
## hard/soft clipping seq length tests
#$b_file = File::Spec->catfile('t', 'data', 'hard_soft.bam');
#ok -e $b_file, 'hard/soft bam file we will test with exists';
#$pb = VertRes::Parser::bam->new(file => $b_file);
#$pb->next_result;
#is_deeply [$pb->get_fields('SEQ', 'SEQ_LENGTH', 'MAPPED_SEQ_LENGTH')], ['CCCATAGCCCTATCCCTAACCCTAACCCGAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC', 76, 72], 'seq length correct both raw and clipped';

exit;
