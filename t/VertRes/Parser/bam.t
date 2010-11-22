#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;

BEGIN {
    use Test::Most tests => 45;
    
    use_ok('VertRes::Parser::bam');
    use_ok('VertRes::Utils::FileSystem');
    use_ok('VertRes::Wrapper::samtools');
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
$pb->get_fields('QNAME', 'FLAG', 'RG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'MRNM', 'MPOS', 'ISIZE', 'SEQ', 'SEQ_LENGTH', 'QUAL', 'XT');
while ($pb->next_result) {
    next;
}
is_deeply $rh, { QNAME => 'IL3_2822:6:1:46:1362', FLAG => 133, RG => '*', RNAME => '*', POS => 0, MAPQ => 0, CIGAR => '*', MRNM => '*', MPOS => 0, ISIZE => 0, SEQ => 'CGTTGTAGCTGAGGCTGACATTGAATTTTCACATCTCGAAAGCTTCATCAGTCT', SEQ_LENGTH => 54, QUAL => '??,<-(\'4+9/&<A>8((2)4<?835?>.9+(.\'%39@<8<3\'6,.)-==./.(', XT => '*'}, 'next_result test on last line';

# test region();
$pb = VertRes::Parser::bam->new(file => $b_file);
my $c = 0;
foreach my $region ('Streptococcus_suis:29888-50000', 'Streptococcus_suis:80000-90000', 'Streptococcus_suis:60000-70000') {
    $pb->region($region);
    while ($pb->next_result) {
        $c++;
    }
}
is $c, 18, 'getting multiple regions test';

# test writing
my $fsu = VertRes::Utils::FileSystem->new();
my $temp_dir = $fsu->tempdir();
my $out_bam1 = File::Spec->catfile($temp_dir, 'out1.bam');
my $out_bam2 = File::Spec->catfile($temp_dir, 'out2.bam');
my $out_bam3 = File::Spec->catfile($temp_dir, 'out3.bam');
$pb->region('');
while ($pb->next_result) {
    $pb->write_result($out_bam1);
}
$pb->close;
is_deeply [get_bam_header($out_bam1)], [get_bam_header($b_file)], 'having set input region to empty string, header of an output file matches the input file';
is scalar(get_bam_body($out_bam1)), 2000, 'output bam has correct number of lines';

$b_file = File::Spec->catfile('t', 'data', '1kg_lane.bam');
ok -e $b_file, '1kg bam file we will test with exists';
$pb = VertRes::Parser::bam->new(file => $b_file);
$c = 0;
$pb->ignore_tags_on_write(qw(OQ XM XG XO));
while ($pb->next_result) {
    $c++;
    if ($c == 1) {
        $pb->write_result($out_bam2);
    }
    elsif ($c == 2) {
        $pb->write_result($out_bam3);
    }
    else {
        last;
    }
}
undef $pb;
my @g1k_lines = grep { s/\tOQ:\S+|\tXM:\S+|\tXG:\S+|\tXO:\S+//g } get_bam_body($b_file);
is_deeply [get_bam_body($out_bam2)], [$g1k_lines[0]], 'ignore_tags_on_write has its effect';
is_deeply [get_bam_body($out_bam3)], [$g1k_lines[1]], 'we can output multiple files in a single next_result loop';

# filtering methods
$pb = VertRes::Parser::bam->new(file => $b_file);
$pb->required_readgroup('SRR035022');
$c = 0;
while ($pb->next_result) {
    $c++;
}
is $c, 3, 'required_readgroup with the correct RG gives all results';
$pb = VertRes::Parser::bam->new(file => $b_file);
$pb->required_readgroup('foo');
$c = 0;
while ($pb->next_result) {
    $c++;
}
is $c, 0, 'required_readgroup with an unused RG gives no results';
$pb = VertRes::Parser::bam->new(file => $b_file);
$c = 0;
while ($pb->next_result) {
    $c++;
}
is $c, 3, 'no required_readgroup in a new instance gives us all results';

$pb = VertRes::Parser::bam->new(file => $b_file);
$pb->required_library('foo');
$c = 0;
while ($pb->next_result) {
    $c++;
}
is $c, 0, 'required_library with an unused lib gives no results';
$pb = VertRes::Parser::bam->new(file => $b_file);
$pb->required_library('Solexa-16652');
$c = 0;
while ($pb->next_result) {
    $c++;
}
is $c, 3, 'required_library with the correct lib gives all results';

$pb = VertRes::Parser::bam->new(file => $b_file);
$pb->required_flag(32);
$c = 0;
while ($pb->next_result) {
    $c++;
}
is $c, 2, 'required_flag gives a subset of results';
$pb = VertRes::Parser::bam->new(file => $b_file);
$pb->filtering_flag(32);
$c = 0;
while ($pb->next_result) {
    $c++;
}
is $c, 1, 'filtering_flag gives the opposite subset of results';

$pb = VertRes::Parser::bam->new(file => $b_file);
$pb->flag_selector(mate_reverse => 1, '1st_in_pair' => 1);
$c = 0;
while ($pb->next_result) {
    $c++;
}
is $c, 2, 'flag_selector worked with two true flag names';
$pb = VertRes::Parser::bam->new(file => $b_file);
$pb->flag_selector(mate_reverse => 0, '1st_in_pair' => 0);
$c = 0;
while ($pb->next_result) {
    $c++;
}
is $c, 0, 'flag_selector worked with two false flag names';
$pb = VertRes::Parser::bam->new(file => $b_file);
$pb->flag_selector(mate_reverse => 0, '1st_in_pair' => 1);
$c = 0;
while ($pb->next_result) {
    $c++;
}
is $c, 1, 'flag_selector worked with a mix of true and false flag names';

$pb = VertRes::Parser::bam->new(file => $b_file);
$pb->minimum_quality(30);
$c = 0;
while ($pb->next_result) {
    $c++;
}
is $c, 0, 'minimum_quality gives a subset of results';
$pb = VertRes::Parser::bam->new(file => $b_file);
$pb->minimum_quality(0);
$c = 0;
while ($pb->next_result) {
    $c++;
}
is $c, 3, 'minimum_quality(0) gives all results';


# hard/soft clipping seq length tests
$b_file = File::Spec->catfile('t', 'data', 'hard_soft.bam');
ok -e $b_file, 'hard/soft bam file we will test with exists';
ok $pb->file($b_file), 'we can set a new file without having to create a new instance';
$pb->get_fields('SEQ', 'SEQ_LENGTH', 'MAPPED_SEQ_LENGTH');
$rh = $pb->result_holder;
$pb->next_result;
is_deeply $rh, { SEQ => 'CCCATAGCCCTATCCCTAACCCTAACCCGAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC', SEQ_LENGTH => 76, MAPPED_SEQ_LENGTH => 72}, 'seq length correct both raw and clipped';

exit;

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
