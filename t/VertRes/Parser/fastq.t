#!/usr/bin/perl -w
use strict;
use warnings;
use File::Spec;

BEGIN {
    use Test::Most tests => 26;
    
    use_ok('VertRes::Parser::fastq');
}

my $fqp = VertRes::Parser::fastq->new();
isa_ok $fqp, 'VertRes::Parser::ParserI';
isa_ok $fqp, 'VertRes::IO';
isa_ok $fqp, 'VertRes::Base';

ok my $rh = $fqp->result_holder(), 'result_holder returned something';
is ref($rh), 'ARRAY', 'result_holder returns an array';
is @{$rh}, 0, 'the result_holder starts off empty';

ok ! $fqp->next_result, 'next_result returns false when we have no file set';

my $fq_file = File::Spec->catfile('t', 'data', 'SRR001629_1.fastq');
ok -e $fq_file, 'file we will test with exists';
my $gz_file = File::Spec->catfile('t', 'data', 'SRR001629_1.fastq.gz');
ok -e $gz_file, 'gz file we will test with exists';
ok $fqp->file($fq_file), 'file set into parser';

# parse the first sequence
$fqp->next_result;
is_deeply $rh, ['SRR001629.5', 'GGGGGCAATGCTGGGGTCGTCCTCCTCAACTCGCTCCAGGGGCCAGGGGATACCGCTCATATCACTAAGGGCGGTGCCCAGGTAGAGGAGCTCGCGATAGTCCCATTCAATGGACGTGTACCGGATGTTTAGGAGAGGCAGGGAGGCGATGATCTGGCATGTGTGCCGCAGGTGTGTCAGGAGGTCGTCAA', '88888>>BBBB>A@@@@BBBBBBBBBAAAAAAAAAAAAAAAAAAABBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB@@@BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB@@@ABBBBBBBB'], 'parsed data for first sequence';

# get info on the 4th sequence
is $fqp->seq('SRR001629.14'), 'AAAAAAGTAGCCAAATCAACAGATCACATTTAGCATT', 'seq test when not yet reached';
is $rh->[0], 'SRR001629.5', 'using seq doesn\'t change our result holder';
$fqp->next_result;
is $rh->[0], 'SRR001629.9', 'using seq doesn\'t mess with next_result';

# get info on the 3rd sequence
is $fqp->quality('SRR001629.13'), '@@@@@@@@@@@@@>>>>@@>>>>>>>>>><BB>@@>@@@>>@@@@@>AAA>@>>>>>>BBB999B>><<<@@@@B@>>99888889888>>>>>;;;999B<<BBBB<<<>>>', 'quality test when allready seen';

# test sequence_ids
my @ids = $fqp->sequence_ids();
my %ids = map { $_ => 1 } @ids;
is $ids{'SRR001629.5'}, 1, 'sequence_ids gave first sequence id';
is $ids{'SRR001629.2599'}, 1, 'sequence_ids gave last sequence id';
is @ids, 1000, 'sequence_ids gave all ids';

# parse the last line
while ($fqp->next_result) { next; };
is_deeply $rh, ['SRR001629.2599', 'CACGTGTCCTCAACCTTGGCAAAATAAACTTTCAAAATTAACTGAGACCTATCTCAGATTTTCGGGGTTCACAGTAGCAAGGAAGTGGGGTCTGAGAGACGCCCC', 'BBBBBBBBBBBBBBBBBBAAAAAA?BBB?AAA?>>>>@@@@?B??BBBBBBBBBBBBB@@@@AAB@@@@AABBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB'], 'parsed data for last sequence';

ok $fqp->exists('SRR001629.13'), 'exists found an existing sequence';
ok ! $fqp->exists('fake'), 'exists didn\'t find a fake sequence';

# needs to work on .gz files as well
ok $fqp->file($gz_file), 'gz file set into parser';
is $fqp->seq('SRR001629.2598'), 'AATGTCTCCTTGTGAACAGACTTTTGAGTATTTGGCTTTGTTATCCCCCAGAGAATACAAATGTCTCTATGGACACCAAGGTCATAATAACTCCACTTCTCCCATCCCCCTCACACCCTTTGGCAGCCTCATATAT', 'seq test on gz compressed fastq - penultimate read';
is $fqp->seq('SRR001629.1683'), 'CAACAAGTTATTTTAATTGAAAATAAATTTTCCTGACCAACTATTCTGTCAAAACCACATTAAATGAAGATAGCTCAGCAGTGACCAAATCACTATAAAAAGCATTACATGTTATGGGAGAAATGAGTGGGA', 'seq test on gz compressed fastq - read in the middle';
is $fqp->seq('SRR001629.2599'), 'CACGTGTCCTCAACCTTGGCAAAATAAACTTTCAAAATTAACTGAGACCTATCTCAGATTTTCGGGGTTCACAGTAGCAAGGAAGTGGGGTCTGAGAGACGCCCC', 'seq test on gz compressed fastq - last read';


exit;
