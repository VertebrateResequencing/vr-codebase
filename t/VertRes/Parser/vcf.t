#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;

BEGIN {
    use Test::Most tests => 25;
    
    use_ok('VertRes::Parser::vcf');
}

my $pvcf = VertRes::Parser::vcf->new();
isa_ok $pvcf, 'VertRes::Parser::ParserI';
isa_ok $pvcf, 'VertRes::IO';
isa_ok $pvcf, 'VertRes::Base';

ok my $rh = $pvcf->result_holder(), 'result_holder returned something';
is ref($rh), 'HASH', 'result_holder returns a hash ref';
is keys %{$rh}, 0, 'the result_holder starts off empty';

ok ! $pvcf->next_result, 'next_result returns false when we have no file set';

my $v_file = File::Spec->catfile('t', 'data', 'vcf.3.3');
ok -e $v_file, 'file we will test with exists';
ok $pvcf->file($v_file), 'file set into parser';

ok $pvcf->next_result, 'next_result now works';
is_deeply $rh, {CHROM => 1,
                POS => 10002,
                ID => undef,
                REF => 'A',
                ALT => ['C'],
                QUAL => 100,
                FILTER => 0,
                INFO => { NS => 75,
                          DP => 988,
                          PHASED => 1,
                          AN => 166,
                          AC => 12 },
                FORMAT => [qw(GT GQ DP HQ)],
                SAMPLES => { NA06985 => { GT => '0|0', GQ => 123, DP => 0, HQ => '123,123' },
                             NA06986 => { GT => '0|1', GQ => 94, DP => 14, HQ => '123,94' },
                             NA06994 => { GT => '1|0', GQ => 94, DP => 8, HQ => '94,121' } } }, 'result_holder contains correct info for first line';
ok $pvcf->next_result, 'next_result worked again';
is $rh->{POS}, 10003, 'result_holder contains correct info for second line';

# check the last line as well
while ($pvcf->next_result) {
    next;
}
is $rh->{POS}, 10234, 'result_holder contains correct POS for last line';
is $rh->{ID}, 'rs001', 'result_holder contains correct ID for last line';
is $rh->{INFO}->{DB}, 1, 'result_holder had DB true for last line';
is $rh->{SAMPLES}->{NA06994}->{DP}, 14, 'result_holder contains correct sample detail for last line';

# test match_filter
$pvcf->_seek_first_result;
is $pvcf->match_filter(0), 0, 'match_filter(0) could be set';
my $count = 0;
while ($pvcf->next_result) {
    $count++;
    next;
}
is $count, 3, '3 results with filter 0';
$pvcf->_seek_first_result;
is $pvcf->match_filter('q10'), 'q10', 'match_filter(q10) could be set';
$count = 0;
my $filter;
while ($pvcf->next_result) {
    $filter = $rh->{FILTER};
    $count++;
    next;
}
is $count, 1, '1 result with filter q10';
is $filter, 'q10', 'that result had a filter of q10';
$pvcf->_seek_first_result;
is $pvcf->match_filter(undef), undef, 'match_filter(undef) could be (un)set';
$count = 0;
while ($pvcf->next_result) {
    $count++;
    next;
}
is $count, 4, '4 results with no filter';

exit;
