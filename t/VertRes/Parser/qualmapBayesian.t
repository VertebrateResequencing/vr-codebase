#!/usr/bin/perl -w
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 15;
    
    use_ok('VertRes::Parser::qualmapBayesian');
}

my $qmb = VertRes::Parser::qualmapBayesian->new();
isa_ok $qmb, 'VertRes::Parser::ParserI';
isa_ok $qmb, 'VertRes::IO';
isa_ok $qmb, 'VertRes::Base';

ok my $rh = $qmb->result_holder(), 'result_holder returned something';
is ref($rh), 'ARRAY', 'result_holder returns an array';
is @{$rh}, 0, 'the result_holder starts off empty';

ok ! $qmb->next_result, 'next_result returns false when we have no file set';

my $qmb_file = $qmb->catfile('t', 'data', 'qualmapBayesian.txt');
ok -e $qmb_file, 'file we will test with exists';
ok $qmb->file($qmb_file), 'file set into parser';

ok $qmb->next_result, 'next_result now works';
is_deeply $rh, [qw(1 1 0 20 0 0)], 'result_holder contains correct info for first line';
ok $qmb->next_result, 'next_result worked again';
is_deeply $rh, [qw(1 1 1 2 68 48)], 'result_holder contains correct info for second line';

# check the last line as well
while ($qmb->next_result) {
    next;
}
is_deeply $rh, [qw(2 54 36 29 36 0)], 'result_holder contains correct info for last line';

exit;
