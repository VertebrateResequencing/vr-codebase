#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;

BEGIN {
    use Test::Most tests => 15;
    
    use_ok('VertRes::Parser::pileup');
}

my $pp = VertRes::Parser::pileup->new();
isa_ok $pp, 'VertRes::Parser::ParserI';
isa_ok $pp, 'VertRes::IO';
isa_ok $pp, 'VertRes::Base';

ok my $rh = $pp->result_holder(), 'result_holder returned something';
is ref($rh), 'ARRAY', 'result_holder returns an array';
is @{$rh}, 0, 'the result_holder starts off empty';

ok ! $pp->next_result, 'next_result returns false when we have no file set';

my $p_file = File::Spec->catfile('t', 'data', 'pileup.txt');
ok -e $p_file, 'file we will test with exists';
ok $pp->file($p_file), 'file set into parser';

ok $pp->next_result, 'next_result now works';
is_deeply $rh, [qw(Streptococcus_suis 4927 c 1 ^h. ?)], 'result_holder contains correct info for first line';
ok $pp->next_result, 'next_result worked again';
is_deeply $rh, [qw(Streptococcus_suis 4928 t 1 . >)], 'result_holder contains correct info for second line';

# check the last line as well
while ($pp->next_result) {
    next;
}
is_deeply $rh, [qw(Streptococcus_suis 4936 g 1 . ?)], 'result_holder contains correct info for last line';

exit;
