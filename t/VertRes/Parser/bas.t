#!/usr/bin/perl -w
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 13;
    
    use_ok('VertRes::Parser::bas');
}

my $basp = VertRes::Parser::bas->new();
isa_ok $basp, 'VertRes::Parser::ParserI';
isa_ok $basp, 'VertRes::IO';
isa_ok $basp, 'VertRes::Base';

ok my $rh = $basp->result_holder(), 'result_holder returned something';
is ref($rh), 'ARRAY', 'result_holder returns an array';
is @{$rh}, 0, 'the result_holder starts off empty';

ok ! $basp->next_result, 'next_result returns false when we have no file set';

my $bas_file = $basp->catfile('t', 'data', 'example2.bas');
ok -e $bas_file, 'file we will test with exists';
ok $basp->file($bas_file), 'file set into parser';

ok $basp->next_result, 'next_result now works';

is_deeply $rh, [qw(NA00001.SLX.bwa.SRP000001.2009.08 SRP000001 NA00001 SLX alib SRR00001 115000 62288 2000 1084 1084 1070 2.05 23.32 286 74.10 275 48)], 'first result was correct';

while ($basp->next_result) {
    next;
}

is_deeply $rh, [qw(NA00003.SLX.bwa.SRP000003.2009.08 SRP000003 NA00003 SLX alib3 SRR00003 115000 62288 2000 1084 1084 1070 2.05 23.32 286 74.10 275 46)], 'last result was correct';

exit;
