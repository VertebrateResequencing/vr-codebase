#!/usr/bin/env perl
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 7;
    
    use_ok('VertRes::Parser::ParserI');
}

my $pi = VertRes::Parser::ParserI->new();
isa_ok $pi, 'VertRes::IO';
isa_ok $pi, 'VertRes::Base';

ok my $rh = $pi->result_holder(), 'result_holder returned something';
is ref($rh), 'ARRAY', 'result_holder returns an array';

can_ok $pi, qw(next_result);
throws_ok { $pi->next_result } qr/just an interface/, 'next_result throws because we are just an interface';

exit;
