#!/usr/bin/perl -w
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 2;
    
    use_ok('VertRes::Wrapper::RecalQual');
}

my $rq = VertRes::Wrapper::RecalQual->new();
isa_ok $rq, 'VertRes::Wrapper::WrapperI';

# just a stub so far...

exit;
