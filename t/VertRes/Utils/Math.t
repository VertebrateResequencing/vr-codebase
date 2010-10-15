#!/usr/bin/env perl
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 3;
    
    use_ok('VertRes::Utils::Math');
}

my $math_util = VertRes::Utils::Math->new();
isa_ok $math_util, 'VertRes::Base';

is $math_util->histogram_median({1 => 5, 2 => 10, 3 => 5}), 2, 'histogram_median test';

exit;
