#!/usr/bin/env perl
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 4;
    
    use_ok('VertRes::Utils::Math');
}

my $math_util = VertRes::Utils::Math->new();
isa_ok $math_util, 'VertRes::Base';

is $math_util->histogram_median({1 => 5, 2 => 10, 3 => 5}), 2, 'histogram_median test';
is_deeply {$math_util->histogram_stats({1 => 3, 2 => 5, 3 => 1, 6 => 2, 10 => 1})}, {'standard_deviation' => '2.7579087378049', 'mean' => '3.16666666666667', 'total' => 12, 'median' => '2'}, 'histogram_stats test';

exit;
