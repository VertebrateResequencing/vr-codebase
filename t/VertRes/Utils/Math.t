#!/usr/bin/env perl
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 6;
    
    use_ok('VertRes::Utils::Math');
}

my $math_util = VertRes::Utils::Math->new();
isa_ok $math_util, 'VertRes::Base';

is $math_util->histogram_median({1 => 5, 2 => 10, 3 => 5}), 2, 'histogram_median test';
is $math_util->histogram_mean({1 => 3, 2 => 5, 3 => 1, 6 => 2, 10 => 1}), 3.16666666666667, 'histogram_mean_test';
is_deeply {$math_util->histogram_quartiles({1 => 3, 2 => 5, 3 => 1, 6 => 2, 10 => 1})}, {q1 => 1, q2 => 2, q3 => 3}, 'histogram_quartiles_test';
is_deeply {$math_util->histogram_stats({1 => 3, 2 => 5, 3 => 1, 6 => 2, 10 => 1})}, {standard_deviation => 2.64049658629248, mean => 3.16666666666667, total => 12, q1 => 1, q2 => 2, q3 => 3}, 'histogram_stats test';

exit;
