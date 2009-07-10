#!/usr/bin/perl -w
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 22;
    
    use_ok('VertRes::Parser::fastqcheck');
}

my $fqc = VertRes::Parser::fastqcheck->new();
isa_ok $fqc, 'VertRes::Parser::ParserI';
isa_ok $fqc, 'VertRes::IO';
isa_ok $fqc, 'VertRes::Base';

ok my $rh = $fqc->result_holder(), 'result_holder returned something';
is ref($rh), 'ARRAY', 'result_holder returns an array';
is @{$rh}, 0, 'the result_holder starts off empty';

ok ! $fqc->next_result, 'next_result returns false when we have no file set';

my $fastqcheck_file = $fqc->catfile('t', 'data', 'fastq.gz.fastqcheck');
ok -e $fastqcheck_file, 'file we will test with exists';
ok $fqc->file($fastqcheck_file), 'file set into parser';
is $fqc->num_sequences, 7156780, 'num_sequences test';
is $fqc->total_length, 364995780, 'total_length test';
is $fqc->avg_length, '51.00', 'avg_length test';
is $fqc->max_length, 51, 'max_length test';
is_deeply [$fqc->standard_deviations], ['0.00', 0.02], 'standard_deviations test';

ok $fqc->next_result, 'next_result now works';
is_deeply $rh, [qw(0 25.4 22.0 20.4 28.7  3.6    0   0   0   0   0  36  19  26
                3 1  16  16   5   7  19  29  10  13  17  17  21  18  17  20  16
                23  22  20  18  23 20  23  27  27  31  41  91 114 115  40  20
                15.2)], 'result_holder contains correct info for total line';
ok $fqc->next_result, 'next_result worked again';
is_deeply $rh, [qw(1 30.3 16.1 21.4 32.0  0.2    0   0   0   0   0   0   0   0
                0   0   0  16   6   6  11  11   0   0   0  11   0   5  17   8
                4  13   0   4   8   7 7  13   5   3  19  11 805   0   0   0   0
                24.5)], 'result_holder contains correct info for first position';

# and check we can get the first row without asking for any header info first
undef $fqc;
$fqc = VertRes::Parser::fastqcheck->new(file => $fastqcheck_file);
$rh = $fqc->result_holder();
ok $fqc->next_result, 'next_result worked again';
is_deeply $rh, [qw(0 25.4 22.0 20.4 28.7  3.6    0   0   0   0   0  36  19  26
                3 1  16  16   5   7  19  29  10  13  17  17  21  18  17  20  16
                23  22  20  18  23 20  23  27  27  31  41  91 114 115  40  20
                15.2)], 'result_holder contains correct info for total line';

# and check the last result
while ($fqc->next_result) {
    next;
}
is_deeply $rh, [qw(51  21.0 19.0 19.9 36.4  3.7    0   0   0   0   0 172   0  94
                0   0   0 102   0   0   0  97   0  46   0 108  40  37   0   0
                0  50  24  20  16  47  36  17   0  52  13  21   0   0   0   0
                0 10.5)], 'result_holder contains correct info for last line';

exit;
