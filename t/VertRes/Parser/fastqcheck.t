#!/usr/bin/perl -w
use strict;
use warnings;
use File::Spec;

BEGIN {
    use Test::Most tests => 26;
    
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

my $fastqcheck_file = File::Spec->catfile('t', 'data', 'fastq.gz.fastqcheck');
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

# test the avg quality methods
undef $fqc;
$fqc = VertRes::Parser::fastqcheck->new(file => $fastqcheck_file);
$rh = $fqc->result_holder();
my $expected_quals = [qw(33.6929292929293 32.4489383215369 33.8459214501511 32.920282542886 33.148743718593 33.204843592331 33.832995951417 34.5202429149798 33.6373737373737 33.9384460141271 34.5242914979757 33.6814964610718 34.1830131445905 35.34375 34.0594758064516 31.3387259858443 32.8558467741936 32.8528225806452 29.8375378405651 32.7464646464646 32.8977732793522 32.4495967741936 32.4118831822759 32.8042381432896 33.0604838709677 32.7993951612903 31.6368209255533 32.3256048387097 22.9908814589666 22.5510616784631 24.2992922143579 22.5556680161943 21.9604863221885 25.3313131313131 24.2096774193548 25.0443101711984 24.1553985872856 23.6963562753036 23.314459049545 21.7979797979798 22.3222222222222 19.2127016129032 20.8827098078868 20.7700101317123 19.837044534413 18.9494438827098 18.0485829959514 17.7537993920973 16.8546922300706 17.1151515151515 17.3699596774194)];
my ($bases, $quals) = $fqc->avg_base_quals();
is_deeply $quals, $expected_quals, 'avg_base_quals worked';
$fqc->next_result;
($bases, $quals) = $fqc->avg_base_quals();
is_deeply $quals, $expected_quals, 'avg_base_quals worked again after calling next_result';
$fqc->next_result;
is_deeply $rh, [qw(1 30.3 16.1 21.4 32.0  0.2    0   0   0   0   0   0   0   0
                0   0   0  16   6   6  11  11   0   0   0  11   0   5  17   8
                4  13   0   4   8   7 7  13   5   3  19  11 805   0   0   0   0
                24.5)], 'result_holder contains correct info after calling next_result after calling avg_base_quals';
is $fqc->avg_qual, 27.8919469928644, 'avg_qual test';

exit;
