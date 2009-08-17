#!/usr/bin/perl -w
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 15;
    
    use_ok('VertRes::Parser::dict');
}

my $pd = VertRes::Parser::dict->new();
isa_ok $pd, 'VertRes::Parser::ParserI';
isa_ok $pd, 'VertRes::IO';
isa_ok $pd, 'VertRes::Base';

ok my $rh = $pd->result_holder(), 'result_holder returned something';
is ref($rh), 'ARRAY', 'result_holder returns an array ref';
is @{$rh}, 0, 'the result_holder starts off empty';

ok ! $pd->next_result, 'next_result returns false when we have no file set';

my $d_file = $pd->catfile('t', 'data', 'human_b36_male.dict');
ok -e $d_file, 'file we will test with exists';
ok $pd->file($d_file), 'file set into parser';

ok $pd->next_result, 'next_result now works';
is_deeply $rh, [1, 247249719, 'file:/lustre/scratch103/sanger/team145/g1k/ref/human_b36_male.fa', '9ebc6df9496613f373e73396d5b3b6b6'], 'result_holder contains correct info for first line';
ok $pd->next_result, 'next_result worked again';
is $rh->[3], 'b12c7373e3882120332983be99aeb18d', 'result_holder contains correct info for second line';

# check the last line as well
while ($pd->next_result) {
    next;
}
is_deeply $rh, ['NC_007605', 171823, 'file:/lustre/scratch103/sanger/team145/g1k/ref/human_b36_male.fa', '6743bd63b3ff2b5b8985d8933c53290a'], 'result_holder contains correct qname for last line';

exit;
