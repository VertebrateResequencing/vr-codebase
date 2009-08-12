#!/usr/bin/perl -w
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 3;
    
    use_ok('VertRes::Utils::Seq');
}

my $seq_util = VertRes::Utils::Seq->new();
isa_ok $seq_util, 'VertRes::Base';

is $seq_util->rev_com('aTnCG'), 'CGNAT', 'rev_com test';

exit;
