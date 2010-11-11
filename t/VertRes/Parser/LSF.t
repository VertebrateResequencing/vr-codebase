#!/usr/bin/env perl 
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 8;
    use_ok('VertRes::Parser::LSF');
}

my $parser;

$parser = VertRes::Parser::LSF->new(file=>'t/data/lsf-timelimit.o');
is $parser->get('status',0), 'RUNLIMIT';
cmp_ok( $parser->get('time',0), '==', 28765.00);
cmp_ok( $parser->get('memory',0), '==', 1407);
is $parser->get('status'), 'OK';
cmp_ok( $parser->get('time'), '==', 71791.51);
cmp_ok( $parser->get('memory'), '==', 4752);

$parser = VertRes::Parser::LSF->new(file=>'t/data/lsf-memlimit.o');
is $parser->get('status'), 'MEMLIMIT';

exit;
