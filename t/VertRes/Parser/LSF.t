#!/usr/bin/env perl 
use strict;
use warnings;
use File::Spec;

BEGIN {
    use Test::Most tests => 18;
    use_ok('VertRes::Parser::LSF');
}

my $parser;

# original-style usage
my $lsf_file = File::Spec->catfile(qw(t data lsf-timelimit.o));
$parser = VertRes::Parser::LSF->new(file => $lsf_file);
is $parser->get('status', 0), 'RUNLIMIT', 'status of first record is RUNLIMIT';
is $parser->get('time', 0), 28801, 'time test';
is $parser->get('memory', 0), 1407, 'memory test';
is $parser->get('status'), 'OK', 'status of last record is OK';
is $parser->get('time'), 71870, 'time of last record is correct';
is $parser->get('memory'), 4752, 'memory of last record is correct';
is $parser->nrecords, 2, 'nrecords test';

my $lsf_file2 = File::Spec->catfile(qw(t data lsf-memlimit.o));
$parser = VertRes::Parser::LSF->new(file => $lsf_file2);
is $parser->get('status'), 'MEMLIMIT', 'status of a MEMLIMIT failed job is correct';

# new-style usage and new methods
$parser = VertRes::Parser::LSF->new(file => $lsf_file);
my $count = 0;
while ($parser->next_result) {
    $count++;
}
is $count, 2, 'next_result gives the correct count';
is $parser->nrecords, 2, 'nrecords still works after using next_result';

my $rh = $parser->result_holder;
is_deeply $rh, ['perl -w _129P2_2:89997751-99997751.pl', 'OK', 4752, 71870, 71791.51, '1.00', 'basement'], 'result_holder contains correct data for last result';

is $parser->status, 'OK', 'status() test';
is $parser->time, 71870, 'time() test';
is $parser->cpu_time, 71791.51, 'cpu_time() test';
is $parser->idle_factor, '1.00', 'idle_factor() test';
is $parser->memory, 4752, 'memory() test';
is $parser->cmd, 'perl -w _129P2_2:89997751-99997751.pl', 'cmd() test';

exit;
