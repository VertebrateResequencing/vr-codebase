#!/usr/bin/env perl
use strict;
use warnings;
use File::Path qw(rmtree);


BEGIN { unshift(@INC, './modules') }
BEGIN {

    use Test::Most;
    use VertRes::Pipelines::Mapping;
}

ok (VertRes::Pipelines::Mapping->create_graphs('bamcheck -q 20', 't/data/S_suis_P17.fa','t/data/simple2.bam'), 'run bamcheck and plot bamcheck');
ok ((-e 't/data/simple2.bam.bc'), 'bam check file created');
ok ((-d 'simple2.bam_graphs/'), 'graphs directory created');
ok ((-e 'simple2.bam_graphs/quals.png'), 'Quality graph created');

unlink('t/data/simple2.bam.bc');
unlink('t/data/S_suis_P17.fa.gc');
rmtree('simple2.bam_graphs');

done_testing();
