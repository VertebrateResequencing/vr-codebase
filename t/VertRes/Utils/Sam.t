#!/usr/bin/perl -w
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 6;
    
    use_ok('VertRes::Utils::Sam');
    use_ok('VertRes::IO');
}

my $sam_util = VertRes::Utils::Sam->new();
isa_ok $sam_util, 'VertRes::Base';

# setup our input files
my $io = VertRes::IO->new();
my $bam1_file = $io->catfile('t', 'data', 'nsort1.bam');
ok -s $bam1_file, 'first bam file ready to test with';
my $bam2_file = $io->catfile('t', 'data', 'nsort2.bam');
ok -s $bam2_file, 'second bam file ready to test with';

# test bams_are_similar
is $sam_util->bams_are_similar($bam1_file, $bam2_file), 1, 'bams are similar';

exit;
