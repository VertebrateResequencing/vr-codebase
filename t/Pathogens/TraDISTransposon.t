#!/usr/bin/env perl
use strict;
use warnings;

BEGIN { unshift(@INC, './modules') }
BEGIN {
    use Test::Most tests => 12;
    use_ok('Pathogens::Parser::TraDISTransposon;');
}

####
ok my $tradis_transposon = Pathogens::Parser::TraDISTransposon->new(
  filename   => 't/data/2822_6_1_1000.fastq',
  tag_length => 7
),'load without supplying tag';
is $tradis_transposon->inferred_tag, 'AAAAAAA', 'most common inferred tag with real data';
is $tradis_transposon->percentage_reads_with_tag, 2.4, 'percentage of reads with tag';

####
ok $tradis_transposon = Pathogens::Parser::TraDISTransposon->new(
  filename   => 't/data/tradis_tags.fastq',
  tag_length => 7
),'load without supplying tag with near homogenous data';
is $tradis_transposon->inferred_tag, 'AAAAAAA', 'most common inferred tag with near homogenous data';
is $tradis_transposon->percentage_reads_with_tag, 75, 'percentage of reads with tag';

####
ok $tradis_transposon = Pathogens::Parser::TraDISTransposon->new(
  filename   => 't/data/tradis_tags.fastq',
  tag_length => 7,
  num_reads_to_sample => 3,
),'load supplying number of reads to sample';
is $tradis_transposon->inferred_tag, 'CCCCCCC', 'most common tag in subsample';
is $tradis_transposon->percentage_reads_with_tag, 25, 'percentage of reads with tag';

###
ok $tradis_transposon = Pathogens::Parser::TraDISTransposon->new(
  filename   => 't/data/tradis_tags.fastq',
  tag_length => 7,
  tag        => 'CCCCCCC'
),'user supplied tag';
is $tradis_transposon->percentage_reads_with_tag, 25, 'percentage of reads with supplied tag';
