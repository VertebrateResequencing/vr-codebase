#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;

BEGIN {
    use Test::Most tests => 13;
    
    use_ok('VertRes::Stats::FastQ');
}

my $sfq = VertRes::Stats::FastQ->new();
isa_ok $sfq, 'VertRes::Stats::FastQ';
isa_ok $sfq, 'VertRes::Base';

my $fqc_file = File::Spec->catfile('t', 'data', 'fastq.gz.fastqcheck');
ok -e $fqc_file, 'file we will test with exists';
is $sfq->num_sequences($fqc_file), 7156780, 'num_sequences test';
is $sfq->total_length($fqc_file), 364995780, 'total_length test';
is $sfq->avg_length($fqc_file), '51.00', 'avg_length test';
is $sfq->max_length($fqc_file), 51, 'max_length test';
is_deeply [$sfq->standard_deviations($fqc_file)], ['0.00', 0.02], 'standard_deviations test';

is $sfq->avg_quality($fqc_file), 27.9, 'avg_quality test';

my $qmb_file = File::Spec->catfile('t', 'data', 'qualmapBayesian_simple.txt');
ok -e $qmb_file, 'second file we will test with exists';
my ($before, $after) = $sfq->changed_quality($qmb_file);
is_deeply $before, {'fastq_1' => ['2.50', '2.50'],
                    'fastq_2' => ['2.50', '2.50']}, 'changed_quality gave the correct before result';
is_deeply $after,  {'fastq_1' => ['3.50', '3.50'],
                    'fastq_2' => ['3.50', '3.50']}, 'changed_quality gave the correct after result';

exit;
