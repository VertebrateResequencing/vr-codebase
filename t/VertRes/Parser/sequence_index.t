#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;

BEGIN {
    use Test::Most tests => 37;
    
    use_ok('VertRes::Parser::sequence_index');
}

my $sip = VertRes::Parser::sequence_index->new(verbose => 1);
isa_ok $sip, 'VertRes::Parser::ParserI';
isa_ok $sip, 'VertRes::IO';
isa_ok $sip, 'VertRes::Base';

ok my $rh = $sip->result_holder(), 'result_holder returned something';
is ref($rh), 'ARRAY', 'result_holder returns an array';
is @{$rh}, 0, 'the result_holder starts off empty';

ok ! $sip->next_result, 'next_result returns false when we have no file set';

my $si_file = File::Spec->catfile('t', 'data', 'sequence.index');
ok -e $si_file, 'file we will test with exists';
ok $sip->file($si_file), 'file set into parser';

my @expected_data = (['data/NA19238/sequence_read/ERR000018.recal.fastq.gz',
                      'ee05c1a260621d8840ddf3028ebb2355',
                      'ERR000018',
                      'SRP000032',
                      '1000Genomes Project Pilot 2',
                      'BGI',
                      'ERA000013',
                      '',
                      'SRS000212',
                      'NA19238',
                      'YRI_1',
                      'ERX000014',
                      'ILLUMINA',
                      'Illumina Genome Analyzer',
                      'HU1000RADCAASE',
                      'BGI-FC307N0AAXX',
                      'BGI-FC307N0AAXX_5',
                      '',
                      'SINGLE',
                      '',
                      '0',
                      '',
                      '',
                      '9612363',
                      '346045068'],
                     ['data/NA12282/sequence_read/SRR015438_2.recal.fastq.gz',
                      '80d7ee75e062bbd76756df1dc94c6539',
                      'SRR015438',
                      'SRP000033',
                      '1000Genomes Project Pilot 3',
                      'WUGSC',
                      'SRA008537',
                      '',
                      'SRS000619',
                      'NA12282',
                      'CEPH - 2',
                      'SRX004024',
                      'ILLUMINA',
                      'Illumina Genome Analyzer II',
                      '2773138721',
                      '28075',
                      'HWI-EAS289_3150M',
                      '260',
                      'PAIRED',
                      'data/NA12282/sequence_read/SRR015438_1.recal.fastq.gz',
                      '0',
                      '',
                      '',
                      '12938797',
                      '659878647']);

# parse the first line
$sip->next_result;
my $expected = shift @expected_data;
is_deeply $rh, $expected, 'parsed data for first line';

# get info on a particular lane from line 10
is $sip->lane_info('ERR000025', 'sample_name'), 'NA19240', 'lane_info test when not yet reached';
is $rh->[0], $expected->[0], 'using lane_info doesn\'t change our result holder';
$sip->next_result;
is $rh->[0], 'data/NA19238/sequence_read/ERR000019.recal.fastq.gz', 'using lane_info doesn\'t mess with next_result';

# get info on a particular lane from line 5
is $sip->lane_info('ERR000020', 'fastq_file'), 'data/NA19240/sequence_read/ERR000020_2.recal.fastq.gz', 'lane_info test when allready seen';

# parse the last line
while ($sip->next_result) { next; };
$expected = shift @expected_data;
is_deeply $rh, $expected, 'parsed data for last line';

# test get_lanes() on headed
my @all_lanes = $sip->get_lanes;
is $all_lanes[0], 'ERR000018', 'got first lane with get_lanes on headed file';
my %all_lanes = map { $_ => 1 } @all_lanes;
ok ! defined $all_lanes{RUN_ID}, 'RUN_ID did not get treated as a lane';

# try parsing a sequence.index with no header line
$si_file = File::Spec->catfile('t', 'data', 'sequence.index.headerless');
ok -e $si_file, 'headerlessfile we will test with exists';
ok $sip->file($si_file), 'headerlessfile set into parser';

is $sip->lane_info('ERR000044', 'sample_name'), 'NA18550', 'headerless file parse worked';
$sip->next_result;
is $rh->[0], 'data/NA18550/sequence_read/ERR000044_1.recal.fastq.gz', 'got first line correctly';

# test get_lanes() on headerless
@all_lanes = $sip->get_lanes;
is $all_lanes[0], 'ERR000044', 'got first lane with get_lanes on headerless file';
is $all_lanes[-1], 'SRR014220', 'got last lane';
is @all_lanes, 8723, 'got all lanes';
is $rh->[0], 'data/NA18550/sequence_read/ERR000044_1.recal.fastq.gz', 'getting all lanes didn\'t alter our result holder';
while ($sip->next_result) {
    next;
}
is $rh->[0], 'data/NA19093/sequence_read/SRR014220_2.recal.fastq.gz', 'got last line correctly';

my @lanes = $sip->get_lanes(sample_name => 'NA11994');
is $lanes[0], 'SRR003428', 'got first lane with a given sample_name';
is $lanes[-1], 'SRR014158', 'got last lane with a given sample_name';
is @lanes, 120, 'got all lanes with a given sample_name';

@lanes = $sip->get_lanes(ignore_withdrawn => 1);
is @lanes, 8508, 'got all non-withdrawn lanes';
@lanes = $sip->get_lanes(ignore_INSTRUMENT_PLATFORM => 'solid');
is @lanes, 4544, 'got all non-solid lanes';

# problem run_id where it is both withdrawn and not withdrawn
my @answers = $sip->lane_info('ERR000061', 'withdrawn');
is_deeply \@answers, [0, 1], 'knew that a given lane was both withdrawn and not';
is $sip->lane_info('ERR000061', 'withdrawn'), 1, 'knew that the lane was noted as withdrawn more times than not';

# in 2010, sequence.index format changed by adding a new ANALYSIS_GROUP column
$si_file = File::Spec->catfile('t', 'data', 'sequence.index.2010');
ok -e $si_file, '2010 file we will test with exists';
ok $sip->file($si_file), 'file set into parser';
is $sip->lane_info('ERR000018', 'analysis_group'), 'high coverage', 'ANALYSIS_GROUP is parsable';

exit;
