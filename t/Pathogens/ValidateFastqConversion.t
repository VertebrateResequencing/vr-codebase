#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Test::MockObject;

BEGIN { unshift(@INC, './modules') }
BEGIN {
    use Test::Most ;
    use_ok('Pathogens::Import::ValidateFastqConversion');

    my $irods_wrapper = Test::MockObject->new();
    $irods_wrapper->fake_module( 'VertRes::Wrapper::iRODS', test => sub{1} );
    $irods_wrapper->fake_new( 'VertRes::Wrapper::iRODS' );
    $irods_wrapper->mock('find_file_by_name', sub{ 1 });
    $irods_wrapper->mock('get_total_reads', sub{ 7157780 });
}

# Create a valid object.
ok my $validate = Pathogens::Import::ValidateFastqConversion->new(
     fastqcheck_filenames => ['t/data/2822_6_1_1000.fastq.fastqcheck','t/data/fastq.gz.fastqcheck'],
     irods_filename       => '1234_5#6.bam'
    ),'Initialise valid object';

# Check values returned by sum_fastq
is $validate->_sum_fastq_reads(),7157780,'Sum of two fastqcheck files.';

# Check is_total_reads_valid
is $validate->is_total_reads_valid, 1,'Check is_total_reads_valid()';


# Check error case 
ok my $bad_validate = Pathogens::Import::ValidateFastqConversion->new(
     fastqcheck_filenames => ['t/data/2822_6_1_1000.fastq.fastqcheck'],
     irods_filename       => '1234_5#6.bam'
    ),'Initialise valid object with one fastq';
is $bad_validate->is_total_reads_valid, 0,'Check is_total_reads_valid() returns 0 on error.';


done_testing();
