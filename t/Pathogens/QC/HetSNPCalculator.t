#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;

BEGIN {
    use Test::Most ;
    use_ok('Pathogens::QC::HetSNPCalculator');
}

my $reference_size = 2194961;
my $lane_path = '/lane/path/';
my $lane = '1234_4#2';
my $sample_dir = 'qc-sample';
ok my $hsc = Pathogens::QC::HetSNPCalculator->new( reference_size => $reference_size, lane_path => $lane_path, lane => $lane, sample_dir => $sample_dir ), 'create instance of object';

is ( $hsc->reference_size, $reference_size, 'Reference size' );
is ( $hsc->lane_path, $lane_path, 'Lane path' );
is ( $hsc->lane, $lane, 'Lane' );
is ( $hsc->sample_dir, $sample_dir, 'Sample dir' );


print Dumper($hsc);

done_testing();
