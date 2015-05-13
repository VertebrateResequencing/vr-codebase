#!/usr/bin/env perl
use strict;
use warnings;

BEGIN {
    use Test::Most ;
    use_ok('Pathogens::QC::HetSNPCalculator');
}

my $reference_size = 2194961;
ok my $hsc = Pathogens::QC::HetSNPCalculator->new(reference_size => $reference_size), 'create instance of object';


done_testing();
