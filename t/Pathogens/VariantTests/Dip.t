use strict;
use warnings;
use Data::Dumper;

use Test::More tests => 2;

BEGIN { use_ok( 'Pathogens::Variant::Event::Dip' ); }

my $object = Pathogens::Variant::Event::Dip->new (
                                              chromosome => 'test'
                                            , position => '1'
                                            , id => '.'
                                            , reference_allele => 'A'
                                            , alternative_allele => 'CGT'
                                            , quality => 50 
                                            , filter => 'something' 
                                            , info => 'INDEL;DP=179;VDB=0.0490;AF1=1;AC1=2;DP4=0,0,80,67;MQ=60;FQ=-290'
                                            );
isa_ok ($object, 'Pathogens::Variant::Event::Dip');
