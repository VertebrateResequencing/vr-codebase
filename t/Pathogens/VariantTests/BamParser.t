use strict;
use warnings;

use Test::More;
use Data::Dumper;

BEGIN { use_ok( 'Pathogens::Variant::Utils::BamParser' ); }

my $object = Pathogens::Variant::Utils::BamParser->new(bam => 'dummy');
isa_ok ($object, 'Pathogens::Variant::Utils::BamParser');


done_testing;
