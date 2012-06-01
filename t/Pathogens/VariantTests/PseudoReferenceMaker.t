use strict;
use warnings;

use Test::More;
use Data::Dumper;

BEGIN { use_ok( 'Pathogens::Variant::Utils::PseudoReferenceMaker' ); }

my $object = Pathogens::Variant::Utils::PseudoReferenceMaker->new;
isa_ok ($object, 'Pathogens::Variant::Utils::PseudoReferenceMaker');


done_testing;
