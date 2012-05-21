use strict;
use warnings;
use Data::Dumper;

use Test::More tests => 2;

BEGIN { use_ok( 'Pathogens::Variant::Evaluator' ); }

my $object = Pathogens::Variant::Evaluator->new ();
isa_ok ($object, 'Pathogens::Variant::Evaluator');
