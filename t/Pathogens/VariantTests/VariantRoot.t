use strict;
use warnings;
use Data::Dumper;

use Test::More tests => 2;

BEGIN { use_ok( 'Pathogens::Variant::Root' ); }

my $object = Pathogens::Variant::Root->new ();
isa_ok ($object, 'Pathogens::Variant::Root');


