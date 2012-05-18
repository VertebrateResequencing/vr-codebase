# -*- perl -*-

# t/001_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'Pathogens::Variant::Event' ); }

my $object = Pathogens::Variant::Event->new ();
isa_ok ($object, 'Pathogens::Variant::Event');
