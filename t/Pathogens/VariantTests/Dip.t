# -*- perl -*-

# t/001_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'Pathogens::Variant::Event::Dip' ); }

my $object = Pathogens::Variant::Event::Dip->new ();
isa_ok ($object, 'Pathogens::Variant::Event::Dip');
