# -*- perl -*-

# t/001_load.t - check module loading and create testing directory

use Test::More tests => 2;

BEGIN { use_ok( 'Pathogens::Variant::Sample' ); }

my $object = Pathogens::Variant::Sample->new (name => 'a dummy sample');
isa_ok ($object, 'Pathogens::Variant::Sample');
