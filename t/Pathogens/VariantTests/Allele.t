use strict;
use warnings;
use Data::Dumper;

use Test::More tests => 2;

BEGIN { use_ok( 'Pathogens::Variant::Allele' ); }

my $object = Pathogens::Variant::Allele->new ();
isa_ok ($object, 'Pathogens::Variant::Allele');
