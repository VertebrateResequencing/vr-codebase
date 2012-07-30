use strict;
use warnings;
use Data::Dumper;

use Test::More tests => 12;

BEGIN { use_ok( 'Pathogens::Variant::Utils::DP4Parser' ); }

my $object = Pathogens::Variant::Utils::DP4Parser->new ();

#format of the $dp4string:
#ref-forward bases, ref-reverse bases, alt-forward bases, alt-reverse bases
isa_ok ($object, 'Pathogens::Variant::Utils::DP4Parser');

#create a fake dp4 string and parse it:
my $dp4string = '1,2,3,4';
$object->parse($dp4string);

#test the parsing results:
is($object->count_reference_forward_bases, 1, "ref forward base count parsed");
is($object->count_reference_reverse_bases, 2, "ref reverse base count parsed");
is($object->count_alternative_forward_bases, 3, "alt forward base count parsed");
is($object->count_alternative_reverse_bases, 4, "alt reverse base count parsed");

is($object->count_referecence_bases, 1+2, "ref total bases count");
is($object->count_alternative_bases, 3+4, "alt total bases count");

is($object->ratio_forward_reference_bases, 1/(1+3), "ratio forward ref bases");
is($object->ratio_forward_alternative_bases, 3/(3+1), "ratio forward alt bases");

is($object->ratio_reverse_reference_bases, 2/(2+4), "ratio reverse ref bases");
is($object->ratio_reverse_alternative_bases, 4/(4+2), "ratio reverse alt bases");



