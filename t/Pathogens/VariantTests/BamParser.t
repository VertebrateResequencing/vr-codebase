use strict;
use warnings;

use Test::More;
use Data::Dumper;

BEGIN { use_ok( 'Pathogens::Variant::Utils::BamParser' ); }

my $object = Pathogens::Variant::Utils::BamParser->new();
isa_ok ($object, 'Pathogens::Variant::Utils::BamParser');

$object->fetch_chromosome_size_into_hash('./data/test.bam');


done_testing;
