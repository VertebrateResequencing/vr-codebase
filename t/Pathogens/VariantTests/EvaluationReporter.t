use strict;
use warnings;
use Data::Dumper;

use Test::More tests => 2;

BEGIN { use_ok( 'Pathogens::Variant::EvaluationReporter' ); }

my $object = Pathogens::Variant::EvaluationReporter->new ();
isa_ok ($object, 'Pathogens::Variant::EvaluationReporter');
