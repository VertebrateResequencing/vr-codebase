use strict;
use warnings;
use Data::Dumper;
use File::Basename; 

use Test::More tests => 3;



BEGIN { use_ok( 'Pathogens::Variant::Utils::BamParser' ); }

#get the directory path to the test file
my (undef, $dir) = fileparse($0);

my $samtools_environment_variable_exists = 1 if exists  $ENV{SAMTOOLS};
is($samtools_environment_variable_exists, 1, "SAMTOOLS environment is present found");

my $object = Pathogens::Variant::Utils::BamParser->new();
isa_ok ($object, 'Pathogens::Variant::Utils::BamParser');

$object->fetch_chromosome_size_into_hash("$dir/data/test.bam");
