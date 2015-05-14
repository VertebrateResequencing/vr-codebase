use strict;
use warnings;
use Data::Dumper;
use File::Basename; 

use Test::More tests => 7;




BEGIN { use_ok( 'Pathogens::Variant::Iterator::Vcf' ); }

my $object;
my $event;

#get the directory path to the test file
my (undef, $dir) = fileparse($0);

#try iterating through a vcf file with a single snp entry
$object = Pathogens::Variant::Iterator::Vcf->new (vcf_file => "$dir/data/single.snp.vcf");
isa_ok ($object, 'Pathogens::Variant::Iterator::Vcf');
$event = $object->next_event;
isa_ok($event, 'Pathogens::Variant::Event');

#try another vcf where a single indel entry is present
$object = Pathogens::Variant::Iterator::Vcf->new (vcf_file => "$dir/data/single.dip.vcf");
$event = $object->next_event;
isa_ok($event, 'Pathogens::Variant::Event');


#try iterating through a vcf thas is gzipped and with single snp
$object = Pathogens::Variant::Iterator::Vcf->new (vcf_file => "$dir/data/single.snp.vcf.gz");
isa_ok ($object, 'Pathogens::Variant::Iterator::Vcf');

my $obj_counter;
while( my $next_event = $object->next_event ) {
    isa_ok($next_event, 'Pathogens::Variant::Event');
    $obj_counter++;
}

is($obj_counter, 1, 'We expect the iterator to return 1 object only');

