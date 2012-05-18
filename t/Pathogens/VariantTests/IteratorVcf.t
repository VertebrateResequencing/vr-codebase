# -*- perl -*-

# t/001_load.t - check module loading and create testing directory

use Test::More tests => 6;
use Data::Dumper;

BEGIN { use_ok( 'Pathogens::Variant::Iterator::Vcf' ); }

my $object;
my $event;

#try iterating through a vcf file with a single snp entry
$object = Pathogens::Variant::Iterator::Vcf->new (vcf_file => '/Users/fy2/sandbox/refactor/VariantCallFramework/t/data/single.snp.vcf');
isa_ok ($object, 'Pathogens::Variant::Iterator::Vcf');
$event = $object->next_event;
isa_ok($event, 'Pathogens::Variant::Event::Snp');

#try another vcf where a single indel entry is present
$object = Pathogens::Variant::Iterator::Vcf->new (vcf_file => '/Users/fy2/sandbox/refactor/VariantCallFramework/t/data/single.dip.vcf');
$event = $object->next_event;
isa_ok($event, 'Pathogens::Variant::Event::Dip');


#try iterating through a vcf thas is gzipped and with single snp
$object = Pathogens::Variant::Iterator::Vcf->new (vcf_file => '/Users/fy2/sandbox/refactor/VariantCallFramework/t/data/single.snp.vcf.gz');
isa_ok ($object, 'Pathogens::Variant::Iterator::Vcf');
$event = $object->next_event;
isa_ok($event, 'Pathogens::Variant::Event::Snp');

