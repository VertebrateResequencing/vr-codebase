use strict;
use warnings;

use Test::More;
use Data::Dumper;

BEGIN { use_ok( 'Pathogens::Variant::Utils::PseudoReferenceMaker' ); }

#set the default argument values
our %args = (
     vcf_file     => './data/single.dip.vcf'
   , out   => 'out.dummy'
   , bam   => 'out.dummy'
   , depth        => 4
   , depth_strand => 2
   , ratio        => 0.8
   , quality      => 50
   , map_quality  => 0
   , af1          => 0.95
   , af1_complement => 0.05
   , ci95         => 0.0
   , strand_bias  => 0.001
   , base_quality_bias => 0.0
   , map_bias     => 0.001
   , tail_bias    => 0.001
);

my $object = Pathogens::Variant::Utils::PseudoReferenceMaker->new(arguments => \%args );
isa_ok ($object, 'Pathogens::Variant::Utils::PseudoReferenceMaker');

unlink('out.dummy');
done_testing;
