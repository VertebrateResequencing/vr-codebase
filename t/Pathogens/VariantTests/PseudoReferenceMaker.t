use strict;
use warnings;

use Test::More 'no_plan';
use Data::Dumper;
use Log::Log4perl qw(get_logger :levels); 

###############################################
#log4perl initialisation settings
my $pseudoseq_logger = get_logger("Pathogens");
$pseudoseq_logger->level($FATAL);
my $layout = Log::Log4perl::Layout::PatternLayout->new("%d %p> %F{1}:%L %M - %m%n");
my $appender = Log::Log4perl::Appender->new("Log::Dispatch::Screen");
$appender->layout($layout);
$pseudoseq_logger->add_appender($appender);
###############################################



BEGIN { use_ok( 'Pathogens::Variant::Utils::PseudoReferenceMaker' ); }

#set the default argument values
our %args = (
     vcf_file     => './data/single.dip.vcf'
   , out   => './data/temporary.test.out'
   , bam   => './data/test.bam'
   , reference => './data/reference.snpAt10.insertionAt30.deletionAt50.fasta'
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

my $name_of_created_file = $object->_generate_vcf_file_with_all_reference_sites;
isnt($name_of_created_file, '', "Managed to create a VCF file with all sites");

unlink($name_of_created_file);
 
#some modules create a dir called "_Inline" in the current dir
#delete it after the test
