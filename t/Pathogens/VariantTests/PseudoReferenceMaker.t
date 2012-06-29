use strict;
use warnings;
use Data::Dumper;
use Log::Log4perl qw(get_logger :levels); 
use File::Basename; 

use Test::More tests => 6;



#get the directory path to the test file
my (undef, $dir) = fileparse($0);

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


my $samtools_environment_variable_exists = 1 if exists  $ENV{SAMTOOLS};
is($samtools_environment_variable_exists, 1, "SAMTOOLS environment is present found");

my $pseudo_reference_output = "$dir/data/temporary.test.out";

#set the default argument values
my %args = (
     out   => $pseudo_reference_output
   , lane_name  => "test_lane"
   , bam   => "$dir/data/test.bam"
   , reference => "$dir/data/reference.snpAt11.insertionAt31.deletionAt51.fasta"
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
unlink("$dir/data/$name_of_created_file");

$object->create_pseudo_reference;
my $file_exists = 1 if (-e $pseudo_reference_output);
is($file_exists, 1, "We have managed to create a pseudoreference file");

#we check if the output is identical to the one that Simon Harris' original script generated
open(FHD_NEW, "<$pseudo_reference_output");

#"simons.pseudo.reference.new.fasta" has been derived from simons ".mfa" file. The mfa was created by
#multiple_mappings_to_bam.py with options "-s -I -c". The .mfa was a multi-fasta file with the original 
#reference as the first, and followed by the pseudo-reference chromosomes.  To derive "simons.pseudo.reference.new.fasta"
#we manually remove the original reference from the mfa, and concatenate the chromosomes of the remaining pseudoreference
#into a single fasta entry. This file is then named "simons.pseudo.reference.new.fasta".
open(FHD_OLD, "<$dir/data/simons.pseudo.reference.fasta");

my $differing_line_count = 0;
while(<FHD_NEW>) {
    my $new_sequence = $_;
    my $simons_sequence = <FHD_OLD>;
    if ($new_sequence ne $simons_sequence) {
        $differing_line_count++;
    }
}
close FHD_NEW;
close FHD_OLD;

is($differing_line_count, 0, "Pseudoreference is identical to the one that would have been generated with Simon's 'multiple_mappings_to_bam.py' script if his is run the with -I (i.e. 'ignore indel positions') option...");
unlink($pseudo_reference_output);