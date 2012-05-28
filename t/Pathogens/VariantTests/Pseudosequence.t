use strict;
use warnings;
use Data::Dumper;
use Log::Log4perl qw(get_logger :levels); 
use Test::More;

###############################################
#log4perl initialisation settings
my $pseudoseq_logger = get_logger("Pathogens::Variant::Evaluator::Pseudosequence");
$pseudoseq_logger->level($FATAL);
my $layout = Log::Log4perl::Layout::PatternLayout->new("%d %p> %F{1}:%L %M - %m%n");
my $appender = Log::Log4perl::Appender->new("Log::Dispatch::Screen");
$appender->layout($layout);
$pseudoseq_logger->add_appender($appender);
###############################################



BEGIN { use_ok( 'Pathogens::Variant::Evaluator::Pseudosequence' ); }
use Pathogens::Variant::EvaluationReporter;
use Pathogens::Variant::Event::Snp;

my $object = Pathogens::Variant::Evaluator::Pseudosequence->new (reporter => Pathogens::Variant::EvaluationReporter->new);
isa_ok ($object, 'Pathogens::Variant::Evaluator::Pseudosequence');


my $event = Pathogens::Variant::Event::Snp->new(
   chromosome => 'FN433596'
 , position => '75182'
 , id => '.'
 , reference_allele => 'A'
 , alternative_allele => 'G,T'
 , quality => 222 
 , filter => '.' 
 , info => 'DP=138;VDB=0.0410;AF1=0.0001;AC1=2;DP4=7,8,56,67;MQ=45;FQ=-238;PV4=1,0.2,1,0.42'

);

my $reporter = Pathogens::Variant::EvaluationReporter->new;

my $evaluator = Pathogens::Variant::Evaluator::Pseudosequence->new (
   minimum_depth => 4
 , minimum_depth_strand => 2
 , minimum_ratio => 0.8
 , minimum_quality => 50
 , minimum_map_quality => 1
 , minimum_af1 => 0.95
 , minimum_ci95 =>0.0
 , minimum_strand_bias => 0.001
 , minimum_base_quality_bias => 0.0
 , minimum_map_bias => 0.001
 , minimum_tail_bias => 0.001
 , reporter => $reporter
);

$evaluator->evaluate($event);
is($event->passed_evaluation, 0, "Evaluation filters the snp out (due to low af1) as expected");

use Pathogens::Variant::Iterator::Vcf;
#this will traverse the VCF file and return Pathogens::Variant::Event::* objects
my $iterator = Pathogens::Variant::Iterator::Vcf->new(vcf_file => './data/single.dip.vcf');
while( $event = $iterator->next_event() ) {
    $evaluator->evaluate($event);
    is($event->passed_evaluation, 0, "Evaluation filters the indel out (no indels allowed) as expected");
}
$iterator->close_vcf();


$event = Pathogens::Variant::Event::Snp->new(
   chromosome => 'FN433596'
 , position => '75182'
 , id => '.'
 , reference_allele => 'A'
 , alternative_allele => 'G'
 , quality => 222 
 , filter => '.' 
 , info => 'DP=138;VDB=0.0410;AF1=0.99;AC1=2;DP4=7,8,56,67;MQ=45;FQ=-238;PV4=1,0.2,1,0.42'
);

$evaluator->evaluate($event);
is($event->passed_evaluation, 1, "Evaluation accepts the snp as expected");


done_testing();
