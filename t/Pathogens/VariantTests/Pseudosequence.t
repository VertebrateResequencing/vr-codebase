use strict;
use warnings;
use Data::Dumper;

use Test::More tests => 2;


BEGIN { use_ok( 'Pathogens::Variant::Evaluator::Pseudosequence' ); }
use Pathogens::Variant::EvaluationReporter;
use Pathogens::Variant::Event::Snp;

my $object = Pathogens::Variant::Evaluator::Pseudosequence->new (reporter => Pathogens::Variant::EvaluationReporter->new);
isa_ok ($object, 'Pathogens::Variant::Evaluator::Pseudosequence');


my $snp_event = Pathogens::Variant::Event::Snp->new(
   chromosome => 'FN433596'
 , position => '75182'
 , id => '.'
 , reference_allele => 'A'
 , alternative_allele => 'G'
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
 , minimum_map_quality => 0
 , minimum_af1 => 0.95
 , minimum_ci95 =>0.0
 , minimum_strand_bias => 0.001
 , minimum_base_quality_bias => 0.0
 , minimum_map_bias => 0.001
 , minimum_tail_bias => 0.001
 , reporter => $reporter
);

my $filter_status = $evaluator->passed_filters($snp_event);

