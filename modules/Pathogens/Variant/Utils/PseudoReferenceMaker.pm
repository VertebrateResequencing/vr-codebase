package Pathogens::Variant::Utils::PseudoReferenceMaker;

use Moose;
use Pathogens::Variant::Iterator::Vcf;
use Pathogens::Variant::Evaluator::Pseudosequence;
use Pathogens::Variant::EvaluationReporter;
use namespace::autoclean;

has 'run_arguments' => (is => 'rw', isa => 'HashRef', required => 1, trigger => \&_initialise);
has '_iterator'     => (is => 'rw', isa => 'Pathogens::Variant::Iterator::Vcf', init_arg => undef );
has '_reporter'     => (is => 'rw', isa => 'Pathogens::Variant::EvaluationReporter', init_arg => undef );
has '_evaluator'    => (is => 'rw', isa => 'Pathogens::Variant::Evaluator::Pseudosequence', init_arg => undef );

has '_reference_lengths' => (
      traits  => ['Hash'],
      is      => 'rw',
      isa     => 'HashRef[Str]',
      default => sub { {} },
      handles => {
          exists_in_reference_lengths => 'exists',
          ids_in_reference_lengths    => 'keys',
          get_reference_length        => 'get',
          set_reference_length        => 'set'
      }
);
 
sub create_pseudo_reference {
    
    my ($self) = @_;
    
    my $c = 0;
    while( my $event = $self->_iterator->next_event() ) {
        
        $c++; print $c, "\n" unless $c % 10000;
        $self->_evaluator->evaluate($event); #note: evaluator will modify heterozygous calls
        if ($event->passed_evaluation) {
            if (not $event->polymorphic) {

            } else {

            }
        }
    }
    $self->_iterator->close_vcf();
}

#initialise all the objects/settings for this run
sub _initialise {
    
    my ($self) = @_;
    
    #complement af1 is needed for non variant site evaluations
    $self->run_arguments->{af1complement} = 1 - $self->run_arguments->{af1};
    
    #total depth has to be at least twice larger than strand-depth
    my $doubled_arg_D_depth_strand = $self->run_arguments->{depth_strand} * 2;
    if ( $self->run_arguments->{depth} < $doubled_arg_D_depth_strand ) {
        print "'-d' (depth) must be larger than '-D' (depth_strand)! Automatically increasing it to " .$doubled_arg_D_depth_strand ."\n";
        $self->run_arguments->{depth} = $doubled_arg_D_depth_strand;
    }
    
    #this will create the vcf file traverser
    my $iterator = Pathogens::Variant::Iterator::Vcf->new(vcf_file => $self->run_arguments->{vcf_file});
    
    #reporter will be used by the evaluator object below to keep counters on evaluation statistics
    my $reporter = Pathogens::Variant::EvaluationReporter->new;
    
    #create an evaluator object to filter bad calls based on user criteria
    my $evaluator = Pathogens::Variant::Evaluator::Pseudosequence->new(
          minimum_depth             => $self->run_arguments->{depth}
        , minimum_depth_strand      => $self->run_arguments->{depth_strand}
        , minimum_ratio             => $self->run_arguments->{ratio}
        , minimum_quality           => $self->run_arguments->{quality}
        , minimum_map_quality       => $self->run_arguments->{map_quality}
        , minimum_af1               => $self->run_arguments->{af1}
        , af1_complement            => $self->run_arguments->{af1_complement}
        , minimum_ci95              => $self->run_arguments->{ci95}
        , minimum_strand_bias       => $self->run_arguments->{strand_bias}
        , minimum_base_quality_bias => $self->run_arguments->{base_quality_bias}
        , minimum_map_bias          => $self->run_arguments->{map_bias}
        , minimum_tail_bias         => $self->run_arguments->{tail_bias}
        , reporter                  => $reporter
    );
    
    #set the objects into appropriate attributes for later access withing this class
    $self->_evaluator($evaluator);
    $self->_iterator($iterator);
    $self->_reporter($reporter);
    
}

sub report_run_statistics {
    my ($self, $fhd) = @_;
    
}

__PACKAGE__->meta->make_immutable;
1;