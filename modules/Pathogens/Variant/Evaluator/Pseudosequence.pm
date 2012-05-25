package Pathogens::Variant::Evaluator::Pseudosequence;
use Moose;
extends 'Pathogens::Variant::Root';
with 'Pathogens::Variant::Role::Evaluator';#forces this class to implement "evaluate" subroutine

use Pathogens::Variant::Utils::DP4Parser;
use Log::Log4perl qw(get_logger);
use Data::Dumper;
use namespace::autoclean;




has 'minimum_depth'        => ( is => 'rw', isa => 'Int', default => 4 );
has 'minimum_depth_strand' => ( is => 'rw', isa => 'Int', default => 2 );
has 'minimum_ratio'        => ( is => 'rw', isa => 'Num', default => 0.8 );
has 'minimum_quality'      => ( is => 'rw', isa => 'Int', default => 0 );
has 'minimum_map_quality'  => ( is => 'rw', isa => 'Int', default => 1 );
has 'minimum_af1'          => ( is => 'rw', isa => 'Num', default => 0.95 );
has 'minimum_ci95'         => ( is => 'rw', isa => 'Num', default => 0.0 );
has 'minimum_strand_bias'  => ( is => 'rw', isa => 'Num', default => 0.001 );
has 'minimum_base_quality_bias' => ( is => 'rw', isa => 'Num', default => 0.0 );
has 'minimum_map_bias'     => ( is => 'rw', isa => 'Num', default => 0.001 );
has 'minimum_tail_bias'    => ( is => 'rw', isa => 'Num', default => 0.001 );

has 'reporter'             => ( is => 'rw', isa => 'Pathogens::Variant::EvaluationReporter', default => sub { return Pathogens::Variant::EvaluationReporter->new } );

has '_dp4_parser'          => ( is => 'ro', isa => 'Pathogens::Variant::Utils::DP4Parser', lazy => 1, default => sub { return Pathogens::Variant::Utils::DP4Parser->new } );

sub evaluate {

    my ($self, $event) = @_;
    my $logger = get_logger("Pathogens::Variant::Evaluator::Pseudosequence");

    $logger->debug("Evaluating event:...\n". Dumper $event);

    if ( $self->_evaluate_vcf_info_field_values($event) ) {
        $event->passed_evaluation(1);
        
        #log if in debug mode
        $logger->debug("Event dump after evaluation:...\n". Dumper $event);
        $logger->debug("Reporter dump after evaluation:...\n". Dumper $self->reporter);
        
        return 1;
    }
    
    $event->passed_evaluation(0); #i.e. the event failed to pass the filters

    #log if in debug mode
    $logger->debug("Event dump after evaluation:...\n". Dumper $event);
    $logger->debug("Reporter dump after evaluation:...\n". Dumper $self->reporter);
    
    return 0;
}

sub _evaluate_vcf_info_field_values {
    
    my ($self, $event) = @_;
    my $logger = get_logger("Pathogens::Variant::Evaluator::Pseudosequence");

    my %param;
    my $passed = 1;
   
    foreach my $param_value ( split(";", $event->info)) {
        my ($parameter, $value) = split('=', $param_value);
        $param{$parameter} = $value;
    }

    #log if in debug mode
    $logger->debug("Evaluating vcf info field values...\n" . Dumper %param);

    if ( exists $param{'MQ'} and $param{'MQ'} < $self->minimum_map_quality) {
        $passed = 0;
        $self->reporter->inc_map_quality_fail;
    }
    
    if ( exists $param{'AF1'} and $param{'AF1'} < $self->minimum_af1) {
        $passed = 0;
        $self->reporter->inc_af1_fail;
    }
    
    if ( exists $param{'PV4'} ) {
        my ( $strand_bias, $base_quality_bias, $map_bias, $tail_distance_bias ) = split( ',', $param{'PV4'} );

        if ( $strand_bias < $self->minimum_strand_bias ) {
            $passed = 0;
            $self->reporter->inc_strand_bias_fail;
        }
        if ( $base_quality_bias < $self->minimum_base_quality_bias) {
            $passed = 0;
            $self->reporter->inc_base_quality_bias_fail;
        }
        if ( $map_bias < $self->minimum_map_bias) {
            $passed = 0;
            $self->reporter->inc_map_bias_fail;
        }
        if ( $tail_distance_bias < $self->minimum_tail_bias) {
            $passed = 0;
            $self->reporter->inc_tail_bias_fail;
        }
    }

    if ( exists $param{'DP4'} ) {
        
        my $parsed = $self->_dp4_parser->parse( $param{'DP4'} );
 
        #depths for ref bases and alt bases
        if ( $parsed->{count_referecence_bases} < $self->minimum_depth) {
            
        }
        if ( $parsed->{count_alternative_bases} < $self->minimum_depth) {
            
        }

        #depths for ref/alt forward bases and ref/alt reverse bases
        if ( $parsed->{count_reference_forward_bases} < $self->minimum_depth_strand) {
            
        }
        if ( $parsed->{count_reference_reverse_bases} < $self->minimum_depth_strand) {
            
        }

        if ( $parsed->{count_alternative_forward_bases} < $self->minimum_depth_strand) {
            
        }
        if ( $parsed->{count_alternative_reverse_bases} < $self->minimum_depth_strand) {
            
        }

        #ratios for ref/alt forward bases and ref/alt reverse bases
        if ( $parsed->{ratio_forward_reference_bases} < $self->minimum_ratio) {
            
        }
        if ( $parsed->{ratio_forward_alternative_bases} < $self->minimum_ratio) {
            
        }
        if ( $parsed->{ratio_reverse_reference_bases} < $self->minimum_ratio) {
            
        }
        if ( $parsed->{ratio_reverse_alternative_bases} < $self->minimum_ratio) {
            
        }
        
    } else {
        
    }
}



=head1 BUGS

=head1 SUPPORT

=head1 AUTHOR

=head1 SEE ALSO

=cut

__PACKAGE__->meta->make_immutable;

1;