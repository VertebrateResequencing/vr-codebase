package Pathogens::Variant::Evaluator::Pseudosequence;
use Moose;


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

has 'reporter' => ( is => 'rw', isa => 'Pathogens::Variant::EvaluationReporter', required => 1 );

sub passed_filters {

    my ($self, $event) = @_;
    if ( $self->_passed_vcf_info_field_filters($event) ) {
        return 1;
    } else {
        return 0;
    }    
}


sub _passed_vcf_info_field_filters {
    my ($self, $event) = @_;
    
    my %param;
    my $passed = 1;
    
    foreach my $param_value ( split(";", $event->info)) {
        my ($parameter, $value) = split('=', $param_value);
        $param{$parameter} = $value;
    }

    if ( exists $param{'MQ'} and $param{'MQ'} < $self->minimum_map_quality) {
        $passed = 0;
        $self->reporter->inc_failed_allele_frequency;
    }
    
    if ( exists $param{'AF1'} and $param{'AF1'} < $self->minimum_af1) {
        $passed = 0;
    }
    
    if ( exists $param{'PV4'} ) {
        my ( $strand_bias, $base_quality_bias, $map_bias, $tail_distance_bias ) = split( ',', $param{'PV4'} );

        if ( $strand_bias < $self->minimum_strand_bias ) {
            $passed = 0;
        }
        if ( $base_quality_bias < $self->minimum_base_quality_bias) {
            $passed = 0;
        }
        if ( $map_bias < $self->minimum_map_bias) {
            $passed = 0;
        }
        if ( $tail_distance_bias < $self->minimum_tail_bias) {
            $passed = 0;
        }
    }
}

sub _dp4_evaluation {

    my $dp4_ratio_forward_reference   = 0;
    my $dp4_ratio_reverse_reference   = 0;
    my $dp4_ratio_forward_alternative = 0;
    my $dp4_ratio_reverse_alternative = 0;

}

=head1 BUGS

=head1 SUPPORT

=head1 AUTHOR

=head1 SEE ALSO

=cut

__PACKAGE__->meta->make_immutable;

1;