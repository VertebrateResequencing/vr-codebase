package Pathogens::Variant::Evaluator::Pseudosequence;
use Moose;
use Switch;

use namespace::autoclean;




has 'depth'        => ( is => 'rw', isa => 'Int', default => 4 );
has 'depth_strand' => ( is => 'rw', isa => 'Int', default => 2 );
has 'ratio'        => ( is => 'rw', isa => 'Num', default => 0.8 );
has 'quality'      => ( is => 'rw', isa => 'Int', default => 0 );
has 'map_quality'  => ( is => 'rw', isa => 'Int', default => 1 );
has 'af1'          => ( is => 'rw', isa => 'Num', required => 1 );
has 'ci95'         => ( is => 'rw', isa => 'Num', required => 1 );
has 'strand_bias'  => ( is => 'rw', isa => 'Num', required => 1 );
has 'quality_bias' => ( is => 'rw', isa => 'Num', required => 1 );
has 'map_bias'     => ( is => 'rw', isa => 'Num', required => 1 );
has 'tail_bias'    => ( is => 'rw', isa => 'Num', required => 1 );

sub evaluate {

    my ($self, $event) = @_;
    if ( $self->_passed_vcf_info_field_related_filters($event) ) {
        
    } else {
        #failed
    }    
}


sub _passed_vcf_info_field_related_filters {
    my ($self, $event) = @_;
    foreach my $param_value ( split(";", $event->info)) {
        my ($param, $value) = split('=', $param_value);

        switch ($param) {
            case 'MQ'     {}
            case 'DP4'    {}
            case 'AF1'    {}
            case 'PV4'    {}
        }
    }
}

sub _dp4_evaluation {
    
}

=head1 BUGS

=head1 SUPPORT

=head1 AUTHOR

=head1 SEE ALSO

=cut

__PACKAGE__->meta->make_immutable;

1;
