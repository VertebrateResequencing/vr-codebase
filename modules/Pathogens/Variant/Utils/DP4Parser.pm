package Pathogens::Variant::Utils::DP4Parser;
use Moose;
use namespace::autoclean;

has 'count_reference_forward_bases' => (is => 'rw', isa => 'Int', default => -1);
has 'count_alternative_forward_bases' => (is => 'rw', isa => 'Int', default => -1);
has 'count_reference_reverse_bases' => (is => 'rw', isa => 'Int', default => -1);
has 'count_alternative_reverse_bases' => (is => 'rw', isa => 'Int', default => -1);
has 'count_alternative_bases' => (is => 'rw', isa => 'Int', default => -1);
has 'count_referecence_bases' => (is => 'rw', isa => 'Num', default => -1);

has 'ratio_forward_reference_bases' => (is => 'rw', isa => 'Num', default => -1);
has 'ratio_forward_alternative_bases' => (is => 'rw', isa => 'Num', default => -1);
has 'ratio_reverse_reference_bases' => (is => 'rw', isa => 'Num', default => -1);
has 'ratio_reverse_alternative_bases' => (is => 'rw', isa => 'Num', default => -1);

sub parse {

    my ($self, $dp4string) = @_;

    #format of the $dp4string:
    #ref-forward bases, ref-reverse bases, alt-forward bases, alt-reverse bases
    my (   $ref_forward                            #1st
         , $ref_reverse                            #2nd
         , $alt_forward                            #3rd
         , $alt_reverse) = split(',', $dp4string); #4th


    #set the counts:
    $self->count_reference_forward_bases($ref_forward);
    $self->count_alternative_forward_bases($alt_forward);
    $self->count_reference_reverse_bases($ref_reverse);
    $self->count_alternative_reverse_bases($alt_reverse);
    
    $self->count_referecence_bases($ref_forward + $ref_reverse);
    $self->count_alternative_bases($alt_forward + $alt_reverse);


    #set the ratios
    if (($ref_forward + $alt_forward) > 0 ) {

        $self->ratio_forward_reference_bases($ref_forward/($ref_forward + $alt_forward));
        $self->ratio_forward_alternative_bases($alt_forward/($ref_forward + $alt_forward));

    } else {
        
        $self->ratio_forward_reference_bases(0);
        $self->ratio_forward_alternative_bases(0);

    }

    if ( ($ref_reverse+$alt_reverse) > 0 ) {

        $self->ratio_reverse_reference_bases($ref_reverse/($ref_reverse + $alt_reverse));
        $self->ratio_reverse_alternative_bases($alt_reverse/($ref_reverse + $alt_reverse));

    } else {

        $self->ratio_reverse_reference_bases(0);
        $self->ratio_reverse_alternative_bases(0);

    }

}

__PACKAGE__->meta->make_immutable;
1;