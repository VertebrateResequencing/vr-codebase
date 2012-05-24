package Pathogens::Variant::Event;
use Moose;
extends 'Pathogens::Variant::Root';
with 'Pathogens::Variant::Role::Evaluable'; #implicitly adds the "passed_evaluation" attribute to this class, with the default "false"

use namespace::autoclean;






#VCF specific 8 fields + meta_lines
has 'chromosome'         => ( is => 'rw', isa => 'Str' );
has 'position'           => ( is => 'rw', isa => 'Int' );
has 'id'                 => ( is => 'rw', isa => 'Str' );
has 'reference_allele'   => ( is => 'rw', isa => 'Str' );
has 'alternative_allele' => ( is => 'rw', isa => 'Str' );
has 'quality'            => ( is => 'rw', isa => 'Num' );
has 'filter'             => ( is => 'rw', isa => 'Str' );
has 'info'               => ( is => 'rw', isa => 'Str' ); #vcf's info field (may vary per entry within the same file)
has 'format'             => ( is => 'rw', isa => 'Str' ); 
has 'meta_lines'         => ( is => 'rw', isa => 'ArrayRef[Str]' );

#Additional non-VCF fields
has 'event_type'         => ( is => 'rw', isa => 'Str' );
has 'polymorphism'       => ( is => 'rw', isa => 'Str' );
has 'samples'            => ( is => 'rw', isa => 'ArrayRef', default => sub {[]} );




sub add_sample {
    my ($self, @samples) = @_;
    push (@{$self->samples}, @samples);
}

=head1 NAME

Pathogens::Variant::Event - A Pathogens::Variant::Event object

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 USAGE

=head1 SUPPORT

=head1 AUTHOR

    Feyruz Yalcin
    CPAN ID: FYALCIN

=head1 SEE ALSO

=cut

__PACKAGE__->meta->make_immutable;

1;
