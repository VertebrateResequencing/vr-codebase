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
has 'info'               => ( is => 'rw', isa => 'Str' ); #may vary per entry
has 'format'             => ( is => 'rw', isa => 'Str' ); 

#non-VCF attributes
has 'event_type'         => ( is => 'rw', isa => 'Str' );
has 'polymorphism'       => ( is => 'rw', isa => 'Str' );

has 'samples' => (
    traits  => ['Array'],
    is      => 'ro',
    isa     => 'ArrayRef[Pathogens::Variant::Sample]',
    default => sub { [] },
    lazy    => 1,
    handles => {
        get_samples      => 'elements', #returns an array, NOT a ref
        add_sample       => 'push',
        has_samples      => 'count',
        clear_samples    => 'clear', 
    },
);



=head1 NAME

Pathogens::Variant::Event - A Pathogens::Variant::Event object

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 USAGE

=head1 SUPPORT

=head1 AUTHOR

=head1 SEE ALSO

=cut

__PACKAGE__->meta->make_immutable;

1;
