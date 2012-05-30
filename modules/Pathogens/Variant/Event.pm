package Pathogens::Variant::Event;
use Moose;
extends 'Pathogens::Variant::Root';
with 'Pathogens::Variant::Role::Evaluable'; #implicitly adds the "passed_evaluation" attribute to this class, with the default "false"

use namespace::autoclean;






#VCF specific 8 fields + meta_lines
has 'chromosome'         => ( is => 'rw', isa => 'Str', lazy => 1, default => 'NoValueSet' );
has 'position'           => ( is => 'rw', isa => 'Int', lazy => 1, default => -1 );
has 'id'                 => ( is => 'rw', isa => 'Str', lazy => 1, default => 'NoValueSet' );
has 'reference_allele'   => ( is => 'rw', isa => 'Str', lazy => 1, default => 'NoValueSet' );
has 'alternative_allele' => ( is => 'rw', isa => 'Str', lazy => 1, default => 'NoValueSet' );
has 'quality'            => ( is => 'rw', isa => 'Num', lazy => 1, default => -1);
has 'filter'             => ( is => 'rw', isa => 'Str', lazy => 1, default => 'NoValueSet' );
has 'info'               => ( is => 'rw', isa => 'Str', lazy => 1, default => 'NoValueSet' ); #may vary per entry
has 'format'             => ( is => 'rw', isa => 'Str', lazy => 1, default => 'NoValueSet' ); 

#non-VCF attributes
has 'type'               => ( is => 'rw', isa => 'Str', lazy => 1, default => 'NoValueSet');
has 'polymorphic'        => ( is => 'rw', isa => 'Bool', lazy => 1, default => 1 );
has 'polymorphism'       => ( is => 'rw', isa => 'Str', lazy => 1, default => '[X/Y]' );

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
