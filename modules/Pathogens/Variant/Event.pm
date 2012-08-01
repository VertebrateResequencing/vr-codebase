=head1 NAME

Pathogens::Variant::Event  - A (non-)polymorphic event. This represents any kind of (non)-polymorphic event, i.e. can correspond to a row in a VCF file for example.

=head1 SYNOPSIS

my $object = Pathogens::Variant::Event->new ( 
                                              chromosome => 'test'
                                            , position => '1'
                                            , id => '.'
                                            , reference_allele => 'A'
                                            , alternative_allele => 'C'
                                            , quality => 50               #could represent variant quality for example...
                                            , filter => 'something' 
                                            , info => 'DP=179;VDB=0.0490;AF1=1;AC1=2;DP4=0,0,80,67;MQ=60;FQ=-290'
                                            );

=cut

package Pathogens::Variant::Event;
use Moose;
extends 'Pathogens::Variant::Root';
with 'Pathogens::Variant::Role::Evaluable'; #implicitly adds the "passed_evaluation" (default "false")

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

__PACKAGE__->meta->make_immutable;

1;
