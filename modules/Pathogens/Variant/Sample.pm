=head1 NAME

Pathogens::Variant::Sample  - This object can old an organism/lane/library etc...

=head1 SYNOPSIS

my $object = Pathogens::Variant::Sample->new;
$object->name('A_friendly_human');
$object->homozygous(1);
$object->genotype('AA');

=cut

package Pathogens::Variant::Sample;
use Moose;
extends 'Pathogens::Variant::Root';

has 'name' => (is => 'rw', isa => 'Str', required => 1);
has 'homozygous' => (is => 'rw', lazy => 1, default => '-1', isa => 'Bool');
has 'genotype' => (is => 'rw', lazy=> 1, default => '-1', isa => 'Str');

has 'alleles' => (
    traits  => ['Array'],
    is      => 'ro',
    isa     => 'ArrayRef[Pathogens::Variant::Allele]',
    default => sub { [] },
    lazy    => 1,
    handles => {
        get_alleles    => 'elements', #returns an array, NOT a ref
        add_alleles    => 'push', 
        clear_alleles  => 'clear', 
        has_alleles    => 'count',
    },
);

__PACKAGE__->meta->make_immutable;

1;
