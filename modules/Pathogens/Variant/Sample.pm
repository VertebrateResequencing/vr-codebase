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
 
=head1 NAME

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 USAGE

=head1 BUGS

=head1 SUPPORT

=head1 AUTHOR

    Feyruz Yalcin
    CPAN ID: FYALCIN

=head1 SEE ALSO

=cut

__PACKAGE__->meta->make_immutable;

1;
