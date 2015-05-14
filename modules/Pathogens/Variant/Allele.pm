=head1 NAME

Pathogens::Variant::Allele  - An allele class

=head1 SYNOPSIS

my $object = Pathogens::Variant::Allele->new;
$object->coverage(10);
$object->base('A');
$object->quality(30);

#etc...

=cut

package Pathogens::Variant::Allele;

use Moose;
extends 'Pathogens::Variant::Root';
use namespace::autoclean;








has 'coverage'   => ( is => 'rw', isa => 'Int');
has 'base'       => ( is => 'rw', isa => 'Str');
has 'quality'    => ( is => 'rw', isa => 'Num');
has 'coverage'   => ( is => 'rw', isa => 'Int');

__PACKAGE__->meta->make_immutable;

1;
