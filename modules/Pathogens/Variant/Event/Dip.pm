package Pathogens::Variant::Event::Dip;
use Moose;
extends 'Pathogens::Variant::Event';

use namespace::autoclean;






has 'insertion_length'       => ( is => 'rw', isa => 'Int', lazy => 1, default => -1 );
has 'deletion_length'        => ( is => 'rw', isa => 'Int', lazy => 1, default => -1 );


=head1 NAME

Pathogens::Variant::Event::Dip - A dip (deletion/insertion polymorphism) class

=head1 SYNOPSIS

  use Pathogens::Variant::Event::Dip;

=head1 DESCRIPTION

The Pathogens::Variant::Event::Dip can hold deletion/insertion (dip) specific information

=head1 USAGE

    my $dip = Pathogens::Variant::Event::Dip->new();

=cut

=head1 SUPPORT

=head1 AUTHOR

    Feyruz Yalcin
    CPAN ID: FYALCIN

=head1 SEE ALSO

L<Pathogens::Variant::Tutorial>, perl(1).

=cut

__PACKAGE__->meta->make_immutable;

1;
