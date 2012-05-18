package Pathogens::Variant::Event::Snp;
use Moose;
extends 'Pathogens::Variant::Event';

use namespace::autoclean;








=head1 NAME

Pathogens::Variant::Event::Snp - A snp object

=head1 SYNOPSIS

  use Pathogens::Variant::Event::Snp;

=head1 DESCRIPTION

The Pathogens::Variant::Event::Snp can hold snp specific information.

=head1 SUPPORT

=head1 AUTHOR

    Feyruz Yalcin
    CPAN ID: FYALCIN

=head1 SEE ALSO

L<Pathogens::Variant::Tutorial>, perl(1).

=cut

__PACKAGE__->meta->make_immutable;

1;
