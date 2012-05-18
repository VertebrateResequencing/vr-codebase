package Pathogens::Variant::Evaluator;
use Moose;

use Pathogens::Variant::Exception;

use namespace::autoclean;

sub evaluate {
    throw Pathogens::Variant::Exception::NotImplemented; 
}

=head1 BUGS

=head1 SUPPORT

=head1 AUTHOR

=cut

__PACKAGE__->meta->make_immutable;

1;
