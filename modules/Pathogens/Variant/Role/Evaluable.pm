=head1 NAME

Pathogens::Variant::Role::Evaluable - A Moose role for variant objects that can be evaluated, i.e. subjected to tests according to some criteria. For example a Pathogens::Variant::Event::Snp object can be "Evaluable", i.e. one can evaluate such an object by looking if fulfills certain thresholds to be considered a real Snp.

=head1 SYNOPSIS

=cut

package Pathogens::Variant::Role::Evaluable;

use Moose::Role;
use namespace::autoclean;


has 'passed_evaluation' => (
      is  => 'rw',
      isa => 'Bool',
      default => undef,
      trigger => \&_set_evaluated
);

has 'was_evaluated' => (
      is  => 'rw',
      isa => 'Bool',
      default => undef
);

sub _set_evaluated {
    my ($self) = @_;
    $self->was_evaluated(1);

}

1;