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