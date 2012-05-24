package Pathogens::Variant::Role::Evaluable;


use Moose::Role;
use namespace::autoclean;


has 'passed_evaluation' => (
      is  => 'rw',
      isa => 'Bool',
);

1;