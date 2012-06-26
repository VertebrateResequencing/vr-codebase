=head1 NAME

Pathogens::Variant::Role::Evaluator - A Moose role for variant objects that evaluate others.

=head1 SYNOPSIS

=cut

package Pathogens::Variant::Role::Evaluator;

use Moose::Role;
use namespace::autoclean;

requires 'evaluate';

1;
