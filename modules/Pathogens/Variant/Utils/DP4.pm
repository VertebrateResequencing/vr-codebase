package DP4;
use Moose;
use namespace::autoclean;


has 'forward_reference_strand_depth'   => (is => 'rw', isa => 'Int', default => -1);
has 'forward_alternative_strand_depth' => (is => 'rw', isa => 'Int', default => -1);

has 'reverse_reference_strand_depth'   => (is => 'rw', isa => 'Int', default => -1);
has 'reverse_alternative_strand_depth' => (is => 'rw', isa => 'Int', default => -1);


has 'ratio_forward_reference_bases'    => (is => 'rw', isa => 'Num', default => -1);
has 'ratio_forward_alternative_bases'  => (is => 'rw', isa => 'Num', default => -1);

has 'ratio_reverse_reference_bases'    => (is => 'rw', isa => 'Num', default => -1);
has 'ratio_reverse_alternative_bases'  => (is => 'rw', isa => 'Num', default => -1);

__PACKAGE__->meta->make_immutable;