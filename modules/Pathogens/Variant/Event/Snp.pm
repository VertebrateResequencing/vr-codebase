=head1 NAME

Pathogens::Variant::Event::Snp  - A Snp class

=head1 SYNOPSIS

my $object = Pathogens::Variant::Event::Snp->new (
                                              chromosome => 'test'
                                            , position => '1'
                                            , id => '.'
                                            , reference_allele => 'A'
                                            , alternative_allele => 'C'     #apparently C/A polymorphism
                                            , quality => 50 
                                            , filter => 'something' 
                                            , info => 'DP=179;VDB=0.0490;AF1=1;AC1=2;DP4=0,0,80,67;MQ=60;FQ=-290'
                                            );

=cut

package Pathogens::Variant::Event::Snp;
use Moose;
extends 'Pathogens::Variant::Event';

use namespace::autoclean;





has "transversion" => (is => 'rw', isa => 'Bool', lazy => 1, default => 0 );
has "transition"   => (is => 'rw', isa => 'Bool', lazy => 1, default => 0 );
has 'synonymous'   => (is => 'rw', isa => 'Bool', lazy => 1, default => 0 );

__PACKAGE__->meta->make_immutable;

1;
