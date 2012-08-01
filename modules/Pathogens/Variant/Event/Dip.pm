=head1 NAME

Pathogens::Variant::Event::Dip  - A dip (deletion/insertion polymorphism) class

=head1 SYNOPSIS

my $object = Pathogens::Variant::Event::Dip->new (
                                              chromosome => 'test'
                                            , position => '1'
                                            , id => '.'
                                            , reference_allele => 'A'
                                            , alternative_allele => 'CGT'     #apparently an insertion relative to the reference.
                                            , quality => 50 
                                            , filter => 'something' 
                                            , info => 'INDEL;DP=179;VDB=0.0490;AF1=1;AC1=2;DP4=0,0,80,67;MQ=60;FQ=-290'
                                            );

=cut

package Pathogens::Variant::Event::Dip;
use Moose;
extends 'Pathogens::Variant::Event';

use namespace::autoclean;






has 'insertion_length'       => ( is => 'rw', isa => 'Int', lazy => 1, default => -1 );
has 'deletion_length'        => ( is => 'rw', isa => 'Int', lazy => 1, default => -1 );

__PACKAGE__->meta->make_immutable;

1;
