=head1 NAME

Pathogens::Variant::Utils::EventManipulator - To manipulate Pathogens::Variant::Event* objects... Currently very simple manipulations.

=head1 SYNOPSIS

my $object = Pathogens::Variant::Utils::EventManipulator->new;
my $event = Pathogens::Variant::Event->new ( 
                                              chromosome => 'test'
                                            , position => '1'
                                            , id => '.'
                                            , reference_allele => 'A'
                                            , alternative_allele => 'C,A'
                                            , quality => 50               #could represent variant quality for example...
                                            , filter => 'bla' 
                                            , info => 'bla'
                                            );

$object->remove_secondary_alternative_heterozygous_alleles($event);
#$event->alternative_allele would now return a 'C' instead of a 'C,A'

=cut

package Pathogens::Variant::Utils::EventManipulator;

use Moose;
use namespace::autoclean;


sub remove_secondary_alternative_heterozygous_alleles {

    my ($self, $event) = @_;
    my ($first_allele) = split(',', $event->alternative_allele); #effectively, all other alleles are discarded
    $event->alternative_allele($first_allele); #overwriting the original value

}

__PACKAGE__->meta->make_immutable;
1;