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