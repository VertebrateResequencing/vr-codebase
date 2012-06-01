package Pathogens::Variant::Utils::BamParser;
use Moose;
extends 'Pathogens::Variant::Root';
use namespace::autoclean;

has 'bam' => (is  => 'rw', isa => 'Str', required => 1);



sub fetch_hashref_chromosome_sizes {
    my ($self) = @_;
    
    #
    #IMPLEMENT A METHOD TO GET A HASHREF OF {REFERENCES => SIZE}
    #
    
    #for tests
    my %chromosome_sizes  = (
        FN433596 => 3043210,
        pTW20_1  => 29585,
        pTW20_2  => 3011,
    );
    
    return \%chromosome_sizes;
}
__PACKAGE__->meta->make_immutable;

1;