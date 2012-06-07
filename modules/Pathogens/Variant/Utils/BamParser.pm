package Pathogens::Variant::Utils::BamParser;
use Moose;
extends 'Pathogens::Variant::Root';
use namespace::autoclean;

has 'bam' => (is  => 'rw', isa => 'Str', required => 1, trigger => \&_populate_chromosome_size_hash);

has 'chromsome_size' => (
      traits    => ['Hash'],
      is        => 'ro',
      isa       => 'HashRef[Str]',
      default   => sub { {} },
      handles   => {
          set_chromosome_size  => 'set',
          get_chromosome_size  => 'get',
          get_chromosome_names => 'keys'
      }
);

sub _populate_chromosome_size_hash {
    my ($self) = @_;

    #TO BE IMPLEMENTED!!!
    #(for test purposes)
    $self->set_chromosome_size('FN433596' => 3043210);
    $self->set_chromosome_size('pTW20_1' => 29585);
    $self->set_chromosome_size('pTW20_2' => 3011);

}
__PACKAGE__->meta->make_immutable;

1;