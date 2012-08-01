=head1 NAME

Pathogens::Variant::Utils::BamParser  - Parses a BAM file to get the reference names and their lengths. 

=head1 SYNOPSIS

#requires SAMTOOLS shell environment variable to point to samtools compiled dir

my $object = Pathogens::Variant::Utils::BamParser->new();
$object->fetch_chromosome_size_into_hash("$dir/data/test.bam");

foreach my $name ($object->get_chromosome_names) {
    my $size_in_bp = $object->get_chromosome_size($name);  
}

=cut

package Pathogens::Variant::Utils::BamParser;
use Moose;
extends 'Pathogens::Variant::Root';
use VertRes::Parser::bam;

use namespace::autoclean;

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

sub fetch_chromosome_size_into_hash {
    my ($self, $bam) = @_;
   
    my $parser = VertRes::Parser::bam->new(file => $bam);
    my %sequence_info = $parser->sequence_info();
    
    foreach my $reference_name (keys %sequence_info) {
        
        #the key will be $reference_name, and the value will be its length as shown in
        #the bam file header with "LN"
        
        $self->set_chromosome_size( $reference_name => $sequence_info{$reference_name}->{"LN"} );
        
    }
    
}
__PACKAGE__->meta->make_immutable;

1;