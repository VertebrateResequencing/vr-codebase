=head1 NAME

GFF.pm   - Represents a GFF from a GFF file

=head1 SYNOPSIS

use Pathogens::RNASeq::GFF;
my $rna_seq_gff = RNASeq::GFF->new(
  filename => 'my_file.gff'
  );
  
$rna_seq_gff->features();

=cut
package Pathogens::RNASeq::GFF;
use Moose;
use Bio::Tools::GFF;
use Pathogens::RNASeq::Feature;

has 'filename'    => ( is => 'rw', isa => 'Str',             required   => 1 );
has 'features'    => ( is => 'rw', isa => 'HashRef',         lazy_build => 1 );
has '_gff_parser' => ( is => 'rw', isa => 'Bio::Tools::GFF', lazy_build => 1 );


sub _build__gff_parser
{
  my ($self) = @_;
  Bio::Tools::GFF->new(-gff_version => 3, -file => $self->filename);
}

sub _build_features
{
  my ($self) = @_;
  my %features;
  
  while( my $raw_feature = $self->_gff_parser->next_feature())
  {
      last unless defined($raw_feature); # No more features
      next unless $raw_feature->primary_tag eq 'CDS'; # Only CDS features.

      my $feature_object = Pathogens::RNASeq::Feature->new(raw_feature => $raw_feature);

      if(defined($features{$feature_object->gene_id}))
      {
        $features{$feature_object->gene_id}->add_discontinuous_feature($raw_feature);
      }
      else
      {
        $features{$feature_object->gene_id} = $feature_object;
      }
      
  }
  
  return \%features;
}

1;