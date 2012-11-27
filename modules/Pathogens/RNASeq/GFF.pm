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
use Bio::Tools::NonStandardGFF;
use Pathogens::RNASeq::Feature;

has 'filename'          => ( is => 'rw', isa => 'Str',             required   => 1 );
has 'features'          => ( is => 'rw', isa => 'HashRef',         lazy_build => 1 );
has '_gff_parser'       => ( is => 'rw', isa => 'Bio::Tools::NonStandardGFF', lazy_build => 1 );
has 'sequence_lengths' => ( is => 'rw', isa => 'HashRef',         lazy_build => 1 );

sub _build__gff_parser
{
  my ($self) = @_;
  Bio::Tools::NonStandardGFF->new(-gff_version => 3, -file => $self->filename);
}

sub _build_features
{
  my ($self) = @_;
  my %features;
  
  while( my $raw_feature = $self->_gff_parser->next_feature())
  {
      last unless defined($raw_feature); # No more features
      next if !($raw_feature->primary_tag eq 'CDS' ||   $raw_feature->primary_tag eq 'polypeptide');
        
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

# create a hash with sequence names and the length of the sequence
sub _build_sequence_lengths
{
	my ($self) = @_;
	my %sequence_lengths;
	while(my $sequence_region = $self->_gff_parser->next_segment())
	{
		$sequence_lengths{$sequence_region->id()} = $sequence_region->end() + 1 - $sequence_region->start();
	}
	return \%sequence_lengths;
}

sub sorted_gene_ids
{
  my ($self) = @_;
  my @sorted_ids = sort _idsort keys(%{$self->features});
  return \@sorted_ids;
}

sub _idsort
{
  my @A = split(/\./,$a,2);
  my @B = split(/\./,$b,2);

  $A[0] cmp $B[0] || $A[1] <=> $B[1];
}

1;