=head1 NAME

IntergenicRegions.pm   - Represents a IntergenicRegions from a GFF file

=head1 SYNOPSIS

use Pathogens::RNASeq::IntergenicRegions;
my $intergenic_regions = Pathogens::RNASeq::IntergenicRegions->new(
  features => $features,
  window_margin => 50,
  );
$intergenic_regions->intergenic_features;

=cut
package Pathogens::RNASeq::IntergenicRegions;
use Moose;
use Pathogens::RNASeq::IntergenicFeature;

has 'features'            => ( is => 'rw', isa => 'HashRef', required => 1 );
has 'window_margin'       => ( is => 'rw', isa => 'Int',     required => 1 );
has 'minimum_size'        => ( is => 'rw', isa => 'Int',     default  => 10 );
has 'sequence_lengths'    => ( is => 'rw', isa => 'HashRef', required => 1 );

has '_exons'               => ( is => 'rw', isa => 'HashRef', lazy_build => 1 );
has 'intergenic_features' => ( is => 'rw',  isa => 'HashRef', lazy_build => 1 );

sub _build__exons
{
  my ($self) = @_;
  my $exons;
  for my $feature (keys %{$self->features})
  {
    for my $exon (@{$self->features->{$feature}->exons})
    {
      $exons->{$self->features->{$feature}->seq_id}->{$exon->[0]} = $exon->[1];
    }
  }
  return $exons;
}

sub _build_intergenic_features
{
  my ($self) = @_;
  my $intergenic_features;
  
  for my $seq_id (keys %{$self->_exons})
  {
    my $previous_feature_start =0;
    my $previous_feature_end =  0 - $self->window_margin;
    
    for my $gene_start (sort keys %{$self->_exons->{$seq_id}})
    {
      my $intergenic_feature = $self->_create_intergenic_feature( $previous_feature_end,$gene_start,$seq_id );
      if(defined($intergenic_feature))
      {
        $intergenic_features->{$intergenic_feature->gene_id} = $intergenic_feature;
      }
      
      $previous_feature_start = $gene_start;
      $previous_feature_end = $self->_exons->{$seq_id}->{$gene_start};
    }
    
    # create an intergenic feature between the last real feature and the end of the sequence
    my $last_intergenic_feature = $self->_create_intergenic_feature( $previous_feature_end, $self->sequence_lengths->{$seq_id} + ($self->sequence_lengths->{$seq_id} - $self->calculate_intergenic_end($self->sequence_lengths->{$seq_id}))   ,$seq_id );
    if(defined($last_intergenic_feature))
    {
      $intergenic_features->{$last_intergenic_feature->gene_id} = $last_intergenic_feature;
    }
    
  }
  
  return $intergenic_features;
}

sub _create_intergenic_feature
{
  my ($self, $previous_feature_end,$gene_start,$seq_id ) = @_;
  my $intergenic_start = $self->calculate_intergenic_start($previous_feature_end);
  my $intergenic_end   = $self->calculate_intergenic_end($gene_start);
  
  if( $intergenic_end - $intergenic_start > $self->minimum_size)
  {
    my $intergenic_feature =  Pathogens::RNASeq::IntergenicFeature->new(
      gene_start => $intergenic_start,
      gene_end   => $intergenic_end,
      seq_id     => $seq_id
      );
      return $intergenic_feature;
  }
  return;
}

sub calculate_intergenic_start
{
   my ($self, $end_of_last_feature) = @_;
   return $end_of_last_feature + $self->window_margin + 1;
}

sub calculate_intergenic_end
{
  my ($self, $start_of_next_feature) = @_;
  return $start_of_next_feature - ( $self->window_margin + 1);
}

1;