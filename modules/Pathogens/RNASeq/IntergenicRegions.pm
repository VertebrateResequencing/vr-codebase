=head1 NAME

IntergenicRegions.pm   - Represents a IntergenicRegions from a GFF file

=head1 SYNOPSIS

use Pathogens::RNASeq::IntergenicRegions;
my $intergenic_regions = RNASeq::IntergenicRegions->new(
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
has 'minimum_size'        => ( is => 'rw', isa => 'Int',     default => 10 );

has '_exons'               => ( is => 'rw', isa => 'HashRef', lazy_build => 1 );
has 'intergenic_features' => ( is => 'rw',  isa => 'HashRef', lazy_build => 1 );

sub _build__exons
{
  my ($self) = @_;
  my %exons;
  for my $feature (keys %{$self->features})
  {
     $exons->{$feature->seq_id}->{$feature->gene_start} = $feature;
  }
  return \%exons;
}

sub _build_intergenic_features
{
  my ($self) = @_;
  my %intergenic_features;
  
  for my $seq_id (keys %{$self->exons})
  {
    my $previous_feature;
    for my $gene_start (sort keys %{$self->exons->{$seq_id}})
    {
      if(defined ($previous_feature))
      {
        if( ($gene_start - $self->window_margin ) - ($previous_feature->gene_end + $self->window_margin) > $self->minimum_size)
        {
          my $intergenic_feature =  Pathogens::RNASeq::IntergenicFeature->new(gene_start => ($previous_feature->gene_end + $self->window_margin) ,
            gene_end => ($gene_start - $self->window_margin ),
            seq_id => $seq_id
            ) ;

          $intergenic_features->{$intergenic_feature->gene_id} = $intergenic_feature;
        }
      }
      $previous_feature = $self->exons->{$seq_id}->{$gene_start};
    }
  }
  
  return \%intergenic_features;
}

1;