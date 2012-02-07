=head1 NAME

IntergenicFeature.pm   - Represents a IntergenicFeature from a GFF file

=head1 SYNOPSIS

use Pathogens::RNASeq::IntergenicFeature;
my $inter = RNASeq::IntergenicFeature->new(
  raw_feature => $feature
  );

=cut
package Pathogens::RNASeq::IntergenicFeature;
use Moose;

has 'gene_start'    => ( is => 'rw', isa => 'Int',                 required => 1 );
has 'gene_end'      => ( is => 'rw', isa => 'Int',                 required => 1 );
has 'seq_id'        => ( is => 'rw', isa => 'Str',                 required => 1 );
has 'gene_strand'   => ( is => 'rw', isa => 'Int',                 default => 1 );

has 'gene_id'       => ( is => 'rw', isa => 'Str',                 lazy_build => 1 );
has 'exon_length'   => ( is => 'rw', isa => 'Int',                 lazy_build => 1 );
has 'exons'         => ( is => 'rw', isa => 'ArrayRef',            lazy_build => 1 );

sub _build_gene_id 
{
   my ($self) = @_;
   return "".$self->seq_id."_intergenic_".$self->gene_start."_".$self->gene_end;
}

sub _build_exon_length
{
  my ($self) = @_;
  return $self->gene_end - $self->gene_start;
}

sub _build_exons
{
  my ($self) = @_;
  my @exons  = [];
  push(@exons, [$self->gene_start, $self->gene_end])
  return \@exons;
}
