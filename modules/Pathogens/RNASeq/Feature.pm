=head1 NAME

Feature.pm   - Represents a Feature from a GFF file

=head1 SYNOPSIS

use Pathogens::RNASeq::Feature;
my $file_meta_data_container = RNASeq::Feature->new(
  raw_feature => $feature
  );

=cut
package Pathogens::RNASeq::Feature;
use Moose;

has 'raw_feature'   => ( is => 'rw', isa => 'Bio::SeqFeature::Generic', required   => 1 );

has 'gene_id'       => ( is => 'rw', isa => 'Str',                 lazy_build => 1 );
has 'seq_id'        => ( is => 'rw', isa => 'Str',                 lazy_build => 1 );
has 'gene_strand'   => ( is => 'rw', isa => 'Int',                 lazy_build => 1 );
has 'gene_start'    => ( is => 'rw', isa => 'Int',                 lazy_build => 1 );
has 'gene_end'      => ( is => 'rw', isa => 'Int',                 lazy_build => 1 );
has 'exon_length'   => ( is => 'rw', isa => 'Int',                 lazy_build => 1 );
has 'exons'         => ( is => 'rw', isa => 'ArrayRef',            lazy_build => 1 );

sub _build_exons
{
  my ($self) = @_;
  my @exons;
  push @exons, [$self->gene_start, $self->gene_end];

  return \@exons;
}

sub _find_feature_id
{
  my ($self) = @_;
  my $gene_id;
  my @junk;
  my @tag_names = ('ID', 'locus_tag', 'Parent');
  
  for my $tag_name (@tag_names)
  {
    if($self->raw_feature->has_tag($tag_name))
    {
      ($gene_id, @junk) = $self->raw_feature->get_tag_values($tag_name);
      return $gene_id;
    }
  }

  return $gene_id;
}

sub _build_gene_id
{
  my ($self) = @_;

  my $gene_id;
  my @junk;
  
  $gene_id = $self->_find_feature_id();
  
  if(! defined($gene_id))
  {
    $gene_id = join("_",($self->seq_id, $self->gene_start, $self->gene_end ));
  }
  $gene_id =~ s/^"|"$//g;

  return $gene_id;
}

sub _build_seq_id
{
  my ($self) = @_;
  $self->raw_feature->seq_id();
}

sub _build_gene_strand
{
  my ($self) = @_;
  $self->raw_feature->strand;
}

sub _build_gene_start
{
  my ($self) = @_;
  $self->raw_feature->start;
}

sub _build_gene_end
{
  my ($self) = @_;
  $self->raw_feature->end;
}

sub _build_exon_length
{
  my ($self) = @_;
  ($self->gene_end - $self->gene_start + 1);
}

sub _filter_out_parent_features_from_exon_list
{
  my ($self) = @_;
  return if(scalar(@{$self->exons}) <= 1);

  my @filtered_exons;

  for my $exon_coords(@{$self->exons})
  {
    if(!($exon_coords->[0] == $self->gene_start && $exon_coords->[1] == $self->gene_end))
    {
      push(@filtered_exons, $exon_coords);
    }
  }
  $self->exons(\@filtered_exons);
}


sub add_discontinuous_feature
{
  my ($self,$raw_feature) = @_;

  my $gene_start = ($raw_feature->start < $self->gene_start) ? $raw_feature->start : $self->gene_start;
  my $gene_end = ($raw_feature->end   > $self->gene_end) ? $raw_feature->end : $self->gene_end;

  my $exon_length = ($raw_feature->end - $raw_feature->start + 1) + $self->exon_length;

  $self->gene_start($gene_start);
  $self->gene_end($gene_end);
  $self->exon_length($exon_length);

  push @{$self->exons}, [$raw_feature->start, $raw_feature->end ];
  $self->_filter_out_parent_features_from_exon_list();
  
  return;
}
1;

