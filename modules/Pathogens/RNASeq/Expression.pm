=head1 NAME

Expression.pm   - Find the expression when given an input aligned file and an annotation file

=head1 SYNOPSIS

use Pathogens::RNASeq::Expression;
my $expression_results = Pathogens::RNASeq::Expression->new(
  sequence_filename => 'my_aligned_sequence.bam',
  annotation_filename => 'my_annotation_file.gff'
  );

$expression_results->output_spreadsheet();

=cut
package Pathogens::RNASeq::Expression;
use Moose;
use Pathogens::RNASeq::SequenceFile;
use Pathogens::RNASeq::GFF;
use Pathogens::RNASeq::AlignmentSlice;
use Pathogens::RNASeq::ExpressionStatsSpreadsheet;

has 'sequence_filename'     => ( is => 'rw', isa => 'Str', required => 1 );
has 'annotation_filename'   => ( is => 'rw', isa => 'Str', required => 1 );
#optional filters
has 'filters'               => ( is => 'rw', isa => 'Maybe[HashRef]'     );

has '_sequence_file'        => ( is => 'rw', isa => 'Pathogens::RNASeq::SequenceFile',               lazy_build  => 1 );
has '_annotation_file'      => ( is => 'rw', isa => 'Pathogens::RNASeq::GFF',                        lazy_build  => 1 );
has '_results_spreadsheet'  => ( is => 'rw', isa => 'Pathogens::RNASeq::ExpressionStatsSpreadsheet', lazy_build  => 1 );
has '_expression_results'   => ( is => 'rw', isa => 'ArrayRef',                                      lazy_build  => 1 );


sub _build__sequence_file
{
  my ($self) = @_;
  Pathogens::RNASeq::SequenceFile->new(filename => $self->sequence_filename);
}

sub _build__annotation_file
{
  my ($self) = @_;
  Pathogens::RNASeq::GFF->new( filename => $self->annotation_filename);
}

sub _build__results_spreadsheet
{
  my ($self) = @_;
  Pathogens::RNASeq::ExpressionStatsSpreadsheet->new( output_filename => $self->sequence_filename."_expression.csv" );
}

sub _build__expression_results
{
  my ($self) = @_;
  my $total_mapped_reads = $self->_sequence_file->total_mapped_reads;
  my @expression_results = ();

  for my $feature_id (keys %{$self->_annotation_file->features})
  {
    my $alignment_slice = Pathogens::RNASeq::AlignmentSlice->new(
      filename           => $self->sequence_filename,
      total_mapped_reads => $total_mapped_reads,
      feature            => $self->_annotation_file->features->{$feature_id},
      filters            => $self->filters
      );
    my $alignment_slice_results = $alignment_slice->rpkm_values;
    $alignment_slice_results->{gene_id} = $feature_id;
    push(@expression_results, $alignment_slice_results);
  }
  
  return \@expression_results;
}


sub output_spreadsheet
{
  my ($self) = @_;
  for my $expression_result (@{$self->_expression_results})
  {
    $self->_results_spreadsheet->add_result($expression_result);
  }
  $self->_results_spreadsheet->build_and_close();
  return 1;
}


1;

