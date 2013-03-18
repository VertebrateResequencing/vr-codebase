=head1 NAME

Expression.pm   - Find the expression when given an input aligned file and an annotation file

=head1 SYNOPSIS

use Pathogens::RNASeq::Expression;
my $expression_results = Pathogens::RNASeq::Expression->new(
  sequence_filename => 'my_aligned_sequence.bam',
  annotation_filename => 'my_annotation_file.gff',
  output_base_filename => 'my_alignement_basename'
  );

$expression_results->output_spreadsheet();

=cut
package Pathogens::RNASeq::Expression;
use Moose;
use Pathogens::RNASeq::SequenceFile;
use Pathogens::RNASeq::GFF;
use Pathogens::RNASeq::AlignmentSlice;
use Pathogens::RNASeq::ExpressionStatsSpreadsheet;
use Pathogens::RNASeq::ValidateInputs;
use Pathogens::RNASeq::Exceptions;
use Pathogens::RNASeq::BitWise;
use Pathogens::RNASeq::IntergenicRegions;
use Pathogens::RNASeq::FeaturesTabFile;

has 'sequence_filename'       => ( is => 'rw', isa => 'Str', required => 1 );
has 'annotation_filename'     => ( is => 'rw', isa => 'Str', required => 1 );
has 'output_base_filename'    => ( is => 'rw', isa => 'Str', required => 1 );

#optional input parameters
has 'filters'                 => ( is => 'rw', isa => 'Maybe[HashRef]'     );
has 'protocol'                => ( is => 'rw', isa => 'Str',  default => 'StandardProtocol' );
has 'samtools_exec'           => ( is => 'rw', isa => 'Str',  default => "samtools" );
has 'window_margin'           => ( is => 'rw', isa => 'Int',  default => 50 );
has 'intergenic_regions'      => ( is => 'rw', isa => 'Bool', default => 0 );
has 'minimum_intergenic_size' => ( is => 'rw', isa => 'Int',  default => 10 );

has '_sequence_file'          => ( is => 'rw', isa => 'Pathogens::RNASeq::SequenceFile',               lazy_build  => 1 );
has '_annotation_file'        => ( is => 'rw', isa => 'Pathogens::RNASeq::GFF',                        lazy_build  => 1 );
has '_results_spreadsheet'    => ( is => 'rw', isa => 'Pathogens::RNASeq::ExpressionStatsSpreadsheet', lazy_build  => 1 );
has '_expression_results'     => ( is => 'rw', isa => 'ArrayRef',                                      lazy_build  => 1 );
has '_alignment_slice_protocol_class'  => ( is => 'rw',                                              lazy_build  => 1 );

sub _build__sequence_file
{
  my ($self) = @_;

	my $validator = Pathogens::RNASeq::ValidateInputs->new( sequence_filename => $self->sequence_filename, annotation_filename => $self->annotation_filename);
	if($validator->are_input_files_valid() == 0)
	{
		Pathogens::RNASeq::Exceptions::FailedToOpenAlignmentSlice->throw( error => "Input files invalid: ".$self->sequence_filename." ".$self->annotation_filename."\n" );
	}
	
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
  Pathogens::RNASeq::ExpressionStatsSpreadsheet->new( output_filename => $self->output_base_filename.".expression.csv", protocol => $self->protocol);
}

sub _corrected_sequence_filename
{
  my ($self) = @_;
  return $self->output_base_filename.".corrected.bam";
}

sub _build__expression_results
{
  my ($self) = @_;
  my $total_mapped_reads = $self->_sequence_file->total_mapped_reads;
  
  Pathogens::RNASeq::BitWise->new(
      filename        => $self->sequence_filename,
      output_filename => $self->_corrected_sequence_filename,
      protocol        => $self->protocol,
      samtools_exec   => $self->samtools_exec
    )->update_bitwise_flags();
  
  my @expression_results = ();
  
  for my $feature_id (keys %{$self->_annotation_file->features})
  {
    my $alignment_slice = $self->_alignment_slice_protocol_class->new(
      filename           => $self->_corrected_sequence_filename,
      total_mapped_reads => $total_mapped_reads,
      feature            => $self->_annotation_file->features->{$feature_id},
      filters            => $self->filters,
      protocol           => $self->protocol,
      samtools_exec      => $self->samtools_exec,
      window_margin      => $self->window_margin
      );
    my $alignment_slice_results = $alignment_slice->rpkm_values;
    
    $alignment_slice_results->{gene_id} = $feature_id;
    $alignment_slice_results->{seq_id}  =  $self->_annotation_file->features->{$feature_id}->seq_id;
    $alignment_slice_results->{locus_tag}  =  $self->_annotation_file->features->{$feature_id}->locus_tag;
    push(@expression_results, $alignment_slice_results);
  }
  
  if(defined($self->intergenic_regions) && $self->intergenic_regions == 1)
  {
    $self->_calculate_values_for_intergenic_regions(\@expression_results,$total_mapped_reads );
  }
  
  return \@expression_results;
}

sub _calculate_values_for_intergenic_regions
{
   my ($self, $expression_results, $total_mapped_reads) = @_;
   # get intergenic regions
   my $intergenic_regions = Pathogens::RNASeq::IntergenicRegions->new(
     features       => $self->_annotation_file->features,
     window_margin  => $self->window_margin,
     minimum_size   => $self->minimum_intergenic_size,
     sequence_lengths => $self->_annotation_file->sequence_lengths
     );

  # print out the features into a tab file for loading into Artemis
     my $tab_file_results = Pathogens::RNASeq::FeaturesTabFile->new(
       output_filename => $self->_corrected_sequence_filename.".intergenic",
       features        => $intergenic_regions->intergenic_features,
       sequence_names  => $intergenic_regions->sequence_names
     );
     $tab_file_results->create_files;

   for my $feature_id (keys %{$intergenic_regions->intergenic_features})
   {
     my $alignment_slice = $self->_alignment_slice_protocol_class->new(
       filename           => $self->_corrected_sequence_filename,
       total_mapped_reads => $total_mapped_reads,
       feature            => $intergenic_regions->intergenic_features->{$feature_id},
       filters            => $self->filters,
       protocol           => $self->protocol,
       samtools_exec      => $self->samtools_exec,
       window_margin      => 0,
       );
     my $alignment_slice_results = $alignment_slice->rpkm_values;

     $alignment_slice_results->{gene_id} = $feature_id;
     $alignment_slice_results->{seq_id}  = $intergenic_regions->intergenic_features->{$feature_id}->seq_id;
     push(@{$expression_results}, $alignment_slice_results);
   }
   
   return $expression_results;
}


sub _build__alignment_slice_protocol_class
{
	my ($self) = @_;
	my $alignment_slice_protocol_class = "Pathogens::RNASeq::".$self->protocol."::AlignmentSlice";
	eval("use $alignment_slice_protocol_class");
  return $alignment_slice_protocol_class;
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

