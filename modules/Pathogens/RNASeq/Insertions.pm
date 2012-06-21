=head1 NAME

Insertions.pm   - Find the expression when given an input aligned file and an annotation file

=head1 SYNOPSIS

use Pathogens::RNASeq::Insertions;
my $expression_results = Pathogens::RNASeq::Insertion->new(
  sequence_filename => 'my_aligned_sequence.bam',
  annotation_filename => 'my_annotation_file.gff',
  output_base_filename => 'my_alignement_basename'
  );

$expression_results->output_spreadsheet();

=cut
package Pathogens::RNASeq::Insertions;
use Moose;
use Pathogens::RNASeq::SequenceFile;
use Pathogens::RNASeq::GFF;
use Pathogens::RNASeq::AlignmentSlice;
use Pathogens::RNASeq::InsertionStatsSpreadsheet;
use Pathogens::RNASeq::ValidateInputs;
use Pathogens::RNASeq::Exceptions;
use Pathogens::RNASeq::BitWise;
use Pathogens::RNASeq::IntergenicRegions;
use Pathogens::RNASeq::FeaturesTabFile;
use Pathogens::RNASeq::InsertSite;

has 'sequence_filename'       => ( is => 'rw', isa => 'Str', required => 1 );
has 'annotation_filename'     => ( is => 'rw', isa => 'Str', required => 1 );
has 'output_base_filename'    => ( is => 'rw', isa => 'Str', required => 1 );

#optional input parameters
has 'filters'                 => ( is => 'rw', isa => 'Maybe[HashRef]'     );
has 'protocol'                => ( is => 'rw', isa => 'Str',  default => 'TradisProtocol' );
has 'samtools_exec'           => ( is => 'rw', isa => 'Str',  default => "samtools" );
has 'intergenic_regions'      => ( is => 'rw', isa => 'Bool', default => 1 );
has 'minimum_intergenic_size' => ( is => 'rw', isa => 'Int',  default => 10 );

has '_sequence_file'          => ( is => 'rw', isa => 'Pathogens::RNASeq::SequenceFile',               lazy_build  => 1 );
has '_annotation_file'        => ( is => 'rw', isa => 'Pathogens::RNASeq::GFF',                        lazy_build  => 1 );
has '_results_spreadsheet'    => ( is => 'rw', isa => 'Pathogens::RNASeq::InsertionStatsSpreadsheet',  lazy_build  => 1 );
has '_insertion_results'      => ( is => 'rw', isa => 'ArrayRef',                                      lazy_build  => 1 );
has '_frequency_of_read_start' => ( is => 'rw', isa => 'HashRef', lazy  => 1, builder => '_build__frequency_of_read_start' );

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
  Pathogens::RNASeq::InsertionStatsSpreadsheet->new( output_filename => $self->output_base_filename.".insertion.csv", protocol => $self->protocol);
}

sub _corrected_sequence_filename
{
  my ($self) = @_;
  return $self->output_base_filename.".corrected.bam";
}

sub _build__frequency_of_read_start
{
  my ($self) = @_;
  my $insertsite_plots_from_bam = Pathogens::RNASeq::InsertSite->new(
     filename => $self->_corrected_sequence_filename,
     output_base_filename => $self->_corrected_sequence_filename
    );
  return $insertsite_plots_from_bam->_frequency_of_read_start();
}

sub _count_insertions_in_feature
{
   my ($self, $feature) = @_;
   my $seqid = $feature->seq_id;
   my %feature_insertion_details ;
   $feature_insertion_details{pos_insert_sites} = 0;
   $feature_insertion_details{neg_insert_sites} = 0;
   
   my $i=0;
   my $insert_sites = 0;
   my $insert_site_reads = 0;
   for($i = $feature->gene_start; $i<= $feature->gene_end ; $i++ )
   {
      if(defined($self->_frequency_of_read_start->{$seqid}) && defined($self->_frequency_of_read_start->{$seqid}{$i}) )
      {
        if(defined($self->_frequency_of_read_start->{$seqid}{$i}{"1"}) && $self->_frequency_of_read_start->{$seqid}{$i}{"1"} > 0)
        {
          $insert_site++;
          $insert_site_reads += $self->_frequency_of_read_start->{$seqid}{$i}{"1"} ;
        }

       if(defined($self->_frequency_of_read_start->{$seqid}{$i}{"-1"}) && $self->_frequency_of_read_start->{$seqid}{$i}{"-1"} > 0)
       {
         $insert_site++;
         $insert_site_reads += $self->_frequency_of_read_start->{$seqid}{$i}{"-1"} ;
       }
     }
   }
   my $length  = (($feature->gene_end +1) - $feature->gene_start);
   
   if($feature->strand == -1)
   {
     $feature_insertion_details{neg_insert_sites} = $insert_sites;
     $feature_insertion_details{neg_insert_site_reads} = $insert_site_reads;
   }
   elsif($feature->strand == 1)
   {
     $feature_insertion_details{pos_insert_sites} = $insert_sites;
     $feature_insertion_details{pos_insert_site_reads} = $insert_site_reads;
   }
   else
   {
     $feature_insertion_details{zero_insert_sites} = $insert_sites;
     $feature_insertion_details{zero_insert_site_reads} = $insert_site_reads;
   }
   
   $feature_insertion_details{normalised_pos_insert_sites} = $feature_insertion_details{pos_insert_sites} /$length;
   $feature_insertion_details{normalised_neg_insert_sites} = $feature_insertion_details{neg_insert_sites} /$length;
   $feature_insertion_details{normalised_zero_insert_sites} = $feature_insertion_details{zero_insert_sites} /$length;
   $feature_insertion_details{total_insert_sites}              = ($feature_insertion_details{neg_insert_sites} + $feature_insertion_details{neg_insert_sites} + $feature_insertion_details{zero_insert_sites});
   $feature_insertion_details{normalised_total_insert_sites}   = $feature_insertion_details{total_insert_sites} /$length;
   
   
   $feature_insertion_details{normalised_pos_insert_site_reads} = $feature_insertion_details{pos_insert_site_reads} /$length;
   $feature_insertion_details{normalised_neg_insert_site_reads} = $feature_insertion_details{neg_insert_site_reads} /$length;
   $feature_insertion_details{normalised_zero_insert_site_reads} = $feature_insertion_details{zero_insert_site_reads} /$length;
   $feature_insertion_details{total_insert_site_reads} =  $feature_insertion_details{pos_insert_site_reads}  + $feature_insertion_details{neg_insert_site_reads}  +  $feature_insertion_details{zero_insert_site_reads};
   $feature_insertion_details{normalised_total_insert_site_reads} = $feature_insertion_details{total_insert_site_reads} /$length;
   
   return \%feature_insertion_details;
}

sub _build__insertion_results
{
  my ($self) = @_;

  Pathogens::RNASeq::BitWise->new(
      filename        => $self->sequence_filename,
      output_filename => $self->_corrected_sequence_filename,
      protocol        => $self->protocol,
      samtools_exec   => $self->samtools_exec
    )->update_bitwise_flags();
  
  my @all_insertion_results = ();
  
  
  for my $feature_id (keys %{$self->_annotation_file->features})
  {
    my $insertion_results_feature  = $self->_count_insertions_in_feature($self->_annotation_file->features->{$feature_id});
    
    $insertion_results_feature->{gene_id} = $feature_id;
    $insertion_results_feature->{seq_id}  =  $self->_annotation_file->features->{$feature_id}->seq_id;
    push(@all_insertion_results, $insertion_results_feature);
  }
  
  if(defined($self->intergenic_regions) && $self->intergenic_regions == 1)
  {
    $self->_calculate_values_for_intergenic_regions(\@all_insertion_results);
  }
  
  return \@all_insertion_results;
}

sub _calculate_values_for_intergenic_regions
{
   my ($self, $all_insertion_results, $total_mapped_reads) = @_;
   # get intergenic regions
   my $intergenic_regions = Pathogens::RNASeq::IntergenicRegions->new(
     features       => $self->_annotation_file->features,
     window_margin  => 0,
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
     
     my $insertion_results_feature  = $self->_count_insertions_in_feature($intergenic_regions->intergenic_features->{$feature_id});

     $insertion_results_feature->{gene_id} = $feature_id;
     $insertion_results_feature->{seq_id}  =  $intergenic_regions->intergenic_features->{$feature_id}->seq_id;
     push(@{$all_insertion_results}, $insertion_results_feature);
   }
   
   return $all_insertion_results;
}


sub output_spreadsheet
{
  my ($self) = @_;
  for my $expression_result (@{$self->_insertion_results})
  {
    $self->_results_spreadsheet->add_result($expression_result);
  }
  $self->_results_spreadsheet->build_and_close();
  return 1;
}


1;

