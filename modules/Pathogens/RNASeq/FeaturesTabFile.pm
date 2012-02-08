=head1 NAME

FeaturesTabFile.pm   - Builds a spreadsheet of expression results

=head1 SYNOPSIS

use Pathogens::RNASeq::FeaturesTabFile;
my $tab_file_results = Pathogens::RNASeq::FeaturesTabFile->new(
  output_filename => '/abc/my_results.tab',
  features => $features;
);
$tab_file_results->create_file;

=cut
package Pathogens::RNASeq::FeaturesTabFile;
use Moose;


has 'output_filename'      => ( is => 'rw', isa => 'Str',      required => 1 );
has 'features'             => ( is => 'rw', isa => 'HashRef',  required => 1 );
has 'sequence_names'       => ( is => 'rw', isa => 'ArrayRef', required => 1 );

has '_output_file_handles' => ( is => 'rw', lazy_build => 1 );

sub _build__output_file_handles
{
  my ($self) = @_;
  my %output_file_handles;
	for my $sequence_name (@{$self->sequence_names} )
  {
	  open($output_file_handles{$sequence_name}, '|-', " gzip >". $self->output_filename.".$sequence_name.gz")  || Pathogens::RNASeq::Exceptions::FailedToOpenFeaturesTabFileForWriting->throw( error => "Cant open ".$self->output_filename." for writing");
	}
  return \%output_file_handles;
}

sub _print_feature
{
  my ($self,$feature) = @_;
  print {$self->_output_file_handles->{$feature->seq_id} } "FT   misc_feature    ".$feature->gene_start."..".$feature->gene_end."\n";
  print {$self->_output_file_handles->{$feature->seq_id} } "FT                   /colour=2\n";
  print {$self->_output_file_handles->{$feature->seq_id} } "FT                   /gene=\"".$feature->gene_id."\"\n";
  print {$self->_output_file_handles->{$feature->seq_id} } "FT                   /seq_id=\"".$feature->seq_id."\"\n";
}

sub create_files
{
  my ($self) = @_;
  for my $feature (sort keys %{$self->features})
  {
    $self->_print_feature($self->features->{$feature});
  }
  
 $self->_close_output_file_handles;
 return 1;
}

sub _close_output_file_handles
{
  my ($self) = @_;
  for my $output_file_handle (values %{$self->_output_file_handles} )
  {
    close($output_file_handle);
  }
  return;
}


1;
