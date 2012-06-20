=head1 NAME

CoveragePlot.pm   - Take in a sequencing file, do a pileup and generate a dot plot of the number of reads at each position.Assumes the BAM is sorted.

=head1 SYNOPSIS

use Pathogens::RNASeq::CoveragePlot;
my $coverage_plots_from_bam = Pathogens::RNASeq::CoveragePlot->new(
   filename => 'my_file.bam',
   output_base_filename => 'my_output_file',
   mpileup_cmd   => "samtools mpileup ",
   mapping_quality => 30
  );
$coverage_plots_from_bam->create_plots();


=cut
package Pathogens::RNASeq::CoveragePlot;
use Moose;
use VertRes::Parser::bam;

has 'filename'                => ( is => 'rw', isa => 'Str',      required  => 1 );
has 'output_base_filename'    => ( is => 'rw', isa => 'Str',      required  => 1 );
# optional inputs
has 'mpileup_cmd'             => ( is => 'rw', isa => 'Str',      default   => "samtools mpileup -d 1000" );
has 'mapping_quality'         => ( is => 'rw', isa => 'Int',      default   => 0 );
                            
has '_input_file_handle'      => ( is => 'rw',                    lazy_build => 1 );
has '_output_file_handles'    => ( is => 'rw', isa => 'HashRef',  lazy_build => 1 );
has '_sequence_names'         => ( is => 'rw', isa => 'ArrayRef', lazy_build => 1 );
has '_sequence_base_counters' => ( is => 'rw', isa => 'HashRef',  lazy_build => 1 );
has '_sequence_information'   => ( is => 'rw', isa => 'HashRef',  lazy_build => 1 );

sub _build__sequence_information
{
  my ($self) = @_;
  my %all_sequences_info = VertRes::Parser::bam->new( file => $self->filename )->sequence_info();
  return \%all_sequences_info;
}

sub _build__sequence_names
{
   my ($self) = @_;
   my @sequence_names = keys %{$self->_sequence_information} ;
   return \@sequence_names;
}

sub _build__sequence_base_counters
{
  my ($self) = @_;
  my %sequence_base_counters;
  for my $sequence_name (@{$self->_sequence_names} )
  {
    $sequence_base_counters{$sequence_name} = 0;
  }
  return \%sequence_base_counters;
}

sub _build__output_file_handles
{
  my ($self) = @_;
  my %output_file_handles;
  for my $sequence_name (@{$self->_sequence_names} )
  {
    open($output_file_handles{$sequence_name}, '|-', " gzip >". $self->output_base_filename.".$sequence_name.coverageplot.gz") || Pathogens::RNASeq::Exceptions::FailedToCreateOutputFileHandle->throw(error => "Couldnt create output file handle for saving coverage plot results for ". $sequence_name. " in ". $self->filename. " and output base ".$self->output_base_filename);
  }
  
  return \%output_file_handles;
}

sub _build__input_file_handle
{
  my ($self) = @_;
  my $input_file_handle; 
  # TODO remove direct call to samtools and allow for filtering options
  # this only extracts the sequence name, base position, bases.
  open($input_file_handle, '-|', $self->mpileup_cmd." -q ".$self->mapping_quality." " $self->filename.' | awk \'{print $1"\t"$2"\t"$5}\'') || Pathogens::RNASeq::Exceptions::FailedToCreateMpileup->throw(error => "Failed to create mpileup for ".$self->filename );
  return $input_file_handle;
}

sub _number_of_forward_reads
{
  my ($self, $read_string) = @_;
  return $self->_number_of_reads("[ACGTN]", $read_string);  
}

sub _number_of_reverse_reads
{
  my ($self, $read_string) = @_;
  return $self->_number_of_reads("[acgtn]", $read_string);
}

sub _number_of_reads
{
  my ($self,$base_regex, $read_string) = @_;
  $_ = $read_string;
  my $base_count = s/$base_regex//g;
  $base_count = ($base_count eq "") ? 0 : $base_count;
  
  while ($read_string =~ /[\+-]([\d]+)$base_regex/g) 
  {
    $base_count -= $1;
  }
  $base_count = ($base_count < 0) ? 0 : $base_count;
  return $base_count;
}

# work out if padding is needed and return it as a formatted string
sub _create_padding_string
{
  my ($self,$previous_counter, $current_counter) = @_;
  my $padding_string = "";
  for(my $i = $previous_counter+1 ; $i < $current_counter; $i++)
  {
    $padding_string .= "0 0\n";
  }
  return $padding_string;
}

sub _print_padding_at_end_of_sequence
{
   my ($self) = @_;
   for my $sequence_name (@{$self->_sequence_names} )
   {
     my $sequence_length = $self->_sequence_information->{$sequence_name}->{'LN'};
     next unless($sequence_length =~ /^[\d]+$/);
     $sequence_length++;
     my $padding_string = $self->_create_padding_string($self->_sequence_base_counters->{$sequence_name}, $sequence_length);
     $self->_sequence_base_counters->{$sequence_name} = $sequence_length;
     print { $self->_output_file_handles->{$sequence_name} } $padding_string;
   }
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

sub create_plots
{
  my ($self) = @_;
  my $input_file_handle = $self->_input_file_handle;

  while(my $line = <$input_file_handle>)
  {
    my($sequence_name, $base_position, $read_string) = split(/\t/, $line);
    my $padding_string = $self->_create_padding_string($self->_sequence_base_counters->{$sequence_name},$base_position);

    $self->_sequence_base_counters->{$sequence_name} = $base_position;
    my $forward_reads = $self->_number_of_forward_reads($read_string);
    my $reverse_reads = $self->_number_of_reverse_reads($read_string);
    
    print { $self->_output_file_handles->{$sequence_name} } $padding_string.$forward_reads." ".$reverse_reads."\n";
  }
  $self->_print_padding_at_end_of_sequence;
  $self->_close_output_file_handles;
  return 1;
}



1;
