=head1 NAME

InsertSite.pm   - Take in a bam file and plot the start position of each read

=head1 SYNOPSIS

use Pathogens::RNASeq::InsertSite;
my $insertsite_plots_from_bam = Pathogens::RNASeq::InsertSite->new(
   filename => 'my_file.bam',
   output_base_filename => 'my_output_file'
  );
$insertsite_plots_from_bam->create_plots();


=cut
package Pathogens::RNASeq::InsertSite;
use Moose;
use VertRes::Parser::bam;
use Bio::DB::Sam;

has 'filename'                => ( is => 'rw', isa => 'Str',      required  => 1 );
has 'output_base_filename'    => ( is => 'rw', isa => 'Str',      required  => 1 );
                            
has '_input_file_handle'      => ( is => 'rw',                    lazy_build => 1 );
has '_output_file_handles'    => ( is => 'rw', isa => 'HashRef',  lazy_build => 1 );
has '_sequence_names'         => ( is => 'rw', isa => 'ArrayRef', lazy_build => 1 );
has '_sequence_base_counters' => ( is => 'rw', isa => 'HashRef',  lazy_build => 1 );
has '_sequence_information'   => ( is => 'rw', isa => 'HashRef',  lazy_build => 1 );

has '_frequency_of_read_start' => ( is => 'rw', isa => 'HashRef',  lazy => 1, builder => '_build__frequency_of_read_start');

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
    open($output_file_handles{$sequence_name}, '|-', " gzip >". $self->output_base_filename.".$sequence_name.insert_site_plot.gz") || Pathogens::RNASeq::Exceptions::FailedToCreateOutputFileHandle->throw(error => "Couldnt create output file handle for saving insertsite plot results for ". $sequence_name. " in ". $self->filename. " and output base ".$self->output_base_filename);
  }
  
  return \%output_file_handles;
}

sub _build__input_file_handle
{
  my ($self) = @_;
  return Bio::DB::Bam->open($self->filename);
}

sub _number_of_forward_reads
{
  my ($self, $sequence_name,$read_coord) = @_;
  return $self->_number_of_reads($sequence_name,$read_coord, 1);
}

sub _number_of_reverse_reads
{
  my ($self, $sequence_name,$read_coord) = @_;
  return $self->_number_of_reads($sequence_name,$read_coord, -1);
}

sub _number_of_reads
{
  my ($self, $sequence_name,$read_coord, $direction) = @_;
  if(defined($self->_frequency_of_read_start->{$sequence_name}{$read_coord}) && defined($self->_frequency_of_read_start->{$sequence_name}{$read_coord}{$direction}))
  {
    return $self->_frequency_of_read_start->{$sequence_name}{$read_coord}{$direction};
  }
  return 0;
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

sub _build__frequency_of_read_start
{
  my ($self) = @_;
  my %frequency_of_read_start ;
  my $header       = $self->_input_file_handle->header;
   my $target_names = $header->target_name;
  while (my $align = $self->_input_file_handle->read1) 
  {
    next if ($align->unmapped);
    my $seqid     = $target_names->[$align->tid];
    if($align->strand == 1)
    {
      $frequency_of_read_start{$seqid}{$align->start}{$align->strand}++;
    }
    else
    {
      $frequency_of_read_start{$seqid}{$align->end}{$align->strand}++;
    }
  }

  return \%frequency_of_read_start;
}

sub create_plots
{
  my ($self) = @_;
  my %read_starts = %{$self->_frequency_of_read_start};
  for my $sequence_name (keys %read_starts)
  {
    my %sequence_read_coords  = %{$read_starts{$sequence_name}};
    for my $read_coord (sort {$a <=> $b} (keys %sequence_read_coords))
    {
      my $padding_string = $self->_create_padding_string($self->_sequence_base_counters->{$sequence_name}, $read_coord);
      $self->_sequence_base_counters->{$sequence_name} = $read_coord;

      my $forward_reads = $self->_number_of_forward_reads($sequence_name,$read_coord);
      my $reverse_reads = $self->_number_of_reverse_reads($sequence_name,$read_coord);

      print { $self->_output_file_handles->{$sequence_name} } $padding_string.$forward_reads." ".$reverse_reads."\n";
    }
  }

  $self->_print_padding_at_end_of_sequence;
  $self->_close_output_file_handles;
  return 1;
}


1;
