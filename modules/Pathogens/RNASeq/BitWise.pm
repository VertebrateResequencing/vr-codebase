=head1 NAME

BitWise.pm   - Manipulate the bitwise flags in a BAM file

=head1 SYNOPSIS

use Pathogens::RNASeq::BitWise;
my $bitwise = Pathogens::RNASeq::BitWise->new(
  filename => '/abc/my_file.bam',
  output_filename => 'output_file.bam'
  );
$bitwise->update_bitwise_flags;

=cut
package Pathogens::RNASeq::BitWise;
use Moose;
use Pathogens::RNASeq::Exceptions;
use VertRes::Wrapper::samtools;

has 'filename'                      => ( is => 'rw', isa => 'Str', required  => 1 );
has 'output_filename'               => ( is => 'rw', isa => 'Str', required  => 1 );
has 'protocol'                      => ( is => 'rw', isa => 'Str', default   => 'StandardProtocol' );

has '_input_sequence_data_filename' => ( is => 'rw', isa => 'Str'); # allow for testing
has '_sequence_data_file_handle'    => ( is => 'rw',               lazy_build => 1 );
has '_read_protocol_class'          => ( is => 'rw',               lazy_build => 1 );
has '_output_file_handle'           => ( is => 'rw',               lazy_build => 1 );


sub update_bitwise_flags
{
	my ($self) = @_;
	my $file_handle = $self->_sequence_data_file_handle;
  
  while(my $line = <$file_handle>)
  {
	  if($line =~ /^@/)
	  {
		  print {$self->_output_file_handle} $line;
		}
		else
		{
		  if($line =~ /^([^\t]+\t)([\d]+)(\t.+)$/ )
		  {
		    my $start_of_line = $1;
		    my $flag = $2;
		    my $end_of_line = $3;

  	    $flag = $self->_read_protocol_class->_calculate_bitwise_flag($flag);
  			print {$self->_output_file_handle} $start_of_line.$flag.$end_of_line."\n" ;
	    }
		  else
		  {
		    print {$self->_output_file_handle} $line;
	    }

	  }
	
	}
	close($self->_output_file_handle);
	$self->_index_output_file;
	return 1;
}

sub _index_output_file
{
  my ($self) = @_;
  my $samtools = VertRes::Wrapper::samtools->new();
  $samtools->index($self->output_filename, $self->output_filename.".bai");
  return 1;
}


sub _build__output_file_handle
{
  my ($self) = @_;
  my $output_file_handle; 
  open($output_file_handle, '|- '," samtools view -bS - > ". $self->output_filename) || Pathogens::RNASeq::Exceptions::FailedToCreateNewBAM->throw(error => "Couldnt open output file for writing ". $self->output_filename);
  return $output_file_handle;
}

sub _build__read_protocol_class
{
	my ($self) = @_;
	my $read_protocol_class = "Pathogens::RNASeq::".$self->protocol."::Read";
	eval("use $read_protocol_class");
  return $read_protocol_class;
}

sub _build__sequence_data_file_handle
{
  my ($self) = @_;
  my $sequence_data_file_handle;
  open($sequence_data_file_handle, $self->_sequence_data_stream ) || Pathogens::RNASeq::Exceptions::FailedToOpenAlignmentSlice->throw( error => "Cant view ".$self->filename."" );
  return $sequence_data_file_handle;
}

sub _sequence_data_stream
{
  my ($self) = @_;
  if($self->_input_sequence_data_filename)
  {
    return $self->_input_sequence_data_filename;
  }
  else
  {
    return "samtools view -h ".$self->filename." |";
  }
}

1;

