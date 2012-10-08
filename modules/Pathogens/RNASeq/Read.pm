=head1 NAME

Read.pm   - Extract a slice of reads for a sequence file within a specific region

=head1 SYNOPSIS

use Pathogens::RNASeq::Read;
my $alignment_slice = Pathogens::RNASeq::Read->new(
  alignment_line => 'xxxxxxx',
  gene_strand => 1,
  exons => [[1,3],[4,5]]
  );
  my %mapped_reads = $alignment_slice->mapped_reads;
  $mapped_reads{sense};
  $mapped_reads{antisense};

=cut
package Pathogens::RNASeq::Read;
use Moose;
use Pathogens::RNASeq::Exceptions;

has 'alignment_line' => ( is => 'rw', isa => 'Str',      required   => 1 );
has 'exons'          => ( is => 'rw', isa => 'ArrayRef', required   => 1 );
has 'gene_strand'    => ( is => 'rw', isa => 'Int',      required   => 1 );

#optional filters
has 'filters'        => ( is => 'rw', isa => 'Maybe[HashRef]'            );

has '_read_details'  => ( is => 'rw', isa => 'HashRef',  lazy_build   => 1 );
has '_read_length'   => ( is => 'rw', isa => 'Int',      lazy_build   => 1 );
has '_read_position' => ( is => 'rw', isa => 'Int',      lazy_build   => 1 );

has 'read_strand'   => ( is => 'rw', isa => 'Int',      lazy_build   => 1 );
has 'mapped_reads' => ( is => 'rw', isa => 'HashRef',  lazy_build   => 1 );

sub _build__read_details
{
  my ($self) = @_;
  
  my($qname, $flag, $rname, $read_position, $mapq, $cigar, $mrnm, $mpos, $isize, $seq, $qual) = split(/\t/,$self->alignment_line);
  
  my $read_details = {
    mapping_quality => $mapq,
    cigar           => $cigar,
    read_position   => $read_position,
    flag            => $flag,
  };
  # hook to allow for reads to be altered for different protocols
  $self->_process_read_details($read_details);
  
  return $read_details;
}

sub _build__read_length
{
  my ($self) = @_;
  my $read_length = 0;
  $self->_read_details->{cigar} =~ s/(\d+)[MIS=X]/$read_length+=$1/eg;
  return $read_length;
}

sub _build_read_strand
{
  my ($self) = @_;
  my $read_strand = $self->_read_details->{flag} & 16 ? -1:1; # $flag bit set to 0 for forward 1 for reverse.
  return $read_strand;
}

sub _build__read_position
{
  my ($self) = @_;
  $self->_read_details->{read_position};
}


sub _build_mapped_reads
{
  my ($self) = @_;
  my %mapped_reads ;
  $mapped_reads{sense} = 0;
  $mapped_reads{antisense} = 0;
  
  return \%mapped_reads unless( $self->_does_read_pass_filters() == 1 );
  
  foreach my $exon (@{$self->exons})
  {
    my($exon_start,$exon_end) = @{$exon};
    if($self->_read_position <= $exon_end && ($self->_read_position + $self->_read_length - 1) >= $exon_start)
    {
      
      if($self->read_strand == $self->gene_strand) 
      {
        $mapped_reads{sense}++;
      }
      else
      {
        $mapped_reads{antisense}++;
      }
      last;
    }
  }
  
  return \%mapped_reads;
}

sub _process_read_details
{
  # hook to allow for reads to be altered for different protocols
}

sub _does_read_pass_filters
{
  my ($self) = @_;
  
  # filter unmapped read
  if( ($self->_read_details->{flag} & 4) == 4)
  {
  	return 0;
  }
  return 1 unless(defined($self->filters));
  
  if(defined($self->filters->{mapping_quality}) && ($self->_read_details->{mapping_quality}  <= $self->filters->{mapping_quality}) )
  {
     return 0;
  }
  
  if(defined($self->filters->{bitwise_flag}) && ( ($self->_read_details->{flag} & $self->filters->{bitwise_flag})  == 0 ))
  {
    return 0;
  }
  
  return 1;
}

# class method
sub _unmark_duplicates
{
  my ($self, $flag) = @_;
  if ($flag & 1024) 
  {
    $flag = $flag - 1024;     
  }
  return $flag;
}

sub _calculate_bitwise_flag
{
	# hook for protocols to update the bitwise flag
	my ($self, $flag) = @_;
	return $flag;
}

1;