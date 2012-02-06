=head1 NAME

Read.pm   - functionality for the strand specific protocol for reads

=head1 SYNOPSIS

use Pathogens::RNASeq::StrandSpecificProtocol::Read;
my $alignment_slice = Pathogens::RNASeq::StrandSpecificProtocol::Read->new(
  alignment_line => 'xxxxxxx',
  gene_strand => 1,
  exons => [[1,3],[4,5]]
  );
  my %mapped_reads = $alignment_slice->mapped_reads;
  $mapped_reads{sense};
  $mapped_reads{antisense};

=cut

package Pathogens::RNASeq::StrandSpecificProtocol::Read;
use Moose;
extends 'Pathogens::RNASeq::Read';

sub _process_read_details
{
  my ($self, $read_details) = @_;
  $read_details->{flag} = Pathogens::RNASeq::StrandSpecificProtocol::Read->_calculate_bitwise_flag($read_details->{flag});
}

sub _calculate_bitwise_flag
{
	my ($self, $flag) = @_;
	if ($flag & 128 && $flag & 2) 
  {
    if ($flag & 16) 
    {
      $flag = $flag - 16;     
    } 
    else 
    {
      $flag = $flag + 16;
    }
  }
  $flag = Pathogens::RNASeq::Read->_unmark_duplicates($flag);
  
  return $flag;
}


1;