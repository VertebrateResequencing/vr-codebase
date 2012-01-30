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

  if ($read_details->{flag} & 128 && $read_details->{flag} & 2) 
  {
    if ($read_details->{flag} & 16) 
    {
      $read_details->{flag} = $read_details->{flag} - 16;     
    } 
    else 
    {
      $read_details->{flag} = $read_details->{flag} + 16;
    }
  } 

}


1;