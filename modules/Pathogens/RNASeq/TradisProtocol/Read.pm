=head1 NAME

Read.pm   - Tradis protocol, just inherits from the base read class

=head1 SYNOPSIS

use Pathogens::RNASeq::TradisProtocol::Read;
my $alignment_slice = Pathogens::RNASeq::TradisProtocol::Read->new(
  alignment_line => 'xxxxxxx',
  gene_strand => 1,
  exons => [[1,3],[4,5]]
  );
  my %mapped_reads = $alignment_slice->mapped_reads;
  $mapped_reads{sense};
  $mapped_reads{antisense};

=cut

package Pathogens::RNASeq::TradisProtocol::Read;
use Moose;
extends 'Pathogens::RNASeq::Read';

sub _process_read_details
{
  my ($self, $read_details) = @_;
  $read_details->{flag} = Pathogens::RNASeq::TradisProtocol::Read->_calculate_bitwise_flag($read_details->{flag});
}

sub _calculate_bitwise_flag
{
	my ($self, $flag) = @_;
  $flag = Pathogens::RNASeq::Read->_unmark_duplicates($flag);
  
  return $flag;
}


1;