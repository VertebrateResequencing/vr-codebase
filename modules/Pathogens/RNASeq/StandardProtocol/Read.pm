=head1 NAME

Read.pm   - Standard protocol, just inherits from the base read class

=head1 SYNOPSIS

use Pathogens::RNASeq::StandardProtocol::Read;
my $alignment_slice = Pathogens::RNASeq::StandardProtocol::Read->new(
  alignment_line => 'xxxxxxx',
  gene_strand => 1,
  exons => [[1,3],[4,5]]
  );
  my %mapped_reads = $alignment_slice->mapped_reads;
  $mapped_reads{sense};
  $mapped_reads{antisense};

=cut

package Pathogens::RNASeq::StandardProtocol::Read;
use Moose;
extends 'Pathogens::RNASeq::Read';

1;