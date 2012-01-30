=head1 NAME

SequenceFile.pm   - Represents a SequenceFile file for use in RNASeq

=head1 SYNOPSIS

use Pathogens::RNASeq::SequenceFile;
my $rna_seq_bam = Pathogens::RNASeq::SequenceFile->new(
  filename => '/abc/my_file.bam'
  );
  $rna_seq_bam->total_mapped_reads;

=cut
package Pathogens::RNASeq::SequenceFile;
use Moose;
extends 'Pathogens::RNASeq::BAMStats';

has 'filename'           => ( is => 'rw', isa => 'Str', required   => 1 );
1;