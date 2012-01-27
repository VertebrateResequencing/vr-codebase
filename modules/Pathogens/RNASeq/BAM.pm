=head1 NAME

BAM.pm   - Represents a BAM file for use in RNASeq

=head1 SYNOPSIS

use Pathogens::RNASeq::BAM;
my $rna_seq_bam = Pathogens::RNASeq::BAM->new(
  filename => '/abc/my_file.bam'
  );
  $rna_seq_bam->total_mapped_reads;

=cut
package Pathogens::RNASeq::BAM;
use Moose;
extends 'Pathogens::RNASeq::BAMStats';

has 'filename'           => ( is => 'rw', isa => 'Str', required   => 1 );
1;