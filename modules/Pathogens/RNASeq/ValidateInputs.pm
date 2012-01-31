=head1 NAME

ValidateInputs.pm   - Validate the input sequence file and the annotation file

=head1 SYNOPSIS

use Pathogens::RNASeq::ValidateInputs;
my $bam_container = Pathogens::RNASeq::ValidateInputs->new(
  sequence_filename => 'my_aligned_sequence.bam',
  annotation_filename => 'my_annotation_file.gff'
  );

=cut
package Pathogens::RNASeq::ValidateInputs;
use Moose;

has 'sequence_filename'     => ( is => 'rw', isa => 'Str', required => 1 );
has 'annotation_filename'   => ( is => 'rw', isa => 'Str', required => 1 );

# are the input files valid? return 1 if true; 0 if false
sub are_input_files_valid
{
	return 1;
}


1;