=head1 NAME

ValidateInputs.pm   - Validate the input sequence file and the annotation file

=head1 SYNOPSIS

use Pathogens::RNASeq::ValidateInputs;
my $validator = Pathogens::RNASeq::ValidateInputs->new(
  sequence_filename => 'my_aligned_sequence.bam',
  annotation_filename => 'my_annotation_file.gff'
  );
$validator->are_input_files_valid();

=cut
package Pathogens::RNASeq::ValidateInputs;
use Moose;
use VertRes::Parser::bam;
use Bio::Tools::GFF;

has 'sequence_filename'           => ( is => 'rw', isa => 'Str',             required   => 1 );
has 'annotation_filename'         => ( is => 'rw', isa => 'Str',             required   => 1 );


has '_actual_sequence_details'    => ( is => 'rw', isa => 'HashRef',         lazy_build => 1 );
has '_annotated_sequence_details' => ( is => 'rw', isa => 'HashRef',         lazy_build => 1 );
has '_gff_parser'                 => ( is => 'rw', isa => 'Bio::Tools::GFF', lazy_build => 1 );

# are the input files valid? return 1 if true; 0 if false
sub are_input_files_valid
{
	my ($self) = @_;
	return 0 unless $self->_sequence_names_match();
	return 0 unless $self->_lengths_match();
	
	return 1;
}

sub _sequence_names_match
{
	my ($self) = @_;
	
	for my $sequence_name (keys %{$self->_actual_sequence_details})
	{
		# cant find the sequence name
		return 0 unless(defined ($self->_annotated_sequence_details->{$sequence_name}));
	}
	
	return 1;
}

sub _lengths_match
{
	my ($self) = @_;
	
	for my $sequence_name (keys %{$self->_actual_sequence_details})
	{
		# cant find the sequence name
		return 0 unless(defined ($self->_annotated_sequence_details->{$sequence_name}));
		if($self->_annotated_sequence_details->{$sequence_name } != $self->_actual_sequence_details->{$sequence_name} )
		{
			return 0;
		}
	}
	
	return 1;
}

sub _build__actual_sequence_details
{
	my ($self) = @_;
	my %actual_sequence_details;
	my %all_sequences_info = VertRes::Parser::bam->new(file => $self->sequence_filename)->sequence_info();
  for my $sequence_name (keys %all_sequences_info)
  {
	  $actual_sequence_details{$sequence_name} = $all_sequences_info{$sequence_name}->{LN};
	}

  return \%actual_sequence_details;
}


# create a hash with sequence names and the length of the sequence
sub _build__annotated_sequence_details
{
	my ($self) = @_;
	my %annotated_sequence_details;
	while(my $sequence_region = $self->_gff_parser->next_segment())
	{
		$annotated_sequence_details{$sequence_region->id()} = $sequence_region->end() + 1 - $sequence_region->start();
	}
	return \%annotated_sequence_details;
}

sub _build__gff_parser
{
  my ($self) = @_;
  Bio::Tools::GFF->new(-gff_version => 3, -file => $self->annotation_filename);
}

1;