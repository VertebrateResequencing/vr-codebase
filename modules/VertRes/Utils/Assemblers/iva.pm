=head1 NAME

VertRes::Utils::Assemblers::iva - assembly utility functions

=head1 SYNOPSIS

use VertRes::Utils::Assemblers::iva;

my $assembly_util = VertRes::Utils::Assemblers::iva->new();

# use any of the utility functions described here, eg.
my $assember = $assembler_class->new(
  assembler_exec => '/path/to/iva',
  optimiser_exec => '/path/to/iva',
);

my \$ok = \$assember->optimise_parameters(

  );

=head1 DESCRIPTION

iva-specific assembly functions

=head1 AUTHOR

path-help@sanger.ac.uk

=cut

package VertRes::Utils::Assemblers::iva;

use strict;
use warnings;
use VertRes::Wrapper::smalt;
use File::Copy;
use File::Spec;
use Utils;

use base qw(VertRes::Utils::Assembly);

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Utils::Assemblers::iva->new();
 Function: Create a new VertRes::Utils::Assemblers::iva object.
 Returns : VertRes::Utils::Assemblers::iva object
 Args    : n/a

=cut

sub new {
    my ($class, @args) = @_;

    my $self = $class->SUPER::new(@args);

    return $self;
}

sub optimise_parameters
{
  my ($self, $num_threads) = @_;
  my $trimming_opts = '';

  if (
      defined($self->{remove_adapters})
      and $self->{remove_adapters}
      and defined($self->{adapter_removal_tool})
      and defined($self->{adapters_file})
      and defined($self->{trimmomatic_jar})
  ) {
    $trimming_opts = "--trimmomatic $self->{trimmomatic_jar} --adapters $self->{adapters_file}";
  }

  if (
      defined($self->{remove_primers})
      and $self->{remove_primers}
      and defined($self->{primer_removal_tool})
      and defined($self->{primers_file})
  ) {
    $trimming_opts .= " --pcr_primers " . $self->{primers_file};
  }

  `$self->{optimiser_exec} $trimming_opts --fr $self->{files_str} --threads $num_threads iva_assembly`;
  File::Copy::move(File::Spec->catfile($self->optimised_directory(), 'contigs.fasta'), File::Spec->catfile($self->optimised_directory(), 'contigs.fa'));
  system("touch ". File::Spec->catfile($self->optimised_directory(), '_iva_optimise_parameters_done'));
  return 1;
}


sub optimised_directory
{
  my ($self) = @_;
  return File::Spec->catfile($self->{output_directory}, 'iva_assembly');
}

sub optimised_assembly_file_path
{
  my ($self) = @_;
  return File::Spec->catfile($self->optimised_directory(), 'contigs.fa');
}



=head2 generate_files_str

 Title   : generate_files_str
 Usage   : my $module = $obj->generate_files_str();
 Function: create the input string for the files to go into the assembler
 Returns : string which can be passed into assembler

=cut
sub generate_files_str
{
  my ($class, $pools_array, $directory) = @_;
  my $files_str = "";

  my $pool_count = 1;
  for my $lane_pool (@$pools_array)
  {
    $files_str .= "$directory/pool_$pool_count.fastq.gz";
  }
  return $files_str;
}


=head2 get_parameters

 Title   : get_parameters
 Usage   : my $params = $obj->get_parameters();
 Function: get parameters. Extend later to read the iva.log file
 Returns : hash containing parameters about final assembly

=cut
sub get_parameters
{
	my ($self, $filename) = @_;
	my %parameters;
	$parameters{assembly_directory} = File::Spec->catfile($self->{output_directory}, 'iva_assembly');
	return \%parameters;
}


=head2 estimate_memory_required

 Title   : estimate_memory_required
 Usage   : my $memory_required_in_kb = $obj->estimate_memory_required();
 Function: estimate the memory required for the assembler in KB
 Returns : integer in kb of memory requirement. IVA usually doesn't use more than ~1GB, so return 1.5GB

=cut
sub estimate_memory_required
{
  return 1500000;
}


sub name
{
  my ($self) = @_;
  return 'iva';
}

1;
