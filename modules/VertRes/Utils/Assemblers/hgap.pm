=head1 NAME

VertRes::Utils::Assemblers::hgap - assembly utility functions

=head1 SYNOPSIS

use VertRes::Utils::Assemblers::hgap;

my $assembly_util = VertRes::Utils::Assemblers::hgap->new();

# use any of the utility functions described here, eg.
my $assember = $assembler_class->new(
  assembler_exec => '/path/to/hgap',
  optimiser_exec => '/path/to/hgap',
);

my \$ok = \$assember->optimise_parameters(

  );

=head1 DESCRIPTION

hgap-specific assembly functions

=head1 AUTHOR

path-help@sanger.ac.uk

=cut

package VertRes::Utils::Assemblers::hgap;

use strict;
use warnings;
use VertRes::Wrapper::smalt;
use File::Copy;
use File::Spec;
use Utils;

use base qw(VertRes::Utils::Assembly);

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Utils::Assemblers::hgap->new();
 Function: Create a new VertRes::Utils::Assemblers::hgap object.
 Returns : VertRes::Utils::Assemblers::hgap object
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
  return 1;
}


sub optimised_directory
{
  my ($self) = @_;
  return File::Spec->catfile($self->{output_directory}, 'hgap_4_0_assembly');
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
  return 1;
}


=head2 get_parameters

 Title   : get_parameters
 Usage   : my $params = $obj->get_parameters();
 Function: get parameters. Extend later to read the hgap.log file
 Returns : hash containing parameters about final assembly

=cut
sub get_parameters
{
	my ($self, $filename) = @_;
	return 1;
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
  return 'hgap';
}

1;
