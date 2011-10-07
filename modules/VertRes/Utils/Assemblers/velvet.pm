=head1 NAME

VertRes::Utils::Assemblers::velvet - assembly utility functions

=head1 SYNOPSIS

use VertRes::Utils::Assemblers::velvet;

my $assembly_util = VertRes::Utils::Assemblers::velvet->new();

# use any of the utility functions described here, eg.
my $assember = $assembler_class->new(
  assembler_exec => '/path/to/velvet', 
  optimiser_exec => '/path/to/VelvetOptimiser.pl',
  min_kmer => 29, 
  max_kmer => 49
);

my \$ok = \$assember->optimise_parameters(

  );

=head1 DESCRIPTION

velvet-specific assembly functions

=head1 AUTHOR

path-help@sanger.ac.uk

=cut

package VertRes::Utils::Assemblers::velvet;

use strict;
use warnings;
use VertRes::Wrapper::velvet;

use base qw(VertRes::Utils::Assembly);

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Utils::Assemblers::velvet->new();
 Function: Create a new VertRes::Utils::Assemblers::velvet object.
 Returns : VertRes::Utils::Assemblers::velvet object
 Args    : n/a

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    return $self;
}

sub optimise_parameters
{
  my ($class, @args) = @_;
  `perl $self->{optimiser_exec} -s $self->{min_kmer} -e $self->{max_kmer} -p '_optimised_parameters' -f $self->{files_str} `;
  
  return 1;
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
    $files_str .= ' -'.$lane_pool->{type}." -fastq.gz $directory/pool_$pool_count.fastq.gz";
  }
  return $files_str;
}

1;
