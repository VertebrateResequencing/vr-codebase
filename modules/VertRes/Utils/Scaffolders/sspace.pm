=head1 NAME

VertRes::Utils::Scaffolders::sspace - scaffolder utility functions

=head1 SYNOPSIS

use VertRes::Utils::Scaffolders::sspace;

my $scaffolder_util = VertRes::Utils::Scaffolders::sspace->new();

# use any of the utility functions described here, eg.
my $scaffolder = $scaffolder_util->new(
  scaffolder_exec => 'sspace'
);
=

=head1 DESCRIPTION

sspace-specific scaffolder functions

=head1 AUTHOR

path-help@sanger.ac.uk

=cut

package VertRes::Utils::Scaffolders::sspace;

use strict;
use warnings;
use VertRes::Wrapper::smalt;

use base qw(VertRes::Utils::Scaffold);

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Utils::Scaffolders::sspace->new();
 Function: Create a new VertRes::Utils::Scaffolders::sspace object.
 Returns : VertRes::Utils::Scaffolders::sspace object
 Args    : n/a

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    return $self;
}

sub create_lib_file
{
  my ($self) = @_;
  my $input_files = join(' ', @{$self->{input_files}} );
  open(LIB_FILE, "+>".$self->{output_directory}."/_scaffolder.config");
  print LIB_FILE "LIB ".$input_files." ".$self->get_insert_size()." 0.3 0";
  close(LIB_FILE);
}

sub scaffold
{
   my ($self) = @_;
   $self->create_lib_file();
   my $config_file = $self->{output_directory}."/_scaffolder.config";
   $self->run_sspace($self->{assembled_file_directory}, $self->{assembled_file}, $config_file, shift(@{$self->{merge_sizes}}), 31, 0);
   
   for my $merge_size (@{$self->{merge_sizes}})
   {
      $self->run_sspace($self->{output_directory}, 'intermediate.scaffolded.fasta', $config_file, $merge_size, 31, 0);
   }
  
  system("mv intermediate.scaffolded.fasta scaffolded_contigs.fa");
}

sub run_sspace
{
   my ($self, $source_directory, $source_name, $config_file, $k, $n, $x) = @_;
   `perl $self->{scaffolder_exec} -l $config_file -n $n -s $source_directory/$source_name -x $x -k $k -b scaffolded`;
   system("mv scaffolded.final.scaffolds.fasta intermediate.scaffolded.fasta");
   $self->cleanup();
}

# assumes that mapping and bamcheck has been run.
sub get_insert_size
{
  my ($self) = @_;
  open(BAMCHECK, $self->{assembled_file_directory}."/contigs.mapped.sorted.bam.bc") or $self->throw("couldnt find bamcheck file to get insert size for scaffolding");
  while(<BAMCHECK>)
  {
    my $line = $_;
    if($line =~ m/insert size average:\t([\d]+)/)
    {
      return $1;
    }
  }
  # sensible default
  return 350;
}


sub cleanup
{
  my ($self) = @_;
  system("rm -rf reads/ bowtieoutput/ pairinfo/ intermediate_results/"); 
}

sub name
{
  my ($self) = @_;
  return 'sspace';
}

1;
