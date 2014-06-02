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
use VertRes::Wrapper::smalt;
use Utils;

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
  my ($self, $num_threads) = @_;
  `perl $self->{optimiser_exec} -t $num_threads -s $self->{min_kmer} -e $self->{max_kmer} -p 'velvet_assembly' -f '$self->{files_str}' `;
  my $params = $self->get_parameters("velvet_assembly_logfile.txt");
  system("mv  $params->{assembly_directory} ".$self->optimised_directory());
  system("touch ".$self->optimised_directory()."/_velvet_optimise_parameters_done");
  
  unlink($self->optimised_directory()."/Sequences");
  unlink($self->optimised_directory()."/PreGraph");
  unlink($self->optimised_directory()."/Graph");
  unlink($self->optimised_directory()."/Graph2");
  unlink($self->optimised_directory()."/stats.txt");

  return 1;
}

sub optimised_directory
{
  my ($self) = @_;
  return "$self->{output_directory}/velvet_assembly";
}

sub optimised_assembly_file_path
{
  my ($self) = @_;
  return join('/',($self->optimised_directory(),'/contigs.fa'));
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


=head2 get_parameters

 Title   : get_parameters
 Usage   : my $params = $obj->get_parameters();
 Function: get parameters for the assembler programs. This comes from velvet optimiser
 Returns : hash containing parameters for each program to be run

=cut
sub get_parameters
{
   my ($self, $filename) = @_;
   open(PARAM_FILE, $self->{output_directory}.'/'.$filename) or die "Couldnt open optimised parameters file";
   
   my %parameters;
   my $found_final_assembly_details = 0;
   while(<PARAM_FILE>)
   {
     my $line = $_;
     if($found_final_assembly_details == 0 && $line =~ m/Final optimised assembly details/)
     {
       $found_final_assembly_details = 1;
       next;
     }
     
     if( $found_final_assembly_details == 1 && $line =~ m/Velveth parameter string: ([a-z_]+_([\d]+)) ([\d]+ .+$)/)
     {
       $parameters{assembly_directory} = $1;
       $parameters{kmer} = $2;
       $parameters{velveth} = $3;
     }
     
     if( $found_final_assembly_details == 1 && $line =~ m/Velvetg parameter string: [a-z_]+_[\d]+ (.+$)/)
     {
       $parameters{velvetg} = $1;
     }
   }
   return \%parameters;
}


=head2 estimate_memory_required

 Title   : estimate_memory_required
 Usage   : my $memory_required_in_kb = $obj->estimate_memory_required();
 Function: estimate the memory required for the assembler in KB
 Returns : integer in kb of memory requirement
 Ram for velvet is estimated based on the total number of reads.
 The memory required for velvet shows a linear relation to the number of reads but with a lot of variation
 between different assemblies. The memory estimate is a rule of thumb based on the observed memory usage.
=cut
sub estimate_memory_required
{
  my ($self, $input_params) = @_;

  my $memory_required = 0.5 * $input_params->{total_number_of_reads} + 1500000;
  $memory_required = $memory_required/2 if $input_params->{error_correct};

  if($memory_required < 1000000)
  {
    $memory_required = 1000000;
  }
  elsif($memory_required > 400000000)
  {
    $memory_required = 400000000;
  }
  
  return int($memory_required);
}

=head2 assembly_directories

 Title   : assembly_directories
 Usage   : my $module = $obj->assembly_directories();
 Function: Find out where the assemblies are located
 Returns : array of paths

=cut
sub assembly_directories
{
  my ($self) = @_;
  my @output_files ;
  opendir(DIR,$self->{output_directory});
  my @files = grep {/^velvet_assembly_/} readdir(DIR);
  for my $file (@files)
  {
    push(@output_files, $self->{output_directory}.'/'.$file) if(-d $self->{output_directory}.'/'.$file);
  }
  return \@output_files;
}

sub name
{
  my ($self) = @_;
  return 'velvet';
}

1;
