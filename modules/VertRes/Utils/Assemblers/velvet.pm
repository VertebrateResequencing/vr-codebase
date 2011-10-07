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
  my ($self, @args) = @_;
  `perl $self->{optimiser_exec} -s $self->{min_kmer} -e $self->{max_kmer} -p '_optimised_parameters' -f '$self->{files_str}' `;
  
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

=head2 do_assembly

 Title   : do_assembly
 Usage   : my $module = $obj->do_assembly();
 Function: run the assembler
 Returns : 1 if successful

=cut
sub do_assembly
{
  my $self = shift;
  
  # get parameters from file
  my %assembler_parameters = $self->get_parameters();
  my $assembly_name = 'assembly_'.time();
  `$self->{assembler_exec}h $assembly_name $assembler_parameters->{velveth}`;
  `$self->{assembler_exec}g $assembly_name $assembler_parameters->{velvetg}`;

  1;
}


=head2 get_parameters

 Title   : get_parameters
 Usage   : my $params = $obj->get_parameters();
 Function: get parameters for the assembler programs. This comes from velvet optimiser
 Returns : hash containing parameters for each program to be run

=cut
sub get_parameters
{
   my $self = shift;
   open(PARAM_FILE, $self->{output_directory}.'/_optimised_parameters_logfile.txt') or die "Couldnt open optimised parameters file";
   
   my %parameters;
   my $found_final_assembly_details = 0; 
   while(<PARAM_FILE>)
   {
     $line = $_;
     if($found_final_assembly_details == 0 && $line =~ m/Final optimised assembly details/)
     {
       $found_final_assembly_details == 1;
       next;
     }
     
     if( $found_final_assembly_details == 1 && $line =~ m/Velveth parameter string: _optimised_parameters_data_[\d]+ ([\d]+ .+$)/)
     {
       $parameters{velveth} = $1;
     }
     
     if( $found_final_assembly_details == 1 && $line =~ m/Velvetg parameter string: _optimised_parameters_data_[\d]+ (.+$)/)
     {
       $parameters{velvetg} = $1;
     }
   }
   return \%parameters;
}


1;
