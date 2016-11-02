=head1 NAME

VertRes::Utils::Assemblers::spades - assembly utility functions

=head1 SYNOPSIS

use VertRes::Utils::Assemblers::spades;

my $assembly_util = VertRes::Utils::Assemblers::spades->new();

# use any of the utility functions described here, eg.
my $assember = $assembler_class->new(
  assembler_exec => '/path/to/spades',
  optimiser_exec => '/path/to/spades', #The spades assembler is like velvetoptimiser in that it runs the assembler for various kmer sizes and then chooses the best assembly
  min_kmer => 29, #odd number
  max_kmer => 101, #odd number
);

my \$ok = \$assember->optimise_parameters(

  );

=head1 DESCRIPTION

spades-specific assembly functions

=head1 AUTHOR

path-help@sanger.ac.uk

=cut

package VertRes::Utils::Assemblers::spades;

use strict;
use warnings;
use VertRes::Wrapper::smalt;
use Utils;

use base qw(VertRes::Utils::Assembly);

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Utils::Assemblers::spades->new();
 Function: Create a new VertRes::Utils::Assemblers::spades object.
 Returns : VertRes::Utils::Assemblers::spades object
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
  
	
  my $kmer_string = $self->_create_kmer_values_string();
	if(defined($self->{spades_kmer_opts}) and ($self->{spades_kmer_opts}) )
	{
		$kmer_string = $self->{spades_kmer_opts};
	}

  my $spades_opts = '';
	if(defined($self->{spades_opts}) and ($self->{spades_opts}) )
	{
		$spades_opts = $self->{spades_opts};
	}
	else
	{
    if (defined($self->{single_cell}) and ($self->{single_cell})){
      $spades_opts = '--sc --careful';
    }
    else {
	  	# Only of use with viruses, dont use with bacteria
      $spades_opts = '--only-assembler';
    }
  }
  
  `python $self->{optimiser_exec} --12 $self->{files_str} $spades_opts --threads $num_threads -k $kmer_string -o spades_assembly`;
  
  my $params = $self->get_parameters("spades.log");

  # Delete any unwanted files that SPAdes produces
  unlink($self->optimised_directory()."/dataset.info");
  unlink($self->optimised_directory()."/contigs.fasta"); #These are the unscaffolded contigs which we are not interested in
  #Should the spades.log file be deleted too? For now, I call it spades_assembly_logfile.txt. It is quite a large file though.
  system("mv ".$self->optimised_directory()."/spades.log ".$self->optimised_directory()."/spades_assembly_logfile.txt");
  system("mv ".$self->optimised_directory()."/scaffolds.fasta ".$self->optimised_directory()."/contigs.fa"); #Move ths scaffolded sequences into a file called contigs.fa
  system("rm -rf ".$self->optimised_directory()."/K*"); #Directories with the intermediate assemblies
  system("rm -rf ".$self->optimised_directory()."/misc");
  system("rm -rf ".$self->optimised_directory()."/mismatch_corrector_contigs");
  system("rm -rf ".$self->optimised_directory()."/configs");
  system("touch ".$self->optimised_directory()."/_spades_optimise_parameters_done");

  return 1;
}

sub _create_kmer_values_string
{
  my ($self) = @_;
  #SPAdes needs a comma separated list of kmer values and all the kmer values should be less than 128
  my @spades_kmer_array;
  my $current_kmer_value = $self->{min_kmer};
  my $max_kmer_value = $self->{max_kmer};
  
  # We make the kmers go from 21-127 effectively. Change once we learn the optimum kmer range for SPAdes
  if($current_kmer_value > 41){
  	$current_kmer_value = 41;
  }
  
  if($max_kmer_value > 127){
  	$max_kmer_value = 127;
  }
  while ($current_kmer_value < $max_kmer_value){
	push(@spades_kmer_array, $current_kmer_value);
	$current_kmer_value = $current_kmer_value + 4; #Step in 4
  }
  my $spades_kmer_string = join (',', @spades_kmer_array);
  return $spades_kmer_string;
}


sub optimised_directory
{
  my ($self) = @_;
  return "$self->{output_directory}/spades_assembly";
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
    $files_str .= "$directory/pool_$pool_count.fastq.gz";
  }
  return $files_str;
}


=head2 get_parameters

 Title   : get_parameters
 Usage   : my $params = $obj->get_parameters();
 Function: get parameters. Extend later to read the spades.log file
 Returns : hash containing parameters about final assembly

=cut
sub get_parameters
{
	my ($self, $filename) = @_;
	my %parameters;
	$parameters{assembly_directory} = "$self->{output_directory}/spades_assembly";
	return \%parameters;
}


=head2 estimate_memory_required

 Title   : estimate_memory_required
 Usage   : my $memory_required_in_kb = $obj->estimate_memory_required();
 Function: estimate the memory required for the assembler in KB
 Returns : integer in kb of memory requirement
 Like for Velvet, Ram for SPAdes is estimated based on the total number of reads.
 The memory required for velvet shows a linear relation to the number of reads but with a lot of variation
 between different assemblies. The memory estimate is a rule of thumb based on the observed memory usage. We
 are yet to do extensive studies into the memory usage of SPAdes. For now, we use the same calculations as
 we do for Velvet.
=cut
sub estimate_memory_required
{
  my ($self, $input_params) = @_;

  my $memory_required = 0.4 * $input_params->{total_number_of_reads} + 1500000;
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

}

sub name
{
  my ($self) = @_;
  return 'spades';
}

1;
