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
  
  #SPAdes needs a comma separated list of kmer values
  my @spades_kmer_array;
  my $kmer_value = $self->{min_kmer};
  while ($kmer_value < $self->{max_kmer}){
	push(@spades_kmer_array, $kmer_value);
	$kmer_value = $kmer_value + 4; #Step in 4
  }
  my $spades_kmer_string = join (',', @spades_kmer_array);	
  
  `python $self->{optimiser_exec} --12 $self->{files_str} --only-assembler -k $spades_kmer_string -o spades_assembly`;
  
  my $params = $self->get_parameters("spades.log"); 
  system("mv  $params->{assembly_directory} ".$self->optimised_directory());
  system("touch ".$self->optimised_directory()."/_spades_optimise_parameters_done");
   
  unlink($self->optimised_directory()."/params.txt");
  unlink($self->optimised_directory()."/dataset.info");
  unlink($self->optimised_directory()."/contigs.fasta");
  #Should the spades.log file be deleted too? For now, I call it spades_assembly_logfile.txt
  system("mv $self->optimised_directory()/spades.log $self->optimised_directory()/spades_assembly_logfile.txt");

  return 1;
}                                                                                                                             

sub optimised_directory
{
  my ($self) = @_;
  return "$self->{output_directory}/spades_assembly";
}

sub optimised_assembly_file_path
{
  my ($self) = @_;
  return join('/',($self->optimised_directory(),'/scaffolds.fa'));
}

sub map_and_generate_stats
{
   my ($self, $directory, $output_directory, $lane_paths) = @_;
   
   my $forward_fastq = '';
   my $reverse_fastq = '';
   
   for my $lane_path ( @$lane_paths)
   {
     my ($base_directory,$base,$suff) = Utils::basename($lane_path);
     opendir(my $lane_dir_handle, $base_directory);
     my @fastq_files  = grep { /\.fastq\.gz$/ } readdir($lane_dir_handle);
     if(@fastq_files >=1 )
     {
       $forward_fastq .= $base_directory.'/'.$fastq_files[0];
     }
     if(@fastq_files >=2 )
     {
       $reverse_fastq .= $base_directory.'/'.$fastq_files[1];
     }
   }

   unless( -e "$output_directory/forward.fastq")
   {
     `gzip -cd $forward_fastq  > $output_directory/forward.fastq`;
   }
   unless(-e "$output_directory/reverse.fastq")
   {
     `gzip -cd $reverse_fastq  > $output_directory/reverse.fastq`;
   }

   my $reference = $self->{reference} || "$directory/contigs.fa";

   my $mapper = VertRes::Wrapper::smalt->new();
   $mapper->setup_custom_reference_index($reference,'-k 13 -s 4','small');

   `smalt map -x -i 3000 -f samsoft -y 0.95 -o $directory/contigs.mapped.sam $reference.small $output_directory/forward.fastq $output_directory/reverse.fastq`;
   $self->throw("Sam file not created") unless(-e "$directory/contigs.mapped.sam");

   `samtools faidx $reference`;
   $self->throw("Reference index file not created") unless(-e "$reference.fai");

   `samtools view -bt $reference.fai $directory/contigs.mapped.sam > $directory/contigs.mapped.bam`;
   $self->throw("Couldnt convert from sam to BAM") unless(-e "$directory/contigs.mapped.bam");
   unlink("$directory/contigs.mapped.sam");

   `samtools sort -m 500000000 $directory/contigs.mapped.bam $directory/contigs.mapped.sorted`;
   $self->throw("Couldnt sort the BAM") unless(-e "$directory/contigs.mapped.sorted.bam");

   `samtools index $directory/contigs.mapped.sorted.bam`;
   $self->throw("Couldnt index the BAM") unless(-e "$directory/contigs.mapped.sorted.bam.bai");

   `bamcheck -c 1,20000,5 -r $reference $directory/contigs.mapped.sorted.bam >  $directory/contigs.mapped.sorted.bam.bc`;
   `plot-bamcheck -s $reference > $reference.gc`;
   `plot-bamcheck -p $directory/qc_graphs/ -r  $reference.gc $directory/contigs.mapped.sorted.bam.bc`;
   
   $self->generate_stats($directory);
   unlink("$directory/contigs.mapped.bam");
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
	my %parameters;
	$parameters{assembly_directory} = "$self{output_directory}/spades";
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

  my $memory_required = 0.5 * $input_params->{total_number_of_reads} + 1500000;

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
