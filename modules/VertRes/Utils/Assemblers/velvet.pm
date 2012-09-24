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
 Function: estimate the momory required for the assembler in KB
 Returns : integer in kb of memory requirement
 Ram required for velvetg = -109635 + 18977*ReadSize + 86326*GenomeSize + 233353*NumReads - 51092*K
 http://listserver.ebi.ac.uk/pipermail/velvet-users/2009-July/000474.html
=cut
sub estimate_memory_required
{
  my ($self, $input_params) = @_;
  my $optimised_params;
  my $kmer_size = $input_params->{kmer_size};
  unless(defined($kmer_size))
  {
    $optimised_params = $self->get_parameters("velvet_assembly_logfile.txt");
    $kmer_size = $optimised_params->{kmer};
  }

  my $memory_required = -109635 + (18977*($input_params->{read_length})) + (86326*($input_params->{genome_size})/1000000) + (233353*($input_params->{total_number_of_reads})/1000000) - (51092*$kmer_size);
  $memory_required *= 2.0;
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
