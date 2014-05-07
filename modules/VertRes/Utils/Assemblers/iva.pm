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
  if (defined($self->{adapeters_file}) and defined($self->{trimmomatic_jar})) {
    $trimming_opts = "--trimmo $self->{trimmomatic_jar} --adapters $self->{adapeters_file}";
  }

  `$self->{optimiser_exec} $trimming_opts --fr $self->{files_str} --threads $num_threads iva_assembly`;
  print "\n-------------\nRUN IVA: $self->{optimiser_exec} $trimming_opts --fr $self->{files_str} --threads $num_threads iva_assembly\n---------\n";
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
