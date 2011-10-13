=head1 NAME

VertRes::Pipelines::Assembly - Assemble genomes

=head1 SYNOPSIS

echo '__Assembly__ assembly.conf' > pipeline.config
# where assembly.conf contains:
root    => '/abs/path/to/root/data/dir',
module  => 'VertRes::Pipelines::Assembly',
prefix  => '_',

db  => {
    database => 'pathogen_prok_track',
    host     => 'mcs6',
    port     => 3347,
    user     => 'pathpipe_rw',
    password => 'xxxx',
},

data => {
    db  => {
        database => 'pathogen_prok_track',
        host     => 'mcs6',
        port     => 3347,
        user     => 'pathpipe_rw',
        password => 'xxxx',
    },

    assembler => 'velvet',
    assembler_exec => '',
    optimiser_exec => '/software/pathogen/external/apps/usr/bin/VelvetOptimiser.pl',
    # rough estimate so we know how much RAM to request
    genome_size => 10000000,
    pools => [
        {
            lanes => ['123_4','456_7#8'],
            type => 'shortPaired',
        },
        {
            lanes => ['543_2','4987_1#7'],
            type => 'longPaired',
        },
    ],
},


=head1 DESCRIPTION


=head1 AUTHOR

path-help@sanger.ac.uk

=cut
package VertRes::Pipelines::Assembly;

use strict;
use warnings;
use VertRes::Utils::Hierarchy;
use VertRes::IO;
use VertRes::Utils::FileSystem;
use VRTrack::VRTrack;
use VRTrack::Lane;
use VRTrack::File;
use File::Basename;
use Time::Format;
use LSF;
use Data::Dumper;
use FileHandle;
use VertRes::Utils::Assembly;

use base qw(VertRes::Pipeline);

our $actions = [ { name     => 'pool_fastqs',
                   action   => \&pool_fastqs,
                   requires => \&pool_fastqs_requires,
                   provides => \&pool_fastqs_provides },
                 { name     => 'optimise_parameters',
                   action   => \&optimise_parameters,
                   requires => \&optimise_parameters_requires,
                   provides => \&optimise_parameters_provides },
                { name     => 'map_back',
                  action   => \&map_back,
                  requires => \&map_back_requires, 
                  provides => \&map_back_provides },
                { name     => 'optimise_parameters_with_reference',
                  action   => \&optimise_parameters_with_reference,
                  requires => \&optimise_parameters_with_reference_requires,
                  provides => \&optimise_parameters_with_reference_provides },
                { name     => 'map_back_with_reference',
                  action   => \&map_back_with_reference,
                  requires => \&map_back_with_reference_requires, 
                  provides => \&map_back_with_reference_provides },
                  
                # { name     => 'statistics',
                #   action   => \&statistics,
                #   requires => \&statistics_requires,
                #   provides => \&statistics_provides },
                # { name     => 'cleanup',
                #   action   => \&cleanup,
                #   requires => \&cleanup_requires,
                #   provides => \&cleanup_provides }
                ];

our %options = (
                do_cleanup => 1);

sub new {
  my ($class, @args) = @_;

  my $self = $class->SUPER::new(%options, actions => $actions, @args);
  if(defined($self->{db}))
  {
    $self->{vrtrack} = VRTrack::VRTrack->new($self->{db}) or $self->throw("Could not connect to the database\n");
  }
  $self->{fsu} = VertRes::Utils::FileSystem->new;

  unless(defined($self->{max_threads}))
  {
    $self->{max_threads} = 8;
  }

  my $assembly_util = VertRes::Utils::Assembly->new(assembler => $self->{assembler});
  my $assembler_class = $assembly_util->find_module();
  eval "require $assembler_class;";
  $self->{assembler_class} = $assembler_class;

  return $self;
}

###########################
# Begin map_back
###########################

sub map_back_requires
{
  my ($self) = @_;
  return $self->optimise_parameters_provides();
}

sub map_back_provides
{
   my ($self) = @_;

   return [$self->{lane_path}."/".$self->{prefix}.$self->{assembler}.'_plot_bamcheck_done'];

}

sub map_back
{
  my ($self, $build_path) = @_;
  my $assembler_class = $self->{assembler_class};
  my $output_directory = $self->{lane_path};
  eval("use $assembler_class; ");
  my $assembler_util= $assembler_class->new( output_directory => qq[$output_directory]);
  my $base_path = $self->{lane_path}.'/../../seq-pipelines';
  
  my $job_name = $self->{prefix}.$self->{assembler}."_map_back";
  my $script_name = $self->{fsu}->catfile($output_directory, $self->{prefix}.$self->{assembler}."_map_back.pl");
  
  my $lane_names = $self->get_all_lane_names($self->{pools});
  my @lane_paths;
  for my $lane_name (@$lane_names)
  {
    push(@lane_paths,$base_path.'/'.$self->{vrtrack}->hierarchy_path_of_lane_name($lane_name).'/'.$lane_name);
  }
  my $lane_paths_str = '("'.join('","', @lane_paths).'")';

  open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
  print $scriptfh qq{
  use strict;
  use $assembler_class;
  use VertRes::Wrapper::smalt;
    
  my \$forward_fastq = '';
  my \$reverse_fastq = '';
  my \@lane_paths = $lane_paths_str;
  
  my \$assembler_util= $assembler_class->new( output_directory => qq[$output_directory]);

  for my \$lane_path ( \@lane_paths)
  {
    \$forward_fastq .= \$lane_path.'_1.fastq.gz ';
    \$reverse_fastq .= \$lane_path.'_2.fastq.gz ';
  }
  
  unless( -e "$output_directory/forward.fastq")
  {
    `gzip -cd \$forward_fastq  > $output_directory/forward.fastq`;
  }
  unless(-e "$output_directory/reverse.fastq")
  {
    `gzip -cd \$reverse_fastq  > $output_directory/reverse.fastq`;
  }
  
  my \$directory = \$assembler_util->optimised_directory();
 
  next unless(-e "\$directory/_$self->{assembler}_optimise_parameters_done");
  next if( -e "\$directory/_$self->{assembler}_plot_bamcheck_done");
  my \$mapper = VertRes::Wrapper::smalt->new();
  \$mapper->setup_reference("\$directory/contigs.fa");
  
  `smalt map -x -i 3000 -f samsoft -o \$directory/contigs.mapped.sam \$directory/contigs.fa.small $output_directory/forward.fastq $output_directory/reverse.fastq`;
  \$assembler_util->throw("Sam file not created") unless(-e "\$directory/contigs.mapped.sam");
  
  `samtools faidx \$directory/contigs.fa`;
  \$assembler_util->throw("Reference index file not created") unless(-e "\$directory/contigs.fa.fai");
  
  `samtools view -bt \$directory/contigs.fa.fai \$directory/contigs.mapped.sam > \$directory/contigs.mapped.bam`;
  \$assembler_util->throw("Couldnt convert from sam to BAM") unless(-e "\$directory/contigs.mapped.bam");
  unlink("\$directory/contigs.mapped.sam");
  
  `samtools sort \$directory/contigs.mapped.bam \$directory/contigs.mapped.sorted`;
  \$assembler_util->throw("Couldnt sort the BAM") unless(-e "\$directory/contigs.mapped.sorted.bam");
  
  `samtools index \$directory/contigs.mapped.sorted.bam`;
  \$assembler_util->throw("Couldnt index the BAM") unless(-e "\$directory/contigs.mapped.sorted.bam.bai");
 
  `bamcheck \$directory/contigs.mapped.sorted.bam >  \$directory/contigs.mapped.sorted.bam.bc`;
  
  `plot-bamcheck -p \$directory/qc_graphs/ \$directory/contigs.mapped.sorted.bam.bc`;
  \$assembler_util->generate_stats(\$directory);
  
  unlink("\$directory/contigs.mapped.bam");
  
  system("touch \$directory/_$self->{assembler}_plot_bamcheck_done");
 
  system("touch _$self->{assembler}_plot_bamcheck_done");
  exit;
                };
  close $scriptfh;
  
  my $action_lock = "$output_directory/$$self{'prefix'}$self->{assembler}_map_back.jids";
  my $memory_required_mb = 5000;

  LSF::run($action_lock, $output_directory, $job_name, {bsub_opts => " -M${memory_required_mb}000 -R 'select[mem>$memory_required_mb] rusage[mem=$memory_required_mb]'"}, qq{perl -w $script_name});

  # we've only submitted to LSF, so it won't have finished; we always return
  # that we didn't complete
  return $self->{No};
  
}

###########################
# End map_back
###########################



###########################
# Begin map_back_with_reference
###########################

sub map_back_with_reference_requires
{
  my ($self) = @_;
  return $self->optimise_parameters_with_reference_provides();
}

sub map_back_with_reference_provides
{
   my ($self) = @_;

   return [$self->{lane_path}."/".$self->{prefix}.$self->{assembler}.'_plot_bamcheck_with_reference_done'];

}

sub map_back_with_reference
{
  my ($self, $build_path) = @_;
  my $assembler_class = $self->{assembler_class};
  my $output_directory = $self->{lane_path};
  eval("use $assembler_class; ");
  my $assembler_util= $assembler_class->new( output_directory => qq[$output_directory]);
  my $base_path = $self->{lane_path}.'/../../seq-pipelines';
  
  my $job_name = $self->{prefix}.$self->{assembler}."_map_back_with_reference";
  my $script_name = $self->{fsu}->catfile($output_directory, $self->{prefix}.$self->{assembler}."_map_back_with_reference.pl");
  
  my $lane_names = $self->get_all_lane_names($self->{pools});
  my @lane_paths;
  for my $lane_name (@$lane_names)
  {
    push(@lane_paths,$base_path.'/'.$self->{vrtrack}->hierarchy_path_of_lane_name($lane_name).'/'.$lane_name);
  }
  my $lane_paths_str = '("'.join('","', @lane_paths).'")';

  open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
  print $scriptfh qq{
  use strict;
  use $assembler_class;
  use VertRes::Wrapper::smalt;
    
  my \$forward_fastq = '';
  my \$reverse_fastq = '';
  my \@lane_paths = $lane_paths_str;
  
  my \$assembler_util= $assembler_class->new( output_directory => qq[$output_directory]);

  for my \$lane_path ( \@lane_paths)
  {
    \$forward_fastq .= \$lane_path.'_1.fastq.gz ';
    \$reverse_fastq .= \$lane_path.'_2.fastq.gz ';
  }
  
  unless( -e "$output_directory/forward.fastq")
  {
    `gzip -cd \$forward_fastq  > $output_directory/forward.fastq`;
  }
  unless(-e "$output_directory/reverse.fastq")
  {
    `gzip -cd \$reverse_fastq  > $output_directory/reverse.fastq`;
  }
  
  my \$directory = \$assembler_util->optimised_with_reference_directory();
 
  next unless(-e "\$directory/_$self->{assembler}_optimise_parameters_with_reference_done");
  next if( -e "\$directory/_$self->{assembler}_plot_bamcheck_with_reference_done");
  my \$mapper = VertRes::Wrapper::smalt->new();
  \$mapper->setup_reference("\$directory/contigs.fa");
  
  `smalt map -x -i 3000 -f samsoft -o \$directory/contigs.mapped.sam \$directory/contigs.fa.small $output_directory/forward.fastq $output_directory/reverse.fastq`;
  \$assembler_util->throw("Sam file not created") unless(-e "\$directory/contigs.mapped.sam");
  
  `samtools faidx \$directory/contigs.fa`;
  \$assembler_util->throw("Reference index file not created") unless(-e "\$directory/contigs.fa.fai");
  
  `samtools view -bt \$directory/contigs.fa.fai \$directory/contigs.mapped.sam > \$directory/contigs.mapped.bam`;
  \$assembler_util->throw("Couldnt convert from sam to BAM") unless(-e "\$directory/contigs.mapped.bam");
  unlink("\$directory/contigs.mapped.sam");
  
  `samtools sort \$directory/contigs.mapped.bam \$directory/contigs.mapped.sorted`;
  \$assembler_util->throw("Couldnt sort the BAM") unless(-e "\$directory/contigs.mapped.sorted.bam");
  
  `samtools index \$directory/contigs.mapped.sorted.bam`;
  \$assembler_util->throw("Couldnt index the BAM") unless(-e "\$directory/contigs.mapped.sorted.bam.bai");
 
  `bamcheck \$directory/contigs.mapped.sorted.bam >  \$directory/contigs.mapped.sorted.bam.bc`;
  
  `plot-bamcheck -p \$directory/qc_graphs/ \$directory/contigs.mapped.sorted.bam.bc`;
  \$assembler_util->generate_stats(\$directory);
  
  unlink("\$directory/contigs.mapped.bam");
  
  system("touch \$directory/_$self->{assembler}_plot_bamcheck_with_reference_done");
 
  system("touch _$self->{assembler}_plot_bamcheck_with_reference_done");
  exit;
                };
  close $scriptfh;
  
  my $action_lock = "$output_directory/$$self{'prefix'}$self->{assembler}_map_back_with_reference.jids";
  my $memory_required_mb = 5000;

  LSF::run($action_lock, $output_directory, $job_name, {bsub_opts => " -M${memory_required_mb}000 -R 'select[mem>$memory_required_mb] rusage[mem=$memory_required_mb]'"}, qq{perl -w $script_name});

  # we've only submitted to LSF, so it won't have finished; we always return
  # that we didn't complete
  return $self->{No};
  
}

###########################
# End map_back_with_reference
###########################



###########################
# Begin optimise
###########################

sub optimise_parameters_provides
{
  my $self = shift;

  return  [$self->{lane_path}."/".$self->{prefix}."$self->{assembler}_optimise_parameters_done"];
}

sub optimise_parameters_requires
{
  my $self = shift;
  return $self->pool_fastqs_provides();
}

sub optimise_parameters
{
      my ($self, $build_path) = @_;

      my $lane_names = $self->get_all_lane_names($self->{pools});
      my $output_directory = $self->{lane_path};

      my $assembler_class = $self->{assembler_class};
      my $optimiser_exec = $self->{optimiser_exec};

      eval("use $assembler_class; ");
      my $assembler_util= $assembler_class->new();
      my $files_str = $assembler_util->generate_files_str($self->{pools}, $output_directory);

      my $job_name = $self->{prefix}.$self->{assembler}.'_optimise_parameters';
      my $script_name = $self->{fsu}->catfile($output_directory, $self->{prefix}.$self->{assembler}."_optimise_parameters.pl");

      my $kmer = $self->calculate_kmer_size();
      
      my $action_lock = "$output_directory/$$self{'prefix'}".$self->{assembler}."_optimise_parameters.jids";
      
      my $memory_required_mb = $self->estimate_memory_required($output_directory, $kmer->{min})/1000;
      my $queue = 'long';
      if($memory_required_mb > 35000)
      {
        $queue = 'hugemem';
      }
      elsif($memory_required_mb < 3000)
      {
        $queue = 'normal';
      }

      my $num_threads = $self->number_of_threads($memory_required_mb);

      open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
      print $scriptfh qq{
use strict;
use $assembler_class;

my \$assembler = $assembler_class->new(
  assembler => qq[$self->{assembler}],
  optimiser_exec => qq[$optimiser_exec],
  min_kmer => $kmer->{min},
  max_kmer => $kmer->{max},
  files_str => qq[$files_str],
  output_directory => qq[$output_directory],
  );

my \$ok = \$assembler->optimise_parameters($num_threads);

\$assembler->throw("optimising parameters for assembler failed - try again?") unless( -e \$assembler->optimised_directory()."/_$self->{assembler}_optimise_parameters_done");
system('touch _$self->{assembler}_optimise_parameters_done');
exit;
              };
              close $scriptfh;

      my $total_memory_mb = $num_threads*$memory_required_mb;
      
      LSF::run($action_lock, $output_directory, $job_name, {bsub_opts => "-n$num_threads -q $queue -M${total_memory_mb}000 -R 'select[mem>$total_memory_mb] rusage[mem=$total_memory_mb] span[hosts=1]'", dont_wait=>1}, qq{perl -w $script_name});

      # we've only submitted to LSF, so it won't have finished; we always return
      # that we didn't complete
      return $self->{No};
}

sub lane_read_length
{
  my ($self) = @_;
  my $lane_names = $self->get_all_lane_names($self->{pools});
  my $vrlane  = VRTrack::Lane->new_by_name($self->{vrtrack}, @$lane_names[0]) or $self->throw("No such lane in the DB: [".@$lane_names[0]."]");
  my $read_length = $vrlane->read_len() || 75;
  return $read_length;
}

# check kmers between 66% and 90% of read size. Choose lane and find out read size
sub calculate_kmer_size
{
  my ($self) = @_;
  my %kmer_size;
  my $read_length = $self->lane_read_length();

  $kmer_size{min} = int($read_length*0.66);
  $kmer_size{max} = int($read_length*0.90);
  
  if($kmer_size{min} % 2 == 0)
  {
    $kmer_size{min}--;
  }
  
  if($kmer_size{max} % 2 == 0)
  {
    $kmer_size{max}--;
  }

  return \%kmer_size;
}

sub number_of_threads
{
  my ($self, $memory_required_mb) = @_;
  my $normal_queue_mem_limit = 35000;
  my $num_threads = 1;
  
  if($normal_queue_mem_limit/$memory_required_mb > 2)
  {
    $num_threads = int($normal_queue_mem_limit/$memory_required_mb);
  }
  if($num_threads  > $self->{max_threads})
  {
    $num_threads = $self->{max_threads};
  }
  
  return $num_threads;
}

###########################
# End optimise
###########################



###########################
# Begin optimise_parameters_with_reference
###########################

sub optimise_parameters_with_reference_provides
{
  my $self = shift;
  return  [$self->{lane_path}."/".$self->{prefix}."$self->{assembler}_optimise_parameters_with_reference_done"];
}

sub optimise_parameters_with_reference_requires
{
  my $self = shift;
  return $self->map_back_provides();
}

sub optimise_parameters_with_reference
{
      my ($self, $build_path) = @_;

      my $lane_names = $self->get_all_lane_names($self->{pools});
      my $output_directory = $self->{lane_path};

      my $assembler_class = $self->{assembler_class};
      my $optimiser_exec = $self->{optimiser_exec};

      eval("use $assembler_class; ");
      my $assembler_util= $assembler_class->new();
      my $files_str = $assembler_util->generate_files_str($self->{pools}, $output_directory);

      my $job_name = $self->{prefix}.$self->{assembler}.'_optimise_parameters_with_reference';
      my $script_name = $self->{fsu}->catfile($output_directory, $self->{prefix}.$self->{assembler}."_optimise_parameters_with_reference.pl");

      my $kmer = $self->calculate_kmer_size();
      
      my $action_lock = "$output_directory/$$self{'prefix'}".$self->{assembler}."_optimise_parameters_with_reference.jids";
      
      my $memory_required_mb = $self->estimate_memory_required($output_directory, $kmer->{min})/1000;
      my $queue = 'long';
      if($memory_required_mb > 35000)
      {
        $queue = 'hugemem';
      }
      elsif($memory_required_mb < 3000)
      {
        $queue = 'normal';
      }

      my $num_threads = $self->number_of_threads($memory_required_mb);

      open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
      print $scriptfh qq{
use strict;
use $assembler_class;

my \$assembler = $assembler_class->new(
  assembler => qq[$self->{assembler}],
  optimiser_exec => qq[$optimiser_exec],
  min_kmer => $kmer->{min},
  max_kmer => $kmer->{max},
  files_str => qq[$files_str],
  output_directory => qq[$output_directory],
  );

my \$ok = \$assembler->optimise_parameters_with_reference($num_threads);

\$assembler->throw("optimising parameters for assembler failed - try again?") unless( -e \$assembler->optimise_parameters_with_reference()."/_$self->{assembler}_optimise_parameters_with_reference_done");
system('touch _$self->{assembler}_optimise_parameters_with_reference_done');
exit;
              };
              close $scriptfh;

      my $total_memory_mb = $num_threads*$memory_required_mb;
      
      LSF::run($action_lock, $output_directory, $job_name, {bsub_opts => "-n$num_threads -q $queue -M${total_memory_mb}000 -R 'select[mem>$total_memory_mb] rusage[mem=$total_memory_mb] span[hosts=1]'", dont_wait=>1}, qq{perl -w $script_name});

      # we've only submitted to LSF, so it won't have finished; we always return
      # that we didn't complete
      return $self->{No};
}


###########################
# End optimise_parameters_with_reference
###########################


###########################
# Begin pool fastqs
###########################


sub pool_fastqs
{
      my ($self, $build_path) = @_;

      my $lane_names = $self->get_all_lane_names($self->{pools});
      my $output_directory = $self->{lane_path};
      my $base_path = $self->{lane_path}.'/../../seq-pipelines';

      my $script_name = $self->{fsu}->catfile($build_path, $self->{prefix}."pool_fastqs.pl");
      open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
      print $scriptfh qq{
use strict;
use VertRes::Pipelines::Assembly;
my \$assembly= VertRes::Pipelines::Assembly->new();
my \@lane_names;
};

   for my $lane_name ( @$lane_names)
   {
     my $lane_path = $self->{vrtrack}->hierarchy_path_of_lane_name($lane_name);
     print $scriptfh qq{\$assembly->shuffle_sequences_fastq_gz("$lane_name", "$base_path/$lane_path", "$output_directory");
};
   }

   my $pool_count = 1;
   for my $lane_pool (@{$self->{pools}})
   {
    my $lane_names_str = '("'.join('.fastq.gz","',@{$lane_pool->{lanes}}).'.fastq.gz")';
    print $scriptfh qq{
\@lane_names = $lane_names_str;
\$assembly->concat_fastq_gz_files(\\\@lane_names, "pool_$pool_count.fastq.gz", "$output_directory", "$output_directory");
};
     for my $lane_name (  @{$lane_pool->{lanes}}  )
     {
        print $scriptfh qq{
unlink("$output_directory/$lane_name.fastq.gz");
};
     }

     $pool_count++;
   }

   print $scriptfh qq{
system("touch $output_directory/_pool_fastqs_done");
exit;
      };
      close $scriptfh;
      my $action_lock = "$output_directory/$$self{'prefix'}pool_fastqs.jids";
      my $job_name = $self->{prefix}.'pool_fastqs';

      LSF::run($action_lock, $output_directory, $job_name, {bsub_opts => '-M500000 -R \'select[mem>500] rusage[mem=500]\''}, qq{perl -w $script_name});

      # we've only submitted to LSF, so it won't have finished; we always return
      # that we didn't complete
      return $self->{No};
}

sub pool_fastqs_requires
{
  [];
}

sub pool_fastqs_provides
{
   my $self = shift;
   my @provided_files ;
   push(@provided_files, $self->{lane_path}.'/'.$self->{prefix}."pool_fastqs_done");

   my $pool_count = 1;
   for my $lane_pool (@{$self->{pools}})
   {
     push(@provided_files, $self->{lane_path}."/pool_$pool_count.fastq.gz");
   }

   return \@provided_files;
}

sub get_all_lane_names
{
  my ($self, $pooled_lanes) = @_;
  my @all_lane_names ;

  for my $lane_pool (@$pooled_lanes)
  {
    for my $lane_name (@{$lane_pool->{lanes}})
    {
      push(@all_lane_names, $lane_name);
    }
  }
  return \@all_lane_names;
}

sub concat_fastq_gz_files
{
  my ($self, $filenames, $outputname, $input_directory, $output_directory) = @_;

  my $filenames_with_paths = '';
  for my $filename (@$filenames)
  {
    $filenames_with_paths = $filenames_with_paths.$input_directory.'/'.$filename.' ';
  }
  `gzip -cd $filenames_with_paths | gzip > $output_directory/$outputname`;
}


# adapted from https://github.com/dzerbino/velvet/blob/master/shuffleSequences_fastq.pl
sub shuffle_sequences_fastq_gz
{
  my ($self, $name, $input_directory, $output_directory) = @_;

  my $filenameA = $name."_1.fastq.gz";
  my $filenameB = $name."_2.fastq.gz";
  my $filenameOut = $name.".fastq.gz";

  open( my $FILEA, "-|",'gunzip -c '.$input_directory.'/'.$filenameA) ;
  open( my $FILEB, "-|",'gunzip -c '.$input_directory.'/'.$filenameB) ;
  open( my $OUTFILE, "|-", "gzip -c  > $output_directory".'/'."$filenameOut");

  # FIXME: if FileB contains more reads than fileA they will get missed! This is also a bug in Velvets shuffleSequences_fastq.pl
  while(<$FILEA>) {
    print $OUTFILE $_;
    $_ = <$FILEA>;
    print $OUTFILE $_;
    $_ = <$FILEA>;
    print $OUTFILE $_;
    $_ = <$FILEA>;
    print $OUTFILE $_;

    my $file_b_line = <$FILEB>;
    next unless(defined($file_b_line));
    print $OUTFILE $file_b_line;
    $_ = <$FILEB>;
    print $OUTFILE $_;
    $_ = <$FILEB>;
    print $OUTFILE $_;
    $_ = <$FILEB>;
    print $OUTFILE $_;
  }
  close($FILEA);
  close($FILEB);
  close($OUTFILE);
}

###########################
# End pool fastqs
###########################



sub total_number_of_reads
{
  my ($self) = @_;
  my $lane_names = $self->get_all_lane_names($self->{pools});
  my $total_reads = 0;

  for my $lane_name (@{$lane_names})
  {
    my $vrlane  = VRTrack::Lane->new_by_name($self->{vrtrack}, $lane_name) or $self->throw("No such lane in the DB: [".$lane_name."]");
    $total_reads += $vrlane->raw_reads();
  }
  return $total_reads;
}


#Gives the answer in kb.
sub estimate_memory_required
{
  my ($self, $output_directory, $kmer_size) = @_;

  my %memory_params ;
  $memory_params{total_number_of_reads} = $self->total_number_of_reads();
  $memory_params{genome_size}           = $self->{genome_size};
  $memory_params{read_length}           = $self->lane_read_length();
  $memory_params{kmer_size}             = $kmer_size;

  my $assembler_class = $self->{assembler_class};
  eval("use $assembler_class; ");
  my $assembler_util= $assembler_class->new(output_directory => $output_directory);
  my $memory_required_in_kb = $assembler_util->estimate_memory_required(\%memory_params);
  return $memory_required_in_kb;
}




sub cleanup
{
 # files
 #_pool_fastqs.o
 #_pool_fastqs.e
 #_pool_fastqs.pl
 #_pool_fastqs.jids
 #_optimise_parameters.pl
 #_optimise_parameters.jids
 #_optimise_parameters.o
 #_optimise_parameters.e
 #_run_assembler.e
 #_run_assembler.pl
 #_run_assembler.jids
 #_run_assembler.o

 # directories
 #_optimise_parameters_data_31
 # some files in assembly_xxxxxxxxxxx   Sequences PreGraph Graph2
}

