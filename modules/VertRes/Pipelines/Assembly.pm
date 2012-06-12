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
    num_threads => 4,
    scaffolder => 'sspace',
    scaffolder_exec => '/software/pathogen/external/apps/usr/local/SSPACE-1.2_linux-x86_64/SSPACE_v1-2.pl',
    
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
use VertRes::Utils::Scaffold;

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
                 { name     => 'scaffold',
                   action   => \&scaffold,
                   requires => \&scaffold_requires, 
                   provides => \&scaffold_provides },
                 { name     => 'map_back_scaffolded',
                   action   => \&map_back_scaffolded,
                   requires => \&map_back_scaffolded_requires, 
                   provides => \&map_back_scaffolded_provides },
                 { name     => 'cleanup',
                   action   => \&cleanup,
                   requires => \&cleanup_requires,
                   provides => \&cleanup_provides }
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
  
  
  my $scaffolder_util = VertRes::Utils::Scaffold->new(scaffolder => $self->{scaffolder});
  my $scaffolder_class = $scaffolder_util->find_module();
  eval "require $scaffolder_class;";
  $self->{scaffolder_class} = $scaffolder_class;
  

  return $self;
}


###########################
# Begin scaffold
###########################

sub scaffold_requires
{
  my ($self) = @_;
  return $self->map_back_provides();
}

sub scaffold_provides
{
   my ($self) = @_;
   return [$self->{lane_path}."/scaffolding_results/contigs.fa"];

}

sub scaffold
{
  my ($self, $build_path, $action_lock) = @_;
  $self->scaffold_with_sspace($build_path, $action_lock, "optimised_with_reference_directory", "scaffold");
  
  return $self->{No};
}

sub scaffold_with_sspace
{
  my ($self, $build_path, $action_lock, $working_directory_method_name, $action_name_suffix) = @_;
  
  my $scaffold_class = $self->{scaffolder_class};
  my $scaffold_exec = $self->{scaffolder_exec};
  my $output_directory = $self->{lane_path};
  my $assembler_class = $self->{assembler_class};
  
  eval("use $assembler_class;");
  my $assembler_util= $assembler_class->new( output_directory => qq[$output_directory]);
  my $working_directory = $assembler_util->${working_directory_method_name}();

  
  my $job_name = $self->{prefix}.$self->{assembler}."_$action_name_suffix";
  my $script_name = $self->{fsu}->catfile($output_directory, $self->{prefix}.$self->{assembler}."_$action_name_suffix.pl");
  
  # TODO: clean up output files, estimate memory usage
  
  open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
  print $scriptfh qq{
  use strict;
  use $scaffold_class;
  
  my \@input_files = ('forward.fastq', 'reverse.fastq');
  my \@merge_sizes = (90, 80, 70, 60, 50, 40, 30, 25, 20, 15, 10, 10, 10, 7, 7, 7, 5, 5, 5, 5 );
  my \$scaffolder_util= $scaffold_class->new( output_directory => qq[$output_directory],
    assembled_file => qq[contigs.fa],
    assembled_file_directory => qq[$working_directory],
    scaffolder_exec => qq[$self->{scaffolder_exec}],
    input_files => \\\@input_files,
    merge_sizes => \\\@merge_sizes
  );
  
  \$scaffolder_util->scaffold();
  
  if ( !-e "scaffolding_results" ) { Utils::CMD("mkdir -p scaffolding_results"); }
  system("mv scaffolded_contigs.fa scaffolding_results/contigs.fa");
  
  exit;
                };
  close $scriptfh;
  
  my $memory_required_mb = 10000;

  LSF::run($action_lock, $output_directory, $job_name, {bsub_opts => " -M${memory_required_mb}000 -R 'select[mem>$memory_required_mb] rusage[mem=$memory_required_mb]'"}, qq{perl -w $script_name});
  
}




###########################
# End scaffold
###########################


###########################
# Begin map_back
###########################

sub map_back_scaffolded_requires
{
  my ($self) = @_;
  return $self->scaffold_provides();
}

sub map_back_scaffolded_provides
{
   my ($self) = @_;

   return [$self->{lane_path}."/".$self->{prefix}.$self->{assembler}.'_map_back_scaffolded_done'];

}

sub map_back_scaffolded
{
  my ($self, $build_path, $action_lock) = @_;
  $self->mapping_and_generate_stats($build_path, $action_lock, "scaffolded_directory", "map_back_scaffolded");
  
  return $self->{No};
}


###########################
# End map_back
###########################




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

   return [$self->{lane_path}."/".$self->{prefix}.$self->{assembler}.'_map_back_done'];

}

sub map_back
{
  my ($self, $build_path, $action_lock) = @_;
  $self->mapping_and_generate_stats($build_path, $action_lock, "optimised_directory", "map_back",$self->{reference});
  
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

   return [$self->{lane_path}."/".$self->{prefix}.$self->{assembler}.'_map_back_with_reference_done'];

}

sub map_back_with_reference
{
  my ($self, $build_path, $action_lock) = @_;
  $self->mapping_and_generate_stats($build_path, $action_lock, "optimised_with_reference_directory", "map_back_with_reference");
  
  return $self->{No};
  
}

sub mapping_and_generate_stats
{
  my ($self, $build_path, $action_lock, $working_directory_method_name, $action_name_suffix, $reference) = @_;
  
  my $assembler_class = $self->{assembler_class};
  my $output_directory = $self->{lane_path};
  eval("use $assembler_class; ");
  my $assembler_util= $assembler_class->new( output_directory => qq[$output_directory]);
  my $base_path = $self->{lane_path}.'/../../seq-pipelines';
  
  my $job_name = $self->{prefix}.$self->{assembler}."_$action_name_suffix";
  my $script_name = $self->{fsu}->catfile($output_directory, $self->{prefix}.$self->{assembler}."_$action_name_suffix.pl");
  
  my $lane_names = $self->get_all_lane_names($self->{pools});
  my @lane_paths;
  for my $lane_name (@$lane_names)
  {
    push(@lane_paths,$base_path.'/'.$self->{vrtrack}->hierarchy_path_of_lane_name($lane_name).'/'.$lane_name);
  }
  my $lane_paths_str = '("'.join('","', @lane_paths).'")';
  my $reference_str = '';
  if($reference)
  {
    $reference_str = qq{,reference => qq[$reference]} ;
  }

  open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
  print $scriptfh qq{
  use strict;
  use $assembler_class;
    
  my \@lane_paths = $lane_paths_str;
  
  my \$assembler_util= $assembler_class->new( output_directory => qq[$output_directory] $reference_str );
  my \$directory = \$assembler_util->${working_directory_method_name}();
  \$assembler_util->map_and_generate_stats(\$directory,qq[$output_directory], \\\@lane_paths );
  
  system("touch \$directory/_$self->{assembler}_${action_name_suffix}_done");
 
  system("touch _$self->{assembler}_${action_name_suffix}_done");
  exit;
                };
  close $scriptfh;
  
  my $memory_required_mb = 8000;

  LSF::run($action_lock, $output_directory, $job_name, {bsub_opts => " -M${memory_required_mb}000 -R 'select[mem>$memory_required_mb] rusage[mem=$memory_required_mb]'"}, qq{perl -w $script_name});
  
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
      my ($self, $build_path,$action_lock) = @_;

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
      
      my $memory_required_mb = int($self->estimate_memory_required($output_directory, $kmer->{min})/1000);
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
      my ($self, $build_path, $action_lock) = @_;

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

      
      my $memory_required_mb = int($self->estimate_memory_required($output_directory, $kmer->{min})/1000);
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
      my ($self, $build_path,$action_lock) = @_;

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
     my $vlane = VRTrack::Lane->new_by_name($self->{vrtrack},$lane_name);
     my $file_names_str ;
     if( @{$vlane->files} == 2)
     {
       my @file_names ;
       for my $file_name (@{$vlane->files})
       {
         push(@file_names, $file_name->name );
       }
       $file_names_str = '("'.join('","',@file_names ).'")';
     }

     print $scriptfh qq{
       my \@filenames_array = $file_names_str;
       \$assembly->shuffle_sequences_fastq_gz("$lane_name", "$base_path/$lane_path", "$output_directory",\\\@filenames_array);
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
  my ($self, $name, $input_directory, $output_directory, $file_names) = @_;

  my $filenameA = $file_names->[0];
  my $filenameB = $file_names->[1];
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


=head2 cleanup_requires

 Title   : cleanup_requires
 Usage   : my $required_files = $obj->cleanup_requires('/path/to/lane');
 Function: Find out what files the cleanup action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub cleanup_requires {
  my ($self) = @_;
  return ["_velvet_map_back_scaffolded_done"];
}

=head2 cleanup_provides

 Title   : cleanup_provides
 Usage   : my $provided_files = $obj->cleanup_provides('/path/to/lane');
 Function: Find out what files the cleanup action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub cleanup_provides {
    return ["_cleanup_done"];
}

=head2 cleanup

 Title   : cleanup
 Usage   : $obj->cleanup('/path/to/lane', 'lock_filename');
 Function: Unlink all the pipeline-related files (_*) in a lane, as well
           as the split directory.
 Returns : $VertRes::Pipeline::Yes
 Args    : lane path, name of lock file to use

=cut

sub cleanup {
  my ($self, $lane_path, $action_lock) = @_;
#  return $self->{Yes} unless $self->{do_cleanup};
  
  my $prefix = $self->{prefix};
  
  # remove job files
  foreach my $file (qw(pool_fastqs 
    velvet_optimise_parameters 
    velvet_map_back 
    velvet_optimise_parameters_with_reference 
    velvet_map_back_with_reference
    velvet_scaffold
    velvet_map_back_scaffolded)) 
    {
      foreach my $suffix (qw(o e pl)) 
      {
        unlink($self->{fsu}->catfile($lane_path, $prefix.$file.'.'.$suffix));
      }
  }
  
  # remove files
  foreach my $file (qw(.RData contigs.fa.png.Rout scaffolded.summaryfile.txt scaffolded.logfile.txt _scaffolder.config scaffolded.final.evidence reverse.fastq forward.fastq pool_1.fastq.gz)) 
  {
    unlink($self->{fsu}->catfile($lane_path, $file));
  }
  Utils::CMD("touch ".$self->{fsu}->catfile($lane_path,"_cleanup_done")   );  
  
  return $self->{Yes};
}

