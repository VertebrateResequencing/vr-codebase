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

use base qw(VertRes::Pipeline);

our $actions = [ { name     => 'pool_fastqs',
                   action   => \&pool_fastqs,
                   requires => \&pool_fastqs_requires, 
                   provides => \&pool_fastqs_provides },
                # { name     => 'optimise_parameters',
                #   action   => \&optimise_parameters,
                #   requires => \&optimise_parameters_requires, 
                #   provides => \&optimise_parameters_provides },
                # { name     => 'run_assembler',
                #   action   => \&run_assembler,
                #   requires => \&run_assembler_requires, 
                #   provides => \&run_assembler_provides },
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
  
  return $self;
}


sub pool_fastqs
{
      my ($self, $build_path) = @_;
    
      my $lane_names = $self->get_all_lane_names($self->{pools});
      my $output_directory = '.';
      

      my $script_name = $self->{fsu}->catfile($build_path, $self->{prefix}."pool_fastqs.pl");
      open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
      print $scriptfh qq{
  use strict;
  use VertRes::Pipelines::Assembly;
  my \$assembly= VertRes::Pipelines::Assembly->new();
  
};

   for my $lane_name ( @$lane_names)
   {
     my $lane_path = $self->{vrtrack}->hierarchy_path_of_lane_name($lane_name);
     print $scriptfh qq{\$assembly->shuffle_sequences_fastq_gz("$self->{root}/$lane_name", "$lane_path", "$output_directory"); 
     };
   }
   
   my $pool_count = 1;
   for my $lane_pool (@{$self->{pools}})
   {
    my $lane_names_str = '("'.join('.fastq.gz","',@{$lane_pool->{lanes}}).'.fastq.gz")';
    print $scriptfh qq{
      my \@lane_names = $lane_names_str;
      \$assembly->concat_fastq_gz_files(\@lane_names, "pool_$pool_count.fastq.gz", "$output_directory", "$output_directory");
     };
     $pool_count++;
   }
 
   print $scriptfh qq{
     
     system("touch _pool_fastqs_done");
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
   ["_pool_fastqs_done"];
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

  while(<$FILEA>) {
    print $OUTFILE $_;
    $_ = <$FILEA>;
    print $OUTFILE $_;
    $_ = <$FILEA>;
    print $OUTFILE $_;
    $_ = <$FILEA>;
    print $OUTFILE $_;

    $_ = <$FILEB>;
    print $OUTFILE $_;
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

