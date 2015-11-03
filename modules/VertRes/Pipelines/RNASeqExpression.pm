=head1 NAME

VertRes::Pipelines::RNASeqExpression -  Pipeline for calculating expression for RNA Seq experiments. Outputs CSV with RPKM values for each gene,
Artemis plot files of the read coverage (broken down into forward and reverse), and it unsets the markduplicates flag in the BAM, replacing the original BAM.

=head1 SYNOPSIS

# make the config files, which specifies the details for connecting to the
# VRTrack g1k-meta database and the data roots:
echo '__VRTrack_RNASeqExpression__ rna_seq_expression.conf' > pipeline.config
# where rna_seq_expression.conf contains:
root    => '/abs/path/to/root/data/dir',
module  => 'VertRes::Pipelines::RNASeqExpression',
prefix  => '_',


 
 
limit => 50,
db  => {
    database => 'g1k_meta',
    host     => 'mcs4a',
    port     => 3306,
    user     => 'vreseq_rw',
    password => 'xxxxxxx',
}

data => {
    sequencing_file_suffix      => '.bam',
    protocol  => "StandardProtocol",
    annotation_filename => "my_reference.gff",
    mapping_quality => 30,
    bitwise_flag => 2,
    mpileup_cmd => 'samtools mpileup',
    window_margin => 50,
    intergenic_regions => 1,
}

# by default __VRTrack_RNASeqExpression__ will pick up lanes that have been both mapped
# and qcd. To override this, set eg:
# vrtrack_processed_flags => { mapped => 1, stored => 0 },

# run the pipeline:
run-pipeline -c pipeline.config

# (and make sure it keeps running by adding that last to a regular cron job)

=head1 DESCRIPTION

Pipeline for calculating expression for RNA Seq experiments.  Outputs CSV with RPKM values for each gene,
Artemis plot files of the read coverage (broken down into forward and reverse), and it unsets the markduplicates flag in the BAM, replacing the original BAM.


=head1 AUTHOR

path-help@sanger.ac.uk

=cut

package VertRes::Pipelines::RNASeqExpression;

use strict;
use warnings;
use VRTrack::VRTrack;
use VRTrack::Lane;
use VertRes::LSF;
use File::Basename;
use VertRes::Parser::bam;
use File::Basename;
use Bio::RNASeq;
use Utils;
use Data::Dumper;

use base qw(VertRes::Pipeline);

our $actions = [ { name     => 'calculate_expression',
                   action   => \&calculate_expression,
                   requires => \&calculate_expression_requires, 
                   provides => \&calculate_expression_provides },
                 { name     => 'update_db',
                   action   => \&update_db,
                   requires => \&update_db_requires, 
                   provides => \&update_db_provides },
                 { name     => 'cleanup',
                   action   => \&cleanup,
                   requires => \&cleanup_requires,
                   provides => \&cleanup_provides }
                   ];

our %options = (bsub_opts => ''
 );



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

sub _find_sequencing_files
{
  my $self = shift;
  opendir(DIR,$self->{lane_path});
  my $sequencing_file_suffix = $self->{sequencing_file_suffix} || ".bam";
  my @sequencing_filenames = grep { /$sequencing_file_suffix$/i }
  readdir(DIR);
  closedir(DIR);

  my $filtered_sequencing_filenames = $self->_filter_out_sym_links(\@sequencing_filenames);
  
  return $filtered_sequencing_filenames;
}

sub _filter_out_sym_links
{
  my $self = shift;
  my $sequencing_filenames = shift;
  my @filtered_sequencing_filenames;
  for my $filename (@{$sequencing_filenames})
  {
    next if( -l $self->{lane_path}."/".$filename);
    push( @filtered_sequencing_filenames, $filename);
  }
  return \@filtered_sequencing_filenames;
}


sub calculate_expression_provides
{
  my $self = shift;

  my @expression_done_files ;
  for my $filename (@{$self->_find_sequencing_files})
  {
    next if($self->does_reference_in_bam_match_annotation($self->{lane_path}."/".$filename,$self->{annotation_file}) == 0);
    
    push(@expression_done_files, $self->{lane_path}."/".$self->{prefix}.$filename."_calculate_expression_done");
  }

  return \@expression_done_files;
}

sub calculate_expression_requires
{
  my $self = shift;
  return [];
}

sub _create_expression_job
{
  my ($self, $build_path,$action_lock, $sequencing_filename) = @_;
  my $output_directory = $self->{lane_path};
  my $umask    = $self->umask_str;

  my $job_name = $self->{prefix}.$sequencing_filename.'_calculate_expression';
  my $script_name = $self->{fsu}->catfile($output_directory, $self->{prefix}.$sequencing_filename.'_calculate_expression.pl');
  my $prefix = $self->{prefix};
  
  my $total_memory_mb = 5000;
  my $queue = 'normal';
  if($self->get_reference_size_from_bam($output_directory.'/'.$sequencing_filename) > 6000000)
  {
    # Large reference so have some sensible defaults to get it to run in a reasonable amount of time
    $total_memory_mb = 20000;
    $self->{parallel_processes} ||= 16;
    $self->{intergenic_regions} ||= 0;
    $self->{no_coverage_plots}  ||= 1;
    $queue = 'long';
  }
  else
  {
    # Small reference so can run more intensive tasks by default
    $self->{intergenic_regions} ||= 1;
    $self->{parallel_processes} ||= 4;
    $self->{no_coverage_plots}  ||= 0;
  }
  
  my($action_lock_filename, $directories, $suffix) = fileparse($action_lock);
  my $sequencing_file_action_lock = $self->{lane_path}.'/'.$self->{prefix}.'calculate_expression.jids';

  my $mpileup_str  = "";
  if(defined ($self->{mpileup_cmd}))
  {
    $mpileup_str = ' mpileup => "'.$self->{mpileup_cmd}.'", ';
  }
  my $window_margin_str = "";
  if(defined ($self->{window_margin}))
  {
    $window_margin_str = ' window_margin => '.$self->{window_margin}.', ';
  }


  my $intergenic_regions_str = "";
  if(defined ($self->{intergenic_regions}))
    {
      $intergenic_regions_str = ' intergenic_regions => '.$self->{intergenic_regions}.', ';
    }

  my $filters = {
		 mapping_quality => $self->{mapping_quality}
		};
	

  if ( defined ($self->{bitwise_flag}) )
    {
      $filters->{bitwise_flag} = $self->{bitwise_flag};
    }

  $Data::Dumper::Terse = 1;

  my $filters_for_script = Dumper($filters);

  my $plots_class = "CoveragePlot";

  my $gtf_file = $self->{annotation_file};
  $gtf_file =~ s!gff$!gtf!i;
  my $feature_counts_cmd = '';
  if(-e $gtf_file)
  {
    my $output_counts_name = $sequencing_filename.".featurecounts.csv";
    $feature_counts_cmd =  'system( "featureCounts -O -T '.$self->{parallel_processes}.' -t exon -g gene_id -a '.$gtf_file.'   -o  '.$output_counts_name.' '.$sequencing_filename.'");';
  }

  open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");

  print $scriptfh qq{
  use strict;
  use Bio::RNASeq;
  use Bio::RNASeq::$plots_class;
  $umask
  
$feature_counts_cmd

  my \$expression_results = Bio::RNASeq->new(
    sequence_filename    => qq[$sequencing_filename],
    annotation_filename  => qq[$self->{annotation_file}],
    filters => $filters_for_script,
    protocol             => qq[$self->{protocol}],
    output_base_filename => qq[$sequencing_filename],
    parallel_processes   => $self->{parallel_processes},
    $window_margin_str
    $intergenic_regions_str
    );

  \$expression_results->output_spreadsheet();
};

  if((!defined($self->{no_coverage_plots})) || defined($self->{no_coverage_plots}) && $self->{no_coverage_plots} != 1 )
  {
        print $scriptfh qq{
        Bio::RNASeq::$plots_class->new(
          filename             => \$expression_results->_corrected_sequence_filename,
          output_base_filename => qq[$sequencing_filename],
          mapping_quality      => $self->{mapping_quality},
          $mpileup_str
        )->create_plots();
      };
  }
  print $scriptfh qq{
  system('touch $prefix${sequencing_filename}_calculate_expression_done');
  exit;
                };
                close $scriptfh;

        VertRes::LSF::run($sequencing_file_action_lock, $output_directory, $job_name, {bsub_opts => "-q $queue -n".$self->{parallel_processes}." -M${total_memory_mb} -R 'span[hosts=1] select[mem>$total_memory_mb] rusage[mem=$total_memory_mb]'", dont_wait=>1}, qq{perl -w $script_name});

        # we've only submitted to LSF, so it won't have finished; we always return
        # that we didn't complete
        return $self->{No};
  
}


sub calculate_expression
{
  my ($self, $build_path,$action_lock) = @_;
  for my $sequencing_filename (@{$self->_find_sequencing_files})
  {
    next if( -e $self->{lane_path}."/".$self->{prefix}.$sequencing_filename."_calculate_expression_done");
    next if($self->does_reference_in_bam_match_annotation($self->{lane_path}."/".$sequencing_filename,$self->{annotation_file}) == 0);
    
    $self->_create_expression_job($build_path,$action_lock, $sequencing_filename);
  }
  return $self->{No};
}

sub get_reference_size_from_bam
{
  my($self, $sequencing_filename) = @_;
  
  my $total_reference_size = 0;
  my $obj = VertRes::Parser::bam->new(file => $sequencing_filename);
  my %all_sequences_info = $obj->sequence_info();
  my @sequence_names = keys(%all_sequences_info);
  
  for my $sequence_line (@sequence_names)
  {
     my $sequence_size = $obj->sequence_info($sequence_line, 'LN');
     next unless(defined($sequence_size));
     $total_reference_size += $sequence_size;
  }
  
  return $total_reference_size;  
}

sub get_reference_from_bam
{
  my($self, $sequencing_filename) = @_;
  
  my $obj = VertRes::Parser::bam->new(file => $sequencing_filename);
  my %all_sequences_info = $obj->sequence_info();
  my @sequence_names = keys(%all_sequences_info);
  my $reference_file = $obj->sequence_info($sequence_names[0], 'UR');
  return undef unless(defined($reference_file));
  $reference_file =~ s/file://;
  return $reference_file;
}

# Assuming the name of the annotation has the same base filename as the reference, we can tell if they matches. Depends on a consistent naming scheme.
sub does_reference_in_bam_match_annotation
{
  my($self, $sequencing_filename, $annotation_filename) = @_;
  my $reference_file = $self->get_reference_from_bam($sequencing_filename);
  return 0 if(! defined($reference_file));
  
   my($reference_basefilename, $directories, $suffix) = fileparse($reference_file, qr/\.[^.]*/);
  
  if($annotation_filename =~ /$reference_basefilename/)
  {
    return 1;
  }
 
  return 0;
}

=head2 update_db_requires

 Title   : update_db_requires
 Usage   : my $required_files = $obj->update_db_requires('/path/to/lane');
 Function: Find out what files the update_db action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub update_db_requires {
    my ($self, $lane_path) = @_;
    return $self->calculate_expression_provides();
}

=head2 update_db_provides

 Title   : update_db_provides
 Usage   : my $provided_files = $obj->update_db_provides('/path/to/lane');
 Function: Find out what files the update_db action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub update_db_provides {
   my ($self) = @_;
    return [ $self->{lane_path}."/".$self->{prefix}."rna_seq_update_db_done"];
}

=head2 update_db

 Title   : update_db
 Usage   : $obj->update_db('/path/to/lane', 'lock_filename');
 Function: Records in the database that the lane has been improved.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut
sub update_db {
    my ($self, $lane_path, $action_lock) = @_;

    # all the done files are there, so just need to update the processed
    # flag in the database (if not already updated)
    my $vrlane = $self->{vrlane};
    my $vrtrack = $vrlane->vrtrack;
    
    return $$self{'Yes'} if $vrlane->is_processed('rna_seq_expression');

    unless($vrlane->is_processed('rna_seq_expression')){
      $vrtrack->transaction_start();
      $vrlane->is_processed('rna_seq_expression',1);
      $vrlane->update() || $self->throw("Unable to set rna_seq_expression status on lane $lane_path");
      $vrtrack->transaction_commit();
    }

    
    my $job_status =  File::Spec->catfile($lane_path, $self->{prefix} . 'job_status');
    Utils::CMD("rm $job_status") if (-e $job_status);
    Utils::CMD("touch ".$self->{fsu}->catfile($lane_path,$self->{prefix}."rna_seq_update_db_done")   );  

    return $$self{'Yes'};
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
  return [ $self->{lane_path}."/".$self->{prefix}."rna_seq_update_db_done"];
}

=head2 cleanup_provides

 Title   : cleanup_provides
 Usage   : my $provided_files = $obj->cleanup_provides('/path/to/lane');
 Function: Find out what files the cleanup action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub cleanup_provides {
  my ($self) = @_;
    return [ $self->{lane_path}."/".$self->{prefix}."rna_seq_cleanup_done"];
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
  foreach my $file (@{$self->_find_sequencing_files()}) 
    {
      foreach my $suffix (qw(o e pl)) 
      {
        unlink($self->{fsu}->catfile($lane_path, $prefix.$file.'_calculate_expression.'.$suffix));
      }
  }
  
  Utils::CMD("touch ".$self->{fsu}->catfile($lane_path,$prefix."rna_seq_cleanup_done")   );  
  $self->update_file_permissions($lane_path);
  return $self->{Yes};
}

1;
