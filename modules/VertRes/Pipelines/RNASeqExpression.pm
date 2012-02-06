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
    mpileup_cmd => 'samtools mpileup'
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
use LSF;
use File::Basename;

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

  my $job_name = $self->{prefix}.$sequencing_filename.'_calculate_expression';
  my $script_name = $self->{fsu}->catfile($output_directory, $self->{prefix}.$sequencing_filename.'_calculate_expression.pl');
  my $total_memory_mb = 3000;
  
  my($action_lock_filename, $directories, $suffix) = fileparse($action_lock);
  my $sequencing_file_action_lock = $self->{lane_path}.$sequencing_filename.$action_lock_filename;

  my $mpileup_str  = "";
  if(defined ($self->{mpileup_cmd}))
  {
    $mpileup_str = ' mpileup => "'.$self->{mpileup_cmd}.'", ';
  }

        open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
        print $scriptfh qq{
  use strict;
  use Pathogens::RNASeq::Expression;
  
  
  my \$expression_results = Pathogens::RNASeq::Expression->new(
    sequence_filename    => qq[$sequencing_filename],
    annotation_filename  => qq[$self->{annotation_file}],
    mapping_quality      => $self->{mapping_quality},
    protocol             => qq[$self->{protocol}],
    output_base_filename => qq[$sequencing_filename]
    );

  \$expression_results->output_spreadsheet();
  
  Pathogens::RNASeq::CoveragePlot->new(
    filename             => \$expression_results->_corrected_sequence_filename,
    output_base_filename => qq[$sequencing_filename],
    $mpileup_str
  )->create_plots();

  system('touch _${sequencing_filename}_calculate_expression_done');
  exit;
                };
                close $scriptfh;

        LSF::run($sequencing_file_action_lock, $output_directory, $job_name, {bsub_opts => "-M${total_memory_mb}000 -R 'select[mem>$total_memory_mb] rusage[mem=$total_memory_mb]'", dont_wait=>1}, qq{perl -w $script_name});

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
    
    $self->_create_expression_job($build_path,$action_lock, $sequencing_filename);
  }
  return $self->{No};
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
    return [ $self->{lane_path}."/".$self->{prefix}."update_db_done"];
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
    Utils::CMD("touch ".$self->{fsu}->catfile($lane_path,"_update_db_done")   );  

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
  return [ $self->{lane_path}."/".$self->{prefix}."update_db_done"];
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
    return [ $self->{lane_path}."/".$self->{prefix}."cleanup_done"];
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
  
  Utils::CMD("touch ".$self->{fsu}->catfile($lane_path,"_cleanup_done")   );  
  
  return $self->{Yes};
}

1;
