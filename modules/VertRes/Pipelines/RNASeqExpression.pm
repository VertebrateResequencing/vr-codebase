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
    filters => {mapping_quality => 30}
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
  $self->{lane_path}."/".
  
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


  my $expression_results = Pathogens::RNASeq::Expression->new(
    sequence_filename    => $sequence_file,
    annotation_filename  => $annotation_file,
    filters              => \%filters,
    protocol             => $protocols{$protocol_name},
    output_base_filename => $output_base_filename
    );


        open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
        print $scriptfh qq{
  use strict;
  use Pathogens::RNASeq::Expression;
  
  
  my \$expression_results = Pathogens::RNASeq::Expression->new(
    sequence_filename    => $sequencing_filename,
    annotation_filename  => $self->{annotation_file},
    filters              => $self->{filters},
    protocol             => $self->{protocol},
    output_base_filename => $sequencing_filename
    );

  \$expression_results->output_spreadsheet();

  system('touch _${sequencing_filename}_calculate_expression_done');
  exit;
                };
                close $scriptfh;

        my $total_memory_mb = $num_threads*$memory_required_mb;

        LSF::run($action_lock, $output_directory, $job_name, {bsub_opts => "-n$num_threads -q $queue -M${total_memory_mb}000 -R 'select[mem>$total_memory_mb] rusage[mem=$total_memory_mb] span[hosts=1]'", dont_wait=>1}, qq{perl -w $script_name});

        # we've only submitted to LSF, so it won't have finished; we always return
        # that we didn't complete
        return $self->{No};
  
}


sub calculate_expression
{
      my ($self, $build_path,$action_lock) = @_;

}















