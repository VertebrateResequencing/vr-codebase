
=head1 NAME

VertRes::Pipelines::AnnotateAssembly -  Pipeline for annotating an assembly

=head1 SYNOPSIS

# make the config files, which specifies the details for connecting to the
# VRTrack g1k-meta database and the data roots:
echo '__VRTrack_AnnotateAssembly__ annotation.conf' > pipeline.config
# where annotation.conf contains:
root    => '/abs/path/to/root/data/dir',
module  => 'VertRes::Pipelines::AnnotateAssembly',
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
    assembler => 'velvet',
    annotation_tool   => 'Prokka',
    dbdir => '/lustre/scratch118/infgen/pathogen/pathpipe/prokka',
    tmp_directory => '/tmp',
    pipeline_version => 1
}

# by default __VRTrack_AnnotateAssembly__ will pick up lanes that have been both mapped
# and qcd. To override this, set eg:
# vrtrack_processed_flags => { mapped => 1, stored => 0 },

# run the pipeline:
run-pipeline -c pipeline.config

# (and make sure it keeps running by adding that last to a regular cron job)

=head1 DESCRIPTION

Pipeline for annotating an assembly

=head1 AUTHOR

path-help@sanger.ac.uk

=cut

package VertRes::Pipelines::AnnotateAssembly;

use strict;
use warnings;
use VRTrack::VRTrack;
use VRTrack::Lane;
use VRTrack::Library;
use VRTrack::Sample;
use VertRes::LSF;
use base qw(VertRes::Pipeline);
use VertRes::Utils::FileSystem;
use VertRes::Utils::Assembly;
use File::Spec;
use Utils;

our $actions = [
    {
        name     => 'annotate_assembly',
        action   => \&annotate_assembly,
        requires => \&annotate_assembly_requires,
        provides => \&annotate_assembly_provides
    },
    {
        name     => 'update_db',
        action   => \&update_db,
        requires => \&update_db_requires,
        provides => \&update_db_provides
    },
    {
        name     => 'cleanup',
        action   => \&cleanup,
        requires => \&cleanup_requires,
        provides => \&cleanup_provides
    }
];

our %options = ( bsub_opts => '' );

sub new {
    my ( $class, @args ) = @_;

    my $self = $class->SUPER::new( %options, actions => $actions, @args );
    if ( defined( $self->{db} ) ) {
        $self->{vrtrack} = VRTrack::VRTrack->new( $self->{db} ) or $self->throw("Could not connect to the database\n");
    }
    $self->{fsu} = VertRes::Utils::FileSystem->new;
    $self->{assembly_file}           = $self->_assembly_for_lane();
    $self->{assembly_base_directory} = $self->_assembly_base_directory_for_lane();


    return $self;
}

=head2 _genus_of_lane

 Title   : _genus_of_lane
 Usage   : my $genus = $obj->_genus_of_lane($vrlane, $vrtrack);
 Function: Given a lane, lookup the database and return the Genus, derived from the Species name
 Returns : String
 Args    : VRTrack::Lane object, VRTrack instance

=cut
sub _genus_of_lane
{
  my ($self,$vrlane, $vrtrack) = @_;
  my $individual = $self->_individual_of_lane($vrlane, $vrtrack);
  return undef unless(defined($individual));
  my $species = $individual->species;
  return undef unless(defined($species));
  return $species->genus();
}

sub _sample_accession_of_lane
{
  my ($self,$vrlane, $vrtrack) = @_;
  my $individual = $self->_individual_of_lane($vrlane, $vrtrack);
  return undef unless(defined($individual));
  return $individual->acc;
}

sub _individual_of_lane
{
   my ($self,$vrlane, $vrtrack) = @_;
   my $lib = VRTrack::Library->new($vrtrack, $vrlane->library_id);
   return undef unless(defined($lib));
   my $sample = VRTrack::Sample->new($vrtrack, $lib->sample_id);
   return undef unless(defined($sample));
   return $sample->individual;
}


sub _assembly_object {
    my ($self) = @_;

    my $assembly_util = VertRes::Utils::Assembly->new( assembler => $self->{assembler} );
    my $assembler_class = $assembly_util->find_module();
    eval "require $assembler_class;";
    my $assembler = $assembler_class->new(
        assembler        => $self->{assembler},
        output_directory => $self->{lane_path},
    );
    return $assembler;
}

sub _assembly_for_lane {
    my ($self) = @_;
    return $self->_assembly_object->optimised_assembly_file_path();
}

sub _assembly_base_directory_for_lane {
    my ($self) = @_;
    return $self->_assembly_object->optimised_directory();
}

sub _annotation_base_name {
    my ($self) = @_;
    return $self->{vrlane}->name();
}

sub _annotation_base_directory {
    my ($self) = @_;
    return join( '/', ( $self->_assembly_base_directory_for_lane, 'annotation' ) );
}

sub _annotation_base_name_with_path {
    my ($self) = @_;
    return join( '/', ( $self->_annotation_base_directory(), $self->_annotation_base_name ) );
}

sub annotate_assembly_provides {
    my $self = shift;
    my @provided_files ;

    for my $final_file (@{$self->cleanup_provides()})
    {
     unless(-e $final_file)
     {
	   push(@provided_files, join( '.', ( $self->_annotation_base_name_with_path, 'gff' ) ) );
  	   return \@provided_files;
     }
    }   

    return $self->cleanup_provides();
}

sub annotate_assembly_requires {
    my ($self) = @_;
    return [ $self->{assembly_file} ];
}

=head2 annotate_assembly

 Title   : annotate_assembly
 Usage   : $obj->annotate_assembly('/path/to/lane', 'lock_filename');
 Function: Take an assembly from the assembly pipeline and automatically annotate it.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub annotate_assembly {
    my ($self, $lane_path, $action_lock) = @_;
    
    my $memory_in_mb = $self->{memory} || 5000;
    
    my $lane_name = $self->{vrlane}->name;
    my $genus = $self->_genus_of_lane($self->{vrlane}, $self->{vrtrack});
    my $sample_accession = $self->_sample_accession_of_lane($self->{vrlane}, $self->{vrtrack}) || '';
    my $umask    = $self->umask_str;
    my $pipeline_version = join('/',($self->_annotation_base_directory,'pipeline_version_'.$self->{pipeline_version}));
    my $kingdom = $self->{kingdom} || "Bacteria";
    
      my $script_name = $self->{fsu}->catfile($lane_path, $self->{prefix}."annotate_assembly.pl");
      open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
      print $scriptfh qq{
  use strict;
  use Bio::AutomatedAnnotation;
  $umask

  my \$obj = Bio::AutomatedAnnotation->new(
    assembly_file => qq[$self->{assembly_file}],
    annotation_tool  => qq[$self->{annotation_tool}],
    sample_name      => qq[$lane_name],
    accession_number => qq[$sample_accession],
    dbdir            => qq[$self->{dbdir}],
    tmp_directory    => qq[$self->{tmp_directory}],
    genus            => qq[$genus],
    outdir           => "annotation",
    kingdom          => qq[$kingdom]
  );
  \$obj->annotate;
  
  system("touch $pipeline_version");
  system("touch $self->{prefix}annotate_assembly_done");
  exit;
      };
      close $scriptfh;

    my $job_name = $self->{prefix}.'annotate_assembly';
      
    #$self->delete_bsub_files($lane_path, $job_name);   #commenting this out to allow reruns with increased resources
    VertRes::LSF::run($action_lock, $lane_path, $job_name, {bsub_opts => "-M${memory_in_mb} -R 'select[mem>$memory_in_mb] rusage[mem=$memory_in_mb]'"}, qq{perl -w $script_name});

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
    return $self->annotate_assembly_provides();
}

=head2 update_db_provides

 Title   : update_db_provides
 Usage   : my $provided_files = $obj->update_db_provides('/path/to/lane');
 Function: Find out what files the update_db action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub update_db_provides {
    my $self = shift;
    my @provided_files ;

    for my $final_file (@{$self->cleanup_provides()})
    {
     unless(-e $final_file)
     {
  	   push(@provided_files, $self->{lane_path}."/".$self->{prefix}."annotate_update_db_done");
  	   return \@provided_files;
     }
    }   

    return $self->cleanup_provides();
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
    return $$self{'Yes'} unless(defined($vrlane));
    
    my $vrtrack = $vrlane->vrtrack;
    
    unless($vrlane->is_processed('annotated')){
      $vrtrack->transaction_start();
      $vrlane->is_processed('annotated',1);
      $vrlane->update() || $self->throw("Unable to set annotated status on lane $lane_path");
      $vrtrack->transaction_commit();
    }

    
    my $job_status =  File::Spec->catfile($lane_path, $self->{prefix} . 'job_status');
    Utils::CMD("rm $job_status") if (-e $job_status);
    Utils::CMD("touch ".$self->{fsu}->catfile($lane_path,"$self->{prefix}annotate_update_db_done")   );  

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
  return [ $self->{lane_path}."/".$self->{prefix}."annotate_update_db_done"];
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
	
    return [$self->{lane_path}."/".$self->{prefix}."annotate_cleanup_done"];
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
  
  my $prefix = $self->{prefix};
  # remove job files
  for my $file (qw(annotate_assembly )) 
    {
      foreach my $suffix (qw(o e pl)) 
      {
        unlink($self->{fsu}->catfile($lane_path, $prefix.$file.'.'.$suffix));
      }
  }
  
  for my $file ($self->_annotation_base_directory().'/'.$self->{vrlane}->name.'.tbl', $self->{lane_path}."/".$self->{prefix}."annotate_update_db_done", $self->{lane_path}."/".$self->{prefix}."annotate_assembly_done")
  {
	  unlink($file) if(-e $file);
  }
  
  Utils::CMD("touch ".$self->{fsu}->catfile($lane_path,"$self->{prefix}annotate_cleanup_done")   );  
  $self->update_file_permissions($lane_path);
  return $self->{Yes};
}

1;

