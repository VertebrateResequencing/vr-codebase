=head1 NAME

VertRes::Pipelines::Import_iRODS_fastq - pipeline for importing fastq files from iRODS


=head1 EXAMPLE CONFIG FILES:

=head2 pipeline.conf

__VRTrack_Import__ importirods.conf

=head2 import_irods.conf

root    => '/abs/path/to/root/data/dir',
module  => 'VertRes::Pipelines::Import_iRODS_fastq',
prefix  => '_',

db =>
   {
        database => 'pathogen_example_track',
        host     => 'web-mii-shap',
        port     => 3303,
        user     => 'pathpipe_ro',
        password => '',
   },
data => 
   {  # Pipeline-specific data
        db  => 
            {
               database => 'pathogen_example_track',
               host     => 'web-mii-shap',
               port     => 3303,
               user     => 'pathpipe_rw',
               password => 'xxx',       
            },
   }, 


=head1 DESCRIPTION

A module for importing fastq files into a tracking database. 

The module gets bam files from iRODS, converts them into fastq files and imports 
them into a tracking database.

The module is based on VertRes::Pipelines::Import_iRODS and uses the get_bams 
method from that module. It adapts a routine from VertRes::Pipelines::Mapping.pm 
to convert fastqs into bams using the bam2fastq method in VertRes::Utils::Sam. 
Bam files are deleted once fastq files have been written. Finally, the tracking 
database is updated using the method from VertRes::Pipelines::Import. 

=head1 AUTHOR

Craig Porter: cp7@sanger.ac.uk

=cut


package VertRes::Pipelines::Import_iRODS_fastq;
use base qw(VertRes::Pipelines::Import_iRODS);

use strict;
use warnings;
use LSF;
use VRTrack::VRTrack;
use VRTrack::Lane;
use VRTrack::File;
use VertRes::Utils::FileSystem;
use VertRes::Pipelines::Import;
use VertRes::Pipelines::Import_iRODS;

our @actions =
(
    # Create the hierarchy path, download and bamcheck the bam files.
    {
        'name'     => 'get_bams',
        'action'   => \&VertRes::Pipelines::Import_iRODS::get_bams,
        'requires' => \&get_bams_requires, 
        'provides' => \&get_bams_provides,
    },

    # Convert to fastq then delete bams.
    { 
	'name'     => 'bam_to_fastq',
	'action'   => \&bam_to_fastq,
	'requires' => \&bam_to_fastq_requires, 
        'provides' => \&bam_to_fastq_provides 
    },

    # If all files downloaded OK, update the VRTrack database.
    {
        'name'     => 'update_db',
        'action'   => \&VertRes::Pipelines::Import::update_db,
        'requires' => \&update_db_requires, 
        'provides' => \&update_db_provides,
    },
);

our $options = 
{
    'bamcheck'        => 'bamcheck -q 20',
    'bsub_opts'       => "-q normal -R 'select[type==X86_64] rusage[thouio=1]'",
};


# --------- OO stuff --------------

sub VertRes::Pipelines::Import_iRODS_fastq::new 
{
    my ($class, %args) = @_;
    my $self = $class->SUPER::new(%$options,'actions'=>\@actions,%args);

    $self->{fsu} = VertRes::Utils::FileSystem->new; 

    # Skip lane without updating db unless only bams in lane.
    foreach my $file (@{$$self{files}})
    {
	unless($file =~ /\.bam$/i)
	{
	    my $verbosity = $self->verbose;
	    $self->verbose(1);
	    $self->debug("Skipping import of lane: Cannot import $file\n");
	    $self->verbose($verbosity);
	    $self->{actions} = []; 
	}
    }

    return $self;
}

#---------- get_bams ---------------------

# Requires nothing
sub get_bams_requires
{
    my ($self) = @_;
    return [];
}

# Return empty file bams.done 
# Imported bam files are deleted at the next stage.
sub get_bams_provides
{
    my ($self, $lane_path) = @_;
    
    unless(-e "$lane_path/$$self{prefix}import_bams.done")
    {
	# Test for imported bams then create bams.done file
	my $files_expected = (@{$$self{files}});
	my $files_found = 0;
	foreach my $file (@{$$self{files}})
	{
	    if(-e "$lane_path/$file"){ $files_found++; }
	}
	
	if($files_found == $files_expected )
	{
	    `touch $lane_path/$$self{prefix}import_bams.done`;
	}
    }

    # Return bams.done file
    return ["$lane_path/$$self{prefix}import_bams.done"];
}

#---------- bam_to_fastq ------------------

sub bam_to_fastq_requires 
{
    my ($self,$lane_path) = @_;
    return ["$lane_path/".$$self{prefix}.'import_bams.done'];
}

sub bam_to_fastq_provides {
  my ($self, $lane_path) = @_;
   
  if( $self->is_paired )
  {
    return ["$self->{lane}_1.fastq.gz", "$self->{lane}_2.fastq.gz", "$self->{lane}_1.fastq.gz.fastqcheck", "$self->{lane}_2.fastq.gz.fastqcheck"];
  }
  else
  {
    return ["$self->{lane}_1.fastq.gz", "$self->{lane}_1.fastq.gz.fastqcheck"];
  }
}

sub is_paired {
  my ($self) = @_;
  return $self->{vrlane}->{is_paired};
}

# Adapted from Mapping.pm 
# Converts from bam to fastq then deletes bam files.
sub bam_to_fastq {
    my ($self, $lane_path, $action_lock) = @_;
    
    my ($bam) = @{$$self{files}};

    my $in_bam = $self->{fsu}->catfile($lane_path, $bam);
    my $fastq_base = $self->{lane};
    
    
    my $memory = $self->{memory};
    if (! defined $memory || $memory < 6900) {
        $memory = 6900;
    }
    my $java_mem = int($memory * 0.9);
    my $queue = $memory >= 30000 ? "hugemem" : "normal";
    
    my $fastqs_str ; 
    if( $self->is_paired )
    {
      $fastqs_str  = qq{ (File::Spec->catfile(\$dir, "$self->{lane}_1.fastq"), File::Spec->catfile(\$dir, "$self->{lane}_2.fastq")) };
    }
    else
    {
      $fastqs_str  = qq{ (File::Spec->catfile(\$dir, "$self->{lane}_1.fastq")) };
    }
    
    # Script to be run by LSF to convert bam to fastq
    # bam2fastq does full sanity checking and safe result file creation
    my $script_name = $self->{fsu}->catfile($lane_path, $self->{prefix}."bam2fastq.pl");
    open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
    print $scriptfh qq{
use strict;
use VertRes::Utils::Sam;
use File::Spec;

my \$dir = '$lane_path';
my \@fastqs = $fastqs_str;
# convert to fastq
VertRes::Utils::Sam->new(verbose => 1, quiet => 0, java_memory => $java_mem)->bam2fastq(qq[$in_bam], qq[$fastq_base]);
};
    unless($self->is_paired){
    print $scriptfh qq{
# rename single-ended fastqs
system("mv $self->{lane}.fastq $self->{lane}_1.fastq");
system("mv $self->{lane}.fastq.fastqcheck $self->{lane}_1.fastq.fastqcheck");
};
    }
    print $scriptfh qq{
foreach my \$fastq (\@fastqs) {
    # compress & checksum fastq
    system("md5sum \$fastq > \$fastq.md5; gzip \$fastq; md5sum \$fastq.gz > \$fastq.gz.md5");
}

# delete bam files
unlink("$bam", "$bam.bai", "$bam.bc", "$bam.md5");

foreach my \$fastq (\@fastqs) {
    # rename the fastqcheck file made by bam2fastq
    system("mv \$fastq.fastqcheck \$fastq.gz.fastqcheck");
}

exit;
    };
    close $scriptfh;
    
    my $job_name = $self->{prefix}.'bam2fastq';
    $self->archive_bsub_files($lane_path, $job_name);
    LSF::run($action_lock, $lane_path, $job_name, {bsub_opts => "-q $queue -M${memory}000 -R 'select[mem>$memory] rusage[mem=$memory]'" }, qq{perl -w $script_name});
    
    # we've only submitted to LSF, so it won't have finished; we always return
    # that we didn't complete
    return $self->{No};
}

#---------- update_db ---------------------

# Requires the gzipped fastq and fastqcheck files.
sub update_db_requires
{
    my ($self, $lane_path) = @_;
    
    if( $self->is_paired )
    {
      return ["$self->{lane}_1.fastq.gz", "$self->{lane}_2.fastq.gz", "$self->{lane}_1.fastq.gz.fastqcheck", "$self->{lane}_2.fastq.gz.fastqcheck"];
    }
    else
    {
      return ["$self->{lane}_1.fastq.gz", "$self->{lane}_1.fastq.gz.fastqcheck"];
    }
}

# This subroutine will check existence of the key 'db'. If present, it is assumed
#   that Import should write the stats and status into the VRTrack database. In this
#   case, 0 is returned, meaning that the task must be run. The task will change the
#   QC status from NULL to pending, therefore we will not be called again.
#
#   If the key 'db' is absent, the empty list is returned and the database will not
#   be written.
#
sub update_db_provides
{
    my ($self) = @_;
    if ( exists($$self{db}) ) { return 0; }
    my @provides = ();
    return \@provides;
}

#---------- Debugging and error reporting -----------------

sub format_msg
{
    my ($self,@msg) = @_;
    return '['. scalar gmtime() ."]\t". join('',@msg);
}

sub warn
{
    my ($self,@msg) = @_;
    my $msg = $self->format_msg(@msg);
    if ($self->verbose > 0) 
    {
        print STDERR $msg;
    }
    $self->log($msg);
}

sub debug
{
    # The granularity of verbose messaging does not make much sense
    #   now, because verbose cannot be bigger than 1 (made Base.pm
    #   throw on warn's).
    my ($self,@msg) = @_;
    if ($self->verbose > 0) 
    {
        my $msg = $self->format_msg(@msg);
        print STDERR $msg;
        $self->log($msg);
    }
}

sub throw
{
    my ($self,@msg) = @_;
    my $msg = $self->format_msg(@msg);
    Utils::error($msg);
}

sub log
{
    my ($self,@msg) = @_;

    my $msg = $self->format_msg(@msg);
    my $status  = open(my $fh,'>>',$self->log_file);
    if ( !$status ) 
    {
        print STDERR $msg;
    }
    else 
    { 
        print $fh $msg; 
    }
    if ( $fh ) { close($fh); }
}


1;

