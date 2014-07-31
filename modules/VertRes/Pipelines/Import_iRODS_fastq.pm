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
use base qw(VertRes::Pipelines::Import_iRODS VertRes::Pipelines::Import);

use strict;
use warnings;
use VertRes::LSF;
use VRTrack::VRTrack;
use VRTrack::Lane;
use VRTrack::File;
use VertRes::Utils::FileSystem;
use VertRes::Pipelines::Import;
use VertRes::Pipelines::Import_iRODS;
use Pathogens::Import::ValidateFastqConversion;


our @actions =
(
    # Create the hierarchy path, download and bamcheck the bam files.
    {
        'name'     => 'get_bams',
        'action'   => \&VertRes::Pipelines::Import_iRODS::get_bams,
        'requires' => \&get_bams_requires, 
        'provides' => \&get_bams_provides,
    },

    # Convert to fastq.
    { 
	'name'     => 'bam_to_fastq',
	'action'   => \&bam_to_fastq,
	'requires' => \&bam_to_fastq_requires, 
        'provides' => \&bam_to_fastq_provides 
    },

    # Compress and validate fastq.
    { 
	'name'     => 'compress_and_validate',
	'action'   => \&compress_and_validate,
	'requires' => \&compress_and_validate_requires, 
        'provides' => \&compress_and_validate_provides 
    },

    # If all files downloaded OK, update the VRTrack database.
    {
        'name'     => 'update_db',
        'action'   => \&update_db,
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
    my $self = $class->VertRes::Pipelines::Import_iRODS::new(%$options,'actions'=>\@actions,%args);

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

# Return bam files for import
sub get_bams_provides
{
    my ($self, $lane_path) = @_;
    return $$self{files};
}

#---------- bam_to_fastq ------------------

sub bam_to_fastq_requires 
{
    my ($self,$lane_path) = @_;
#    return ["$lane_path/".$$self{prefix}.'import_bams.done'];
    return $$self{files};
}

sub bam_to_fastq_provides {
  my ($self, $lane_path) = @_;
   
  if( $self->is_paired )
  {
    return ["$self->{lane}_1.fastq", "$self->{lane}_2.fastq", "$self->{lane}_1.fastq.fastqcheck", "$self->{lane}_2.fastq.fastqcheck"];
  }
  else
  {
    return ["$self->{lane}_1.fastq", "$self->{lane}_1.fastq.fastqcheck"];
  }
}

sub is_paired {
  my ($self) = @_;
  return $self->{vrlane}->{is_paired};
}

# Adapted from Mapping.pm 
# Converts from bam to fastq.
sub bam_to_fastq {
    my ($self, $lane_path, $action_lock) = @_;
    
    my ($bam) = @{$$self{files}};

    my $in_bam = $self->{fsu}->catfile($lane_path, $bam);
    my $fastq_base = $self->{lane};
    
    
    my $memory = $self->{memory};
    if (! defined $memory || $memory < 8000) {
        $memory = 8000;
    }
    
    ### We need to check old jobs here to see if it bummed out because of memory and increase the java memory accordingly
    ### if we know whats its going to be increased to we can set java to 90% of it.
    
    my $java_mem = int($memory * 0.95);
    my $queue = $memory >= 30000 ? "hugemem" : "normal";
    my $samtools_sorting_memory = 300000000;
    
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
use VertRes::Wrapper::samtools;
use Bio::Tradis::DetectTags;
use Bio::Tradis::AddTagsToSeq;

my \$dir = '$lane_path';
my \@fastqs = $fastqs_str;

my \$is_tradis =
  Bio::Tradis::DetectTags->new( bamfile => qq[$in_bam] )->tags_present;
if (defined(\$is_tradis) && \$is_tradis == 1) {
	my \$trbam = qq[$in_bam].".tratmp.bam";
	my \$add_tag_obj =
      Bio::Tradis::AddTagsToSeq->new( bamfile => qq[$in_bam], outfile => \$trbam);
	\$add_tag_obj->add_tags_to_seq();
	system("mv \$trbam $in_bam");

	# Remove previous runs of bamcheck without tags - forces rerun
	my \$bc = qq[$in_bam] . ".bc";
	if(-e \$bc){
		system("rm \$bc");
	}
}

# Remove output files from failed runs.
for my \$fastq (\@fastqs)
{
    unlink(\$fastq,\$fastq.'.fastqcheck');
}

# output reads with PF pass only
system("samtools view -F 0x200 $in_bam > pf_pass.sam");

# Reads with PF fail have bases set to N and quality scores to 0.
system("samtools view -f 0x200 $in_bam | awk -F '\\t'  'BEGIN{OFS=\\"\\t\\";} {gsub(/[ACGT]/,\\"N\\",\\\$10) }; {gsub(/./,\\"!\\",\\\$11) };   1' > pf_fail.sam");

#  remove secondary alignments
system("samtools view -H $in_bam | cat - pf_pass.sam pf_fail.sam | samtools view -F 0x100 -b -S - > $in_bam");
unlink("pf_pass.sam");
unlink("pf_fail.sam");

VertRes::Wrapper::samtools->new()->sort(qq[$in_bam], qq[sorted], n => 1, m => $samtools_sorting_memory);
system("mv sorted.bam $in_bam");

VertRes::Utils::Sam->new(verbose => 1, quiet => 0, java_memory => $java_mem )->bam2fastq(qq[$in_bam], qq[$fastq_base]);
};
    unless($self->is_paired){
    print $scriptfh qq{
# rename single-ended fastqs
system("mv $self->{lane}.fastq $self->{lane}_1.fastq");
system("mv $self->{lane}.fastq.fastqcheck $self->{lane}_1.fastq.fastqcheck");
};
    }
    print $scriptfh qq{
exit;
};
    close $scriptfh;

    my $job_name = $self->{prefix}.'bam2fastq';
    $self->archive_bsub_files($lane_path, $job_name);
    VertRes::LSF::run($action_lock, $lane_path, $job_name, {bsub_opts => "-q $queue -M${memory} -R 'select[mem>$memory] rusage[mem=$memory]'", dont_wait=>1 }, qq{perl -w $script_name});
 
    # we've only submitted to LSF, so it won't have finished; we always return
    # that we didn't complete
    return $self->{No};
}


#---------- compress_and_validate ---------

sub compress_and_validate_requires 
{
  my ($self, $lane_path) = @_;
   
  if( $self->is_paired )
  {
    return ["$self->{lane}_1.fastq", "$self->{lane}_2.fastq", "$self->{lane}_1.fastq.fastqcheck", "$self->{lane}_2.fastq.fastqcheck"];
  }
  else
  {
    return ["$self->{lane}_1.fastq", "$self->{lane}_1.fastq.fastqcheck"];
  }
}

sub compress_and_validate_provides {
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

# Compress and validate fastq files
sub compress_and_validate {
    my ($self, $lane_path, $action_lock) = @_;
    
    my $fastq_base = $self->{lane};
    
    my $memory = $self->{memory};
    if (! defined $memory || $memory < 70) {
        $memory = 70;
    }
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
    
    # Script to be run by LSF
    my $script_name = $self->{fsu}->catfile($lane_path, $self->{prefix}."compressfastq.pl");
    open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
    print $scriptfh qq{
use strict;
use File::Spec;
use Pathogens::Import::CompressAndValidate;

my \$irods  = '${$$self{files}}[0]'; 
my \$dir    = '$lane_path';
my \@fastqs = $fastqs_str;

# Compress, validate and checksum fastqs
my \$validator = Pathogens::Import::CompressAndValidate->new( irods_filename => \$irods, fastq_filenames => \\\@fastqs );
\$validator->is_compressed_and_validated() || die("Compress and validate failed\\n");

exit;
};
    close $scriptfh;

    my $job_name = $self->{prefix}.'compressfastq';
    $self->archive_bsub_files($lane_path, $job_name);
    VertRes::LSF::run($action_lock, $lane_path, $job_name, {bsub_opts => "-q $queue -M${memory} -R 'select[mem>$memory] rusage[mem=$memory]'", dont_wait=>1 }, qq{perl -w $script_name});

    
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

# Update Database and clean large files
#
sub update_db
{
    my ($self,$lane_path,$lock_file) = @_;

    # Update database
    $self->VertRes::Pipelines::Import::update_db($lane_path,$lock_file);

    # Remove Large Files
    my @bam_suffix   = ('bam','bam.bai','bam.md5','bam.bc');
    my @fastq_suffix = ('fastq','fastq.fastqcheck');

    my $bam = $self->{files}->[0];
    my $nonhuman = ($bam =~ /_nonhuman.bam$/) ? '_nonhuman':''; # set for nonhuman bams

    for my $suffix (@bam_suffix)
    {
	# Remove bams
	Utils::CMD(qq[rm $lane_path/$$self{lane}$nonhuman.$suffix]);
    }

    for my $suffix (@fastq_suffix)
    {
	# Remove fastqs
	Utils::CMD(qq[rm $lane_path/$$self{lane}_1.$suffix]);
	Utils::CMD(qq[rm $lane_path/$$self{lane}_2.$suffix]) if $self->is_paired;
    }

    return $$self{'Yes'};
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

