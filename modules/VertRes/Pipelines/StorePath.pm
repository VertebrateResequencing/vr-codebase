=head1 NAME

VertRes::Pipelines::StorePath - pipeline for storing paths data on tier2 disks

=head1 SYNOPSIS

# An example config setup:

# Make the pipeline config, e.g.
echo '__Storing__ storeSample.conf' > pipeline.config

# Where storeSample.conf contains:

root => '/path/to/hierarchy/root', # required if paths in fod are relative, otherwise can be anything
module  => 'VertRes::Pipelines::StorePath',
prefix  => '_',
limit => 50,

fod => 'sample_paths.fod',

data => {
    storage_base => 'hashed_samples/vrtrack_dbname',
}

# Notes:

# Paths to be stored must be supplied in a file given by the 'fod' option.

# Paths will be stored on the Vertebrate Resequencing nfs disks in hashed 
# directories. The 'storage_base' option in the data section can be set to 
# separate data at the top level of these disks. In the examples above, we 
# are storing sample release data from a project associated with a VRTrack
# database. If not supplied, 'storage_base' will default to the value 
# 'hashed_paths'.

# the 'limit' option specifies the maximum number of storage jobs to run at
# once, to protect IO being too great. Defaults to 50.

# run the pipeline:
run-pipeline -c pipeline.config

# (and make sure it keeps running by adding that last to a regular cron job)

=head1 DESCRIPTION

A module for storing paths from fast discs where they were built to slower discs
where they can be archived.

=head1 AUTHOR

Shane McCarthy: sm15@sanger.ac.uk

=cut

package VertRes::Pipelines::StorePath;

use strict;
use warnings;
use VertRes::Utils::FileSystem;
use VertRes::Utils::Hierarchy;
use File::Basename;
use File::Spec;
use Cwd qw(abs_path cwd);
use VertRes::LSF;
use Utils;

use base qw(VertRes::Pipeline);

our $actions = [ { name     => 'store_nfs',
                   action   => \&store_nfs,
                   requires => \&store_nfs_requires, 
                   provides => \&store_nfs_provides },
                 { name     => 'cleanup',
                   action   => \&cleanup,
                   requires => \&cleanup_requires, 
                   provides => \&cleanup_provides } ];

our %options = (bsub_opts => '',
                queue => 'normal',
                do_cleanup => 0,
                storage_base => 'hashed_paths');

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Pipelines::StorePath->new();
 Function: Create a new VertRes::Pipelines::StorePath object;
 Returns : VertRes::Pipelines::StoreSample object
 Args    : lane => readgroup id
           optional args as per VertRes::Pipeline

=cut

sub new 
{
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(%options, actions => $actions, @args);
    
    $$self{fsu} = VertRes::Utils::FileSystem->new();
    $$self{hu} = VertRes::Utils::Hierarchy->new();
    
    return $self;
}

=head2 store_nfs_requires

 Title   : store_nfs_requires
 Usage   : my $required_files = $obj->store_nfs_requires('/path/to/lane');
 Function: Find out what files the store_nfs action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub store_nfs_requires 
{
    my ($self, $path) = @_;
    return [];
}

=head2 store_nfs_provides

 Title   : store_nfs_provides
 Usage   : my $provided_files = $obj->store_nfs_provides('/path/to/lane');
 Function: Find out what files the store_nfs action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub store_nfs_provides 
{
    my ($self, $path) = @_;
    return ['.store_nfs_done'];
}

=head2 store_nfs

 Title   : store_nfs
 Usage   : $obj->store_nfs('/path/to/store', 'lock_filename');
 Function: Moves the directory to an nfs disk, replacing original path with
            a symlink to this location. Creates file '.store_nfs_done' if successful.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : path, name of lock file to use

=cut

sub store_nfs 
{
    my ($self, $path, $action_lock) = @_;
    
    my $orig_bsub_opts = $$self{bsub_opts};
    $$self{bsub_opts} = "-q $$self{queue}";
    
    my $done_file = File::Spec->catfile($path, '.store_nfs_done');
    my $job_name = File::Spec->catfile($path, "$self->{prefix}store_nfs");
    
    my $storage_path = $self->storage_path($path);
    $self->delete_bsub_files($path, File::Basename::basename($job_name));
    VertRes::LSF::run($action_lock, $path, $job_name, $self, 
             qq{perl -MVertRes::Pipelines::StorePath -Mstrict -e "VertRes::Pipelines::StorePath->new()->store_path(qq[$path], qq[$storage_path]) || die qq[Failed to store path, $path in $storage_path]; system(qq[touch $done_file]);"});
    $self->debug("Storing '$path' to '$storage_path'\n");
    
    $$self{bsub_opts} = $orig_bsub_opts;
    
    return $$self{No};
}

=head2 cleanup_requires

 Title   : cleanup_requires
 Usage   : my $required_files = $obj->cleanup_requires('/path/to/lane');
 Function: Find out what files the cleanup action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub cleanup_requires 
{
    my ($self, $path) = @_;
    return ['.store_nfs_done'];
}

=head2 cleanup_provides

 Title   : cleanup_provides
 Usage   : my $provided_files = $obj->cleanup_provides('/path/to/lane');
 Function: Find out what files the cleanup action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub cleanup_provides 
{
    my ($self, $path) = @_;
    return [];
}

=head2 cleanup

 Title   : cleanup
 Usage   : $obj->cleanup('/path/to/lane', 'lock_filename');
 Function: Remove pipeline associated files if the do_cleanup config option is set.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub cleanup {
    my ($self, $path, $action_lock) = @_;
    return $$self{Yes} unless $$self{do_cleanup};
    
    my $prefix = $$self{prefix};
    my $base = File::Spec->catfile($path, $prefix);
    
    foreach my $file (qw(log job_status)) {
        unlink($base.$file) if (-e $base.$file);
    }
    
    foreach my $suffix (qw(.o .e .jids .o.previous .e.previous)) {
        unlink($base.'store_nfs'.$suffix) if (-e $base.'store_nfs'.$suffix);
    }
    
    return $$self{Yes};
}

sub is_finished {
    my ($self, $path, $action) = @_;
    
    if ($$action{name} eq 'cleanup') {
        return $self->{No};
    }
    
    return $self->SUPER::is_finished($path, $action);
}

=head2 storage_path

 Title   : storage_path
 Usage   : $obj->storage_path('/path/to/lane');
 Function: Determines the best nfs storage location for the data.
 Returns : path to storage location.
 Args    : path to be stored

=cut

sub storage_path 
{
    my ($self, $path) = @_;
    
    my $storage_base = File::Spec->catfile($$self{storage_base}, $$self{fsu}->hashed_path($path));
    my $storage_path;
    # Check to see if the storage path or temp starage path already exists on one of 
    # the nfs disks
    my @nfs_disks = $$self{hu}->nfs_disks;
    foreach my $disk (@nfs_disks) 
    {
        $storage_path = File::Spec->catfile($disk, $storage_base);
        return $storage_path if (-d $storage_path || -d $storage_path."_store_nfs_temp");
    }
    # Otherwise choose the best available disk for the storage path
    $storage_path = File::Spec->catfile($$self{hu}->nfs_disk, $storage_base);
    $self->throw("Could not determine storage path for $path\n") unless $storage_path;
    return $storage_path;
}


=head2 store_path

 Title   : store_path
 Usage   : $obj->store_path('/path/to/store', '/nfs/storage/path');
 Function: Moves the directory to an nfs disk, replacing original path with
           a symlink to this location.
 Returns : boolean. True if storage successful, false otherwise.
 Args    : path that will be stored, path to storage location on nfs

=cut

sub store_path 
{
    my ($self, $path, $storage_path) = @_;
    
    my $storage_path_temp = $storage_path."_store_nfs_temp";
    
    # Check that everything is ready for the move to happen
    if (-l $path && -d $storage_path) 
    {
        if (abs_path($path) eq abs_path($storage_path)) 
        {
            # Path already stored
            return 1;
        }
        else
        {
            $self->warn("Path '$path' is a symlink that does not point to the storage path '$storage_path'");
            return 0;
        }
    }
    elsif (-d $path && (-d $storage_path || -d $storage_path_temp)) 
    {
        $self->warn("Storage path '$storage_path' or '$storage_path_temp' already exists");
        return 0;
    }
    
    unless (-d $path) 
    {
        $self->warn("Path '$path' wasn't a directory, can't move it to store it on nfs");
        return 0;
    }
    
    # Do the move to temp storage path
    my $cwd = cwd();
    my $moved = $$self{fsu}->move($path, $storage_path_temp);
    
    unless ($moved) 
    {
        $self->warn("Failed to move $path to storage path '$storage_path_temp'");
        return 0;
    }
    
    # Replace original path with a symlink to the storage path
    symlink($storage_path, $path) || $self->throw("Failed to create symlink from $storage_path to $path");
    rename($storage_path_temp, $storage_path) || $self->throw("Could not rename $storage_path_temp to $storage_path");
    chdir($cwd) || $self->throw("Could not change directory to $cwd");
    
    return 1;
}

1;