=head1 NAME

VertRes::Pipelines::StoreLane - pipeline for storing lane-level data on tier2
                                disks

=head1 SYNOPSIS

# make the config files, which specifies the details for connecting to the
# VRTrack g1k-meta database and the data roots:
echo '__VRTrack_Storing__ g1k_storelane.conf' > pipeline.config
# where g1k_storelane.conf contains:
root    => '/abs/path/to/root/data/dir',
module  => 'VertRes::Pipelines::StoreLane',
prefix  => '_',
limit => 50,
db  => {
    database => 'g1k_meta',
    host     => 'mcs4a',
    port     => 3306,
    user     => 'vreseq_rw',
    password => 'xxxxxxx',
}

# by default __VRTrack_Storing__ will pick up lanes that have been both mapped
# and qcd. To override this, set eg:
# vrtrack_processed_flags => { mapped => 1, stored => 0 },
#
# the 'limit' option specifies the maximum number of storage jobs to run at
# once, to protect IO being too great. Defaults to 50.

# run the pipeline:
run-pipeline -c pipeline.config

# (and make sure it keeps running by adding that last to a regular cron job)

=head1 DESCRIPTION

A module for moving lane directories in the Vertebrate Resequencing Informatics 
data hierarchy from fast discs where they were mapped and qc'd to slower discs
where they can be archived.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Pipelines::StoreLane;

use strict;
use warnings;
use VertRes::Utils::FileSystem;
use VRTrack::VRTrack;
use VRTrack::Lane;
use VertRes::LSF;
use Utils;

use base qw(VertRes::Pipeline);

our $actions = [ { name     => 'store_nfs',
                   action   => \&store_nfs,
                   requires => \&store_nfs_requires, 
                   provides => \&store_nfs_provides } ];

our %options = (bsub_opts => '');

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Pipelines::StoreLane->new(lane => '/path/to/lane');
 Function: Create a new VertRes::Pipelines::StoreLane object;
 Returns : VertRes::Pipelines::StoreLane object
 Args    : lane => readgroup id
           optional args as per VertRes::Pipeline

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(%options, actions => $actions, @args);
    
    my $lane = $self->{lane} || $self->throw("lane readgroup not supplied, can't continue");
    
    # if we've been supplied a list of lane paths to work with, instead of
    # getting the lanes from the db, we won't have a vrlane object; make one
    if (! $self->{vrlane}) {
        $self->throw("db option was not supplied in config") unless $self->{db};
        my $vrtrack = VRTrack::VRTrack->new($self->{db}) or $self->throw("Could not connect to the database\n");
        my $vrlane  = VRTrack::Lane->new_by_name($vrtrack, $lane) or $self->throw("No such lane in the DB: [$lane]");
        $self->{vrlane} = $vrlane;
        
        my $files = $vrlane->files();
        foreach my $file (@{$files}) {
            push @{$self->{files}}, $file->hierarchy_name;
        }
    }
    $self->{vrlane} || $self->throw("vrlane object missing");
    
    $self->{fsu} = VertRes::Utils::FileSystem->new;
    
    return $self;
}

=head2 store_nfs_requires

 Title   : update_db_requires
 Usage   : my $required_files = $obj->store_nfs_requires('/path/to/lane');
 Function: Find out what files the store_nfs action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub store_nfs_requires {
    my ($self, $lane_path) = @_;
    return [];
}

=head2 store_nfs_provides

 Title   : store_nfs_provides
 Usage   : my $provided_files = $obj->store_nfs_provides('/path/to/lane');
 Function: Find out what files the store_nfs action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub store_nfs_provides {
    return [];
}

=head2 store_nfs

 Title   : store_nfs
 Usage   : $obj->store_nfs('/path/to/lane', 'lock_filename');
 Function: Moves the lane directory to an nfs disk. Records in the database that
           the lane has been stored.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub store_nfs {
    my ($self, $lane_path, $action_lock) = @_;
    
    my $vrlane = $self->{vrlane};
    return $self->{Yes} if $vrlane->is_processed('stored');
    my $vrtrack = $vrlane->vrtrack;
    my $db = $vrtrack->database_params;
    my $lane_name = $vrlane->name;
	my $umask    = $self->umask_str;
    
    my $host = $db->{host} || $self->throw("db params missing host");
    my $database = $db->{database} || $self->throw("db params missing database");
    my $user = $db->{user} || $self->throw("db params missing database");
    my $password = $db->{password} || $self->throw("db params missing password");
    my $port = $db->{port} || $self->throw("db params missing port");
    
    my $job_name = $self->{prefix}.'store_nfs';
    my $script_name = $self->{fsu}->catfile($lane_path, $job_name.'.pl');
    
    open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
    print $scriptfh qq{
use strict;
use warnings;
use VertRes::Utils::Hierarchy;
use VRTrack::VRTrack;
use VRTrack::Lane;
$umask

my \$vrtrack = VRTrack::VRTrack->new({host => '$host',
                                     port => '$port',
                                     database => '$database',
                                     user => '$user',
                                     password => '$password'});

my \$vrlane = VRTrack::Lane->new_by_name(\$vrtrack, '$lane_name');

my \$hu = VertRes::Utils::Hierarchy->new();
my \$ok = \$hu->store_lane('$lane_path', \$vrlane);

\$ok || die "Failed to store lane\n";

exit;
    };
    close $scriptfh;
    
    $self->archive_bsub_files($lane_path, $job_name);
    
    VertRes::LSF::run($action_lock, $lane_path, $job_name, {bsub_opts => "-M80 -R 'select[mem>80] rusage[mem=80]'" }, qq{perl -w $script_name});
    
    return $self->{No};
}

sub is_finished {
    my ($self, $lane_path, $action) = @_;
    
    if ($action->{name} eq 'store_nfs') {
        my $vrlane = $self->{vrlane};
        return $self->{Yes} if $vrlane->is_processed('stored');
        return $self->{No};
    }
    
    return $self->SUPER::is_finished($lane_path, $action);
}

1;

