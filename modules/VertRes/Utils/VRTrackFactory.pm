=head1 NAME

VertRes::Utils::VRTrackFactory - factory class for getting VRTrack objects

=head1 SYNOPSIS

use VertRes::Utils::VRTrackFactory;

my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => 'mouse',
                                                          mode => 'r');

my @database_names = VertRes::Utils::VRTrackFactory->databases();

=head1 DESCRIPTION

This is largely deprecated in favour of VRTrack::Factory. Internally we just
call VRTrack::Factory methods.

=head1 AUTHOR

Thomas Keane tk2@sanger.ac.uk

=cut

package VertRes::Utils::VRTrackFactory;
use base qw(VertRes::Base);

use strict;
use warnings;
use VRTrack::Factory;

my $FSU_FILE_EXISTS_DB_NAME = $ENV{FSU_FILE_EXISTS_DB_NAME} || 'vrtrack_fsu_file_exists';
my $NFS_DISC_BASENAME = $ENV{NFS_DISC_BASENAME} || '/nfs/vertres';


=head2 new

 Title   : instantiate
 Usage   : my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(
                                                            database => 'mouse',
                                                            mode => 'r');
 Function: Ask the factory to return a VRTrack object to the database specified.
 Returns : VRTrack object
 Args    : database name: a valid VRTrack database,
           mode: either 'r' or 'rw' connection

=cut

sub instantiate {
    return VRTrack::Factory::instantiate(@_);
}

=head2 connection_details

 Title   : connection_details
 Usage   : my %details = VertRes::Utils::VRTrackFactory->connection_details('r');
 Function: Find out what connection details are being used to instantiate().
 Returns : hash with keys: host, port, user, password
 Args    : mode string (r|rw)

=cut

sub connection_details {
    return VRTrack::Factory::connection_details(@_);
}

=head2 databases

 Title   : databases
 Usage   : my @db_names = VertRes::Utils::VRTrackFactory->databases();
 Function: Find out what databases are available to instantiate. Excludes any
           test databases by default.
 Returns : list of strings
 Args    : boolean, which if true will also return test databases and databases
           with old schema versions (default false). Optionally a second
           boolean, which if true will only return databases that have old
           schema versions (the first boolean must be true for this to work)

=cut

sub databases {
    return VRTrack::Factory::databases(@_);
}

=head2 fsu_file_exists_db_name

 Title   : fsu_file_exists_db_name
 Usage   : my @db_names = VertRes::Utils::VRTrackFactory->fsu_file_exists_db_name();
 Function: Name of database where details of files on disk exist
 Returns : name of database

=cut

sub fsu_file_exists_db_name {
    return $FSU_FILE_EXISTS_DB_NAME;
}

=head2 nfs_disc_basename

 Title   : nfs_disc_basename
 Usage   : my $nfs_disk_basename = VertRes::Utils::VRTrackFactory->nfs_disc_basename();
 Function: Base directory name for archived data
 Returns : directory name

=cut

sub nfs_disc_basename {
    return $NFS_DISC_BASENAME;
}

1;
