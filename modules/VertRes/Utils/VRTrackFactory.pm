=head1 NAME

VertRes::Utils::VRTrackFactory - factory class for getting VRTrack objects

=head1 SYNOPSIS

use VertRes::Utils::VRTrackFactory;

my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => 'mouse',
                                                          mode => 'r');

my @database_names = VertRes::Utils::VRTrackFactory->databases();

=head1 DESCRIPTION

A simple factory class that returns VRTrack objects to centralise the database
connection information.

Database host, port, username and password are set by environment variables:
VRTRACK_HOST
VRTRACK_PORT
VRTRACK_RO_USER  (for the 'r' mode read-only capable username)
VRTRACK_RW_USER  (for the 'rw' mode read-write capable username)
VRTRACK_PASSWORD

=head1 AUTHOR

Thomas Keane tk2@sanger.ac.uk

=cut

package VertRes::Utils::VRTrackFactory;
use base qw(VertRes::Base);

use strict;
use warnings;

use DBI;
use VRTrack::VRTrack;

my $HOST = $ENV{VRTRACK_HOST};
my $PORT = $ENV{VRTRACK_PORT} || 3306;
my $READ_USER = $ENV{VRTRACK_RO_USER};
my $WRITE_USER = $ENV{VRTRACK_RW_USER};
my $WRITE_PASS = $ENV{VRTRACK_PASSWORD};


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
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    
    my $database = $self->{database} || $self->throw("A database name must be provided!");
    my $mode = lc($self->{mode}) || $self->throw("A connection mode name must be provided!");
    
    my %details = VertRes::Utils::VRTrackFactory->connection_details($mode);
    $details{database} = $database;
    
    my $vrtrack = VRTrack::VRTrack->new({%details});
    
    return $vrtrack;
}

=head2 connection_details

 Title   : connection_details
 Usage   : my %details = VertRes::Utils::VRTrackFactory->connection_details('r');
 Function: Find out what connection details are being used to instantiate().
 Returns : hash with keys: host, port, user, password
 Args    : mode string (r|rw)

=cut

sub connection_details {
    my ($class, $mode) = @_;
    my $self = $class->SUPER::new(@_);
    
    $mode = lc($mode) || $self->throw("A connection mode name must be provided!");
    
    $self->throw("Invalid connection mode (r or rw valid): $mode\n") unless $mode =~ /^(?:r|rw)$/;
    
    my $user = $mode eq 'rw' ? $WRITE_USER : $READ_USER;
    my $pass = $mode eq 'rw' ? $WRITE_PASS : '';
    
    return (host => $HOST, port => $PORT, user => $user, password => $pass);
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
    my $class = shift;
    my $include_test_and_old_dbs = shift;
    my $only_old = shift;
    my $self = $class->SUPER::new(@_);
    
    my %dbparams = VertRes::Utils::VRTrackFactory->connection_details('r');
    
    my @databases = grep(s/^DBI:mysql://, DBI->data_sources("mysql", \%dbparams));
    
    # we skip information_schema and any test databases
    @databases = grep(!/^information_schema/, @databases);
    unless ($include_test_and_old_dbs) {
        @databases = grep(!/test/, @databases);
    }
    
    # we have to actually check that these databases are vrtrack databases with
    # the correct schema version
    my $schema_version = VRTrack::VRTrack::SCHEMA_VERSION;
    my %expected_tables;
    foreach (VRTrack::VRTrack->schema()) {
        if (/CREATE TABLE `(.+?)`/i || /create view (\S+)/i) {
            $expected_tables{$1} = 1;
        }
    }
    my @vr_dbs;
    DB: foreach my $db (@databases) {
        my %cd = $self->connection_details('r');
        my $dbh = DBI->connect("dbi:mysql:$db;host=$cd{host};port=$cd{port}", $cd{user}, $cd{password}, { RaiseError => 0 });
        unless ($dbh) {
            $self->warn("Could not connect to database $db to check if it was a VRTrack database");
            next;
        }
        
        my %tables = map { s/`//g; s/^$db\.//; $_ => 1 } $dbh->tables();
        
        foreach my $etable (keys %expected_tables) {
            next DB unless exists $tables{$etable};
        }
        foreach my $table (keys %tables) {
            next DB unless exists $expected_tables{$table};
        }
        
        my $sql = qq[ select * from schema_version ];
        my $rows = $dbh->selectall_arrayref($sql);
        my $is_old = $rows->[0]->[0] < $schema_version;
        unless ($include_test_and_old_dbs) {
            next unless $rows->[0]->[0] == $schema_version;
        }
        
        if (! $only_old || $only_old && $is_old) {
            push(@vr_dbs, $db);
        }
    }
    
    return @vr_dbs;
}

1;
