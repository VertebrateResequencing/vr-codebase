package VRTrack::VRTrack;
# author: jws
=head1 NAME

VRTrack::VRTrack - Sequence Tracking container

=head1 SYNOPSIS
    my $track = VRTrack::VRTrack->new();

    #get arrayref of projects being tracked for traversing hierarchy
    my $projects = $track->projects();

    #also provides accessors for arbitrary objects in hierarchy
    my $lane = $track->get_lane_by_id
    my $lane = $track->get_lane_by_filename

=head1 DESCRIPTION

Retrieves/adds projects in the sequencing tracking database.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;
no warnings 'uninitialized';
use DBI;
use VRTrack::Project;
use VRTrack::Sample;
use VRTrack::Library;
use VRTrack::Lane;
use VRTrack::File;
use VRTrack::Core_obj;

use constant SCHEMA_VERSION => '5';

=head2 new

  Arg [1]    : hashref of database, host, port, user, password connection
               details
  Example    : my $track = VRTrack::VRTrack->new()
  Description: Returns VRTrack object if can connect to database
  Returntype : VRTrack::VRTrack object

=cut

sub new {
    my ($class, $dbparams) = @_;

    my $self = {};
    bless ($self, $class);
    $dbparams->{'port'} ||= 3306;
    my $dbh = DBI->connect( "DBI:mysql:host=$dbparams->{host};".
                            "port=$dbparams->{port};".
                            "database=$dbparams->{database};",
                                      $dbparams->{'user'},
                                      $dbparams->{'password'},
                             {'RaiseError' => 0, 'PrintError'=>0}
                           );
    #my ($dbname) = $dbh->selectrow_array("select DATABASE()");
    #warn $dbname;
    if ($DBI::err){
          warn(sprintf('DB connection failed: %s', $DBI::errstr));
          return undef;
    }

    $self->{_dbh} = $dbh;

    # Check version is OK.
    my $schema_version = $self->schema_version();
    unless ($schema_version == SCHEMA_VERSION){
          warn(sprintf('wrong schema version.  API is %s, DB is %s', SCHEMA_VERSION,$schema_version));
          return undef;
    }

    return $self;
}


=head2 schema_version

  Arg [1]    : None
  Example    : my $schema_version = $project->schema_version();
  Description: Returns database schema_version
  Returntype : int

=cut

sub schema_version {
    my ($self) = @_;

    my $sql = qq[select schema_version from schema_version];
    my $sth = $self->{_dbh}->prepare($sql);

    my $schema_version;
    if ($sth->execute()){
        $schema_version = $sth->fetchrow_array();
    }
    else{
        die(sprintf('Cannot retrieve schema_version: %s', $DBI::errstr));
    }

    return $schema_version;
}


=head2 projects

  Arg [1]    : None
  Example    : my $projects = $track->projects();
  Description: Returns a ref to an array of the project objects that are being tracked
  Returntype : ref to array of VRTrack::Project objects

=cut

sub projects {
    my ($self) = @_;

    unless ($self->{'projects'}){
        my @projects;
        foreach my $id (@{$self->project_ids()}){
            my $obj = VRTrack::Project->new($self->{_dbh},$id);
            push @projects, $obj;
        }
        $self->{'projects'} = \@projects;
    }

    return $self->{'projects'};
}


=head2 project_ids

  Arg [1]    : None
  Example    : my $project_ids = $project->project_ids();
  Description: Returns a ref to an array of the project IDs that are being tracked
  Returntype : ref to array of integer project IDs

=cut

sub project_ids {
    my ($self) = @_;

    unless ($self->{'project_ids'}){
        my $sql = qq[select project_id from project where latest=true];
        my @projects;
        my $sth = $self->{_dbh}->prepare($sql);

        if ($sth->execute()){
            foreach(@{$sth->fetchall_arrayref()}){
                push @projects, $_->[0];
            }
        }
        else{
            die(sprintf('Cannot retrieve projects: %s', $DBI::errstr));
        }

        $self->{'project_ids'} = \@projects;
    }

    return $self->{'project_ids'};
}


=head2 add_project

  Arg [1]    : project name
  Example    : my $newproj = $track->add_project('NOD mouse');
  Description: create a new project, and if successful, return the object
  Returntype : VRTrack::Project object

=cut

sub add_project {
    my ($self, $name) = @_;
    my $dbh = $self->{_dbh};

    # Project names should not be added twice
    my $obj = $self->get_project_by_name($name);
    if ($obj){
        warn "Project $name is already present in the database\n";
        return undef;
    }

    $obj = VRTrack::Project->create($self->{_dbh},$name);
    delete $self->{'project_ids'};
    delete $self->{'projects'};
    return $obj;
}


=head2 get_project_by_name

  Arg [1]    : project name
  Example    : my $project = $track->get_project_by_name('My project');
  Description: retrieve project object by name
  Returntype : VRTrack::Project object

=cut

sub get_project_by_name {
    my ($self, $name) = @_;
    my $obj = VRTrack::Project->new_by_name($self->{_dbh},$name);
    return $obj;
}


=head2 get_project_by_id

  Arg [1]    : project internal id
  Example    : my $project = $track->get_project_by_id(140);
  Description: retrieve project object by internal id
  Returntype : VRTrack::Project object

=cut

sub get_project_by_id {
    my ($self, $id) = @_;
    my $obj = VRTrack::Project->new($self->{_dbh},$id);
    return $obj;
}


=head2 get_project_by_ssid

  Arg [1]    : project sequencescape id
  Example    : my $project = $track->get_project_by_ssid(140);
  Description: retrieve project object by sequencescape id
  Returntype : VRTrack::Project object

=cut

sub get_project_by_ssid {
    my ($self, $id) = @_;
    my $obj = VRTrack::Project->new_by_ssid($self->{_dbh},$id);
    return $obj;
}


=head2 hierarchy_path_of_lane_name

  Arg [1]    : lane name
  Example    : my $lane_hier = $track->hierarchy_path_of_lane_name('2404_1');
  Description: retrieve the hierarchy path for a lane, to the root of the hierarchy.  Does not check the filesystem.
               Returns undef if hierarchy can not be built.
  Returntype : string

=cut

sub hierarchy_path_of_lane_name {
    my ($self, $lane_name) = @_;
    my $hier_path;
    my $lane = VRTrack::Lane->new_by_name($self->{_dbh},$lane_name);
    if ($lane){
        my $lib = VRTrack::Library->new($self->{_dbh},$lane->library_id);
        if ($lib && $lib->seq_tech){
            my $samp = VRTrack::Sample->new($self->{_dbh},$lib->sample_id);
            if ($samp){
                my $proj = VRTrack::Project->new($self->{_dbh},$samp->project_id);
                if ($proj){
                    $hier_path = join '/', ($proj->hierarchy_name,
                                            $samp->hierarchy_name,
                                            $lib->seq_tech->name,
                                            $lib->hierarchy_name,
                                            $lane->hierarchy_name
                                            );
                }
            }
        }
    }

    return $hier_path;
}


=head2 qc_filtered_lane_names

  Arg [1]    : [optional] list of qc_status filters
  Example    : my $all_lanes = $track->qc_filtered_lane_names();
               my $pend_pass_lanes = $track->qc_filtered_lane_names('pending','passed');
  Description: retrieves a (optionally filtered) list of all lane hierarchy names, ordered by project, sample, library names.
               This is a helper function for the qc web interface for speed.
  Returntype : arrayref

=cut

sub qc_filtered_lane_names {
    my ($self,@filter) = @_;
    my $filterclause;
    if (@filter){
        # input validation
        my %allowed = map {$_ => 1} @{VRTrack::Core_obj::list_enum_vals($self,'lane','qc_status')};
        my @goodfilters = grep {$allowed{lc($_)}} @filter;

	$filterclause = 'and lane.qc_status in (';
	$filterclause .= join (",", map {"'$_'"} @goodfilters).')';
    }
    my @lane_names;
    my $sql =qq[select lane.hierarchy_name 
                from latest_project as project,
                    latest_sample as sample,
                    latest_library as library,
                    latest_lane as lane 
                where lane.library_id = library.library_id 
                      and library.sample_id = sample.sample_id 
                      and sample.project_id = project.project_id 
                $filterclause 
                order by project.hierarchy_name, 
                        sample.name, 
                        library.hierarchy_name, 
                        lane.hierarchy_name];
    my $sth = $self->{_dbh}->prepare($sql);

    my $tmpname;
    if ($sth->execute()){
        $sth->bind_columns ( \$tmpname );
        push @lane_names, $tmpname while $sth->fetchrow_arrayref;
    }
    else{
        die(sprintf('Cannot retrieve projects: %s', $DBI::errstr));
    }

    return \@lane_names;
}


=head2 qc_filtered_file_names

  Arg [1]    : [optional] list of qc_status filters
  Example    : my $all_files = $track->qc_filtered_file_names();
               my $files_to_import = $track->qc_filtered_file_names('NULL');

  Description: retrieves a (optionally filtered) list of all file names,
               ordered by project, sample, library, lane names.

               A QC status of NULL for a file indicates that this file has been imported into the database from
               sequence tracking, but is not been imported into a disk hierarchy.
  Returntype : arrayref

=cut

sub qc_filtered_file_names {
    my ($self,@filter) = @_;
    my $filterclause;
    if (@filter){
        # input validation
        my %allowed = map {$_ => 1} @{VRTrack::Core_obj::list_enum_vals($self,'file','qc_status')};
        my @goodfilters = grep {$allowed{lc($_)}} @filter;
        my $get_null = grep {lc($_) eq 'null'} @filter;

        # NULL complicates things, because the SQL syntax is different.  Can't do = NULL or in (NULL)
        if ($get_null){
            if (scalar @goodfilters){
                $filterclause = 'and (file.qc_status is NULL or file.qc_status in (';
                $filterclause .= join (",", map {"'$_'"} @goodfilters).')';
                $filterclause .= ')';
            }
            else {  # only checking for NULL
                $filterclause = 'and file.qc_status is NULL';
            }
        }
        else {
            $filterclause = 'and file.qc_status in (';
            $filterclause .= join (",", map {"'$_'"} @goodfilters).')';
        }
    }
    my @file_names;
    my $sql =qq[select file.name 
                from latest_project as project,
                    latest_sample as sample,
                    latest_library as library,
                    latest_lane as lane,
                    latest_file as file 
                where file.lane_id = lane.lane_id 
                      and lane.library_id = library.library_id 
                      and library.sample_id = sample.sample_id 
                      and sample.project_id = project.project_id 
                $filterclause 
                order by project.hierarchy_name, 
                        sample.name, 
                        library.hierarchy_name, 
                        lane.hierarchy_name, 
                        file.name];
    my $sth = $self->{_dbh}->prepare($sql);

    my $tmpname;
    if ($sth->execute()){
        $sth->bind_columns ( \$tmpname );
        push @file_names, $tmpname while $sth->fetchrow_arrayref;
    }
    else{
        die(sprintf('Cannot retrieve projects: %s', $DBI::errstr));
    }

    return \@file_names;
}
1;
