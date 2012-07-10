package VRTrack::VRTrack;
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

==head1 NOTES

A mysql database required. The schema is mysql specific, so other drivers
cannot be used instead.

=head1 AUTHOR

jws@sanger.ac.uk

=head1 METHODS

=cut

# author: jws
use strict;
use warnings;
use Carp qw(confess croak cluck);
no warnings 'uninitialized';
use DBI;
use File::Spec;
use VRTrack::Project;
use VRTrack::Sample;
use VRTrack::Library;
use VRTrack::Seq_request;
use VRTrack::Lane;
use VRTrack::File;
use VRTrack::Core_obj;
use VRTrack::History;

use constant SCHEMA_VERSION => '19';

our $DEFAULT_PORT = 3306;

our %platform_aliases = (ILLUMINA => 'SLX',
                         Illumina => 'SLX',
                         LS454 => '454');

=head2 new

  Arg [1]    : hashref of {database, host, port, user, password}
               connection details. port defaults to 3306.
  Example    : my $track = VRTrack::VRTrack->new()
  Description: Returns VRTrack object if can connect to database
  Returntype : VRTrack::VRTrack object

=cut

sub new {
    my ($class, $dbparams) = @_;

    my $self = {};
    bless ($self, $class);
    $dbparams->{port} ||= $DEFAULT_PORT;
    my $dbh = DBI->connect( "DBI:mysql:host=$dbparams->{host};".
                            "port=$dbparams->{port};".
                            "database=$dbparams->{database};",
                                      $dbparams->{'user'},
                                      $dbparams->{'password'},
                             {'RaiseError' => 0, 'PrintError'=>0}
                           );
    
    if ($DBI::err){
          warn(sprintf('DB connection failed: %s', $DBI::errstr));
          return undef;
    }
    
    $self->{_db_params} = $dbparams;
    $self->{_dbh} = $dbh;
    $self->{transaction} = 0;

    $dbh->{RaiseError} = 1;
    $dbh->{PrintError} = 1;

    # Check version is OK.
    my $schema_version = $self->schema_version();
    unless ($schema_version == SCHEMA_VERSION){
          warn(sprintf('wrong schema version.  API is %s, DB is %s', SCHEMA_VERSION,$schema_version));
          return undef;
    }

    return $self;
}

=head2 schema

  Arg [1]    : n/a
  Example    : foreach (VRTrack::VRTrack->schema()) { print }
  Description: Get an array of sql lines suitable for printing out and streaming
               into your database to (drop and then) create all the VRTrack
               tables. WARNING: using these sql lines on an existing database
               will DESTROY ALL DATA!
  Returntype : array of ;\n termianted sql strings

=cut

sub schema {
    my @sql;
    
    my $line = '';
    while (<DATA>) {
        chomp;
        next if /^--/;
        next unless /\S/;
        $line .= $_;
        if (/;\s*$/) {
            push(@sql, $line."\n");
            $line = '';
        }
    }
    if ($line =~ /;\s*$/) {
        push(@sql, $line);
    }
    
    return @sql;
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

    # removed cache here, otherwise we would have a reference to a bunch
    # of objects that had a reference back to us.  Bad thing.
    my @projects;
    foreach my $id (@{$self->project_ids()}){
        my $obj = VRTrack::Project->new($self,$id);
        push @projects, $obj;
    }
    return \@projects;
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
        my $history_sql = VRTrack::Core_obj->_history_sql;
        my $sql = qq[select distinct(project_id) from project where 1=1 $history_sql];
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

    $obj = VRTrack::Project->create($self,$name);
    delete $self->{'project_ids'};
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
    my $obj = VRTrack::Project->new_by_name($self,$name);
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
    my $obj = VRTrack::Project->new($self,$id);
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
    my $obj = VRTrack::Project->new_by_ssid($self,$id);
    return $obj;
}


=head2 hierarchy_path_of_lane_hname

  Arg [1]    : lane name
  Example    : my $lane_hier = $track->hierarchy_path_of_lane_hname('2404_1');
  Description: retrieve the hierarchy path for a lane, to the root of the hierarchy.  Does not check the filesystem.
               Returns undef if hierarchy can not be built.
  Returntype : string

=cut

sub hierarchy_path_of_lane_hname {
    my ($self, $lane_name) = @_;
    my $lane = VRTrack::Lane->new_by_hierarchy_name($self,$lane_name);
    return $self->hierarchy_path_of_lane($lane);
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
    my $lane = VRTrack::Lane->new_by_name($self,$lane_name);
    if ( !$lane ) { return undef; }
    return $self->hierarchy_path_of_lane($lane);
}

=head2 hierarchy_path_of_lane

  Arg [1]    : VRTrack::Lane object
  Arg [2]    : hierarchy template
  Example    : my $lane_hier = $track->hierarchy_path_of_lane($vrlane);
  Description: Retrieve the hierarchy path for a lane according to the template
               defined by environment variable 'DATA_HIERARCHY', to the root of
               the hierarchy. Template defaults to:
               'project:sample:technology:library:lane'
               Possible terms are 'genus', 'species-subspecies', 'strain',
               'individual', 'project', 'projectid', 'sample', 'technology',
               'library', 'lane'. ('strain' and 'individual' are synonymous)
               Does not check the filesystem.
               Returns undef if hierarchy can not be built.
  Returntype : string

=cut

sub hierarchy_path_of_lane {
    my ($self, $lane, $template) = @_;
    ($lane && ref($lane) && $lane->isa('VRTrack::Lane')) || confess "A VRTrack::Lane must be supplied\n";

    return $self->hierarchy_path_of_object($lane,$template);
}

=head2 hierarchy_path_of_object

  Arg [1]    : VRTrack::Project, VRTrack::Sample, VRTrack::Library or VRTrack::Lane object
  Arg [2]    : hierarchy template (Optional)
  Example    : my $lane_hier = $track->hierarchy_path_of_object($vrlane);
  Description: Retrieve the hierarchy path for a project, sample, library or lane according 
               to the template defined by environment variable 'DATA_HIERARCHY', to the
               root of the hierarchy. Template defaults to:
               'project:sample:technology:library:lane'
               Possible terms are 'genus', 'species-subspecies', 'strain',
               'individual', 'project', 'projectid', 'sample', 'technology',
               'library', 'lane'. ('strain' and 'individual' are synonymous)
               Does not check the filesystem.
               Returns undef if hierarchy cannot be built.
  Returntype : string

=cut

sub hierarchy_path_of_object {
    my ($self, $object, $template) = @_;

    # Object types
    my %object_type = ('VRTrack::Project' => 'project',
		       'VRTrack::Sample'  => 'sample',
		       'VRTrack::Library' => 'library',
		       'VRTrack::Lane'    => 'lane');

    (defined($object) && exists($object_type{ref($object)})) || confess "A recognised object type must be supplied\n";


    # For all acceptable terms, we generate the corresponding word, but for
    # others we just append the term itself to the hierarchy. This allows for
    # words like DATA or TRACKING to be injected into the hierarchy without much
    # fuss. This does however also allow misspelt words and typos (like spcies)
    # to creep into the heirarchy.
    # Since not all possible terms might be used in the template, we will only
    # create objects if necessary (accessing db is expensive).
    my %objs;

    $objs{$object_type{ref($object)}} = $object; # set input object

    my $get_lane = sub { return $objs{lane}; };
    my $get_lib = sub { $objs{library} ||= eval { VRTrack::Library->new($self, &{$get_lane}->library_id) }; return $objs{library}; };
    my $get_sample = sub { $objs{sample} ||= eval { VRTrack::Sample->new($self, &{$get_lib}->sample_id) };  return $objs{sample}; };
    my $get_project = sub { $objs{project} ||= eval { VRTrack::Project->new($self, &{$get_sample}->project_id) }; return $objs{project}; };
    my $get_individual = sub { $objs{individual} ||= eval { VRTrack::Individual->new($self, &{$get_sample}->individual_id) }; return $objs{individual}; };
    my $get_species = sub { $objs{species} ||= eval { VRTrack::Species->new($self, &{$get_individual}->species_id) }; return $objs{species}; };
    my %terms = (genus => $get_species,
                 'species-subspecies' => $get_species,
                 strain => $get_individual,
                 individual => $get_individual,
                 project => $get_project,
                 projectid => $get_project,
                 projectssid => $get_project,  
                 sample => $get_sample,
                 technology => $get_lib,
                 library => $get_lib,
                 lane => $get_lane);
    
    if ( not defined ($template) ) {
        $template = $ENV{DATA_HIERARCHY} || 'project:sample:technology:library:lane';
    }
    my @path = split(/:/, $template);

    my @hier_path_bits;
    foreach my $term (@path) {
        my $get_method = $terms{$term};
        if ($get_method) {
            my $obj = &{$get_method};
            if ($obj) {
                if ($term eq 'genus') {
                    push(@hier_path_bits, $obj->genus); # This is currently the first word of the species name
                }
                elsif ($term eq 'species-subspecies') {
                    my $species_subspecies = $obj->species_subspecies; # For now this is everything after the first space in species name
                    $species_subspecies =~ s/\W/_/g; # replace any non-word char with underscores
                    $species_subspecies =~ s/_+/_/g; # but don't any more than one underscore at a time
                    push(@hier_path_bits, $species_subspecies);
                }
                elsif ($term eq 'projectid') {
                    push(@hier_path_bits, $obj->id);
                }
                elsif ($term eq 'projectssid') {
                    push(@hier_path_bits, $obj->ssid);
                }
                elsif ($term eq 'technology') {
                    push(@hier_path_bits, $obj->seq_tech->name);
                }
                else {
                    push(@hier_path_bits, $obj->hierarchy_name);
                }
            }
            else {
                # the POD says returns undef, but wouldn't it make more sense
                # to warn or die?
                return;
            }
        }
        else {
            push(@hier_path_bits, $term);
        }

	last if $term eq $object_type{ref($object)}; # Finish at object directory
    }
    
    return File::Spec->catdir(@hier_path_bits);
}


=head2 processed_lane_hnames

  Arg [1]    : list of flags and values
  Example    : my $all_lanes   = $track->processed_lane_hnames();
               my $qc_lanes    = $track->processed_lane_hnames('qc'=>1);
               my $no_qc_lanes = $track->processed_lane_hnames('qc'=>0);
  Description: retrieves a (optionally filtered) list of all lane hierarchy names, ordered by project, sample, library names.
               This is a helper function for the qc web interface for speed.
  Returntype : arrayref

=cut

sub processed_lane_hnames {
    my ($self,@filter) = @_;
    if ( scalar @filter % 2 ) { croak "Expected list of keys and values.\n"; }
    my %flags = VRTrack::Core_obj->allowed_processed_flags();
    my @goodfilters;
    my $filterclause = '';
    while (my ($key,$value)=splice(@filter,0,2))
    {
        if ( !exists($flags{$key}) ) { croak qq[The flag "$key" not recognised.\n]; }
        push @goodfilters, ($value ? '' : '!') . qq[(lane.processed & $flags{$key})];
    }
    if ( scalar @goodfilters )
    {
        $filterclause = ' AND ' . join(' AND ', @goodfilters);
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


=head2 qc_filtered_lane_hnames

  Arg [1]    : [optional] list of qc_status filters
  Example    : my $all_lanes = $track->qc_filtered_lane_hnames();
               my $pend_pass_lanes = $track->qc_filtered_lane_hnames('pending','passed');
  Description: retrieves a (optionally filtered) list of all lane hierarchy names, ordered by project, sample, library names.
               This is a helper function for the qc web interface for speed.
  Returntype : arrayref

=cut

sub qc_filtered_lane_hnames {
    my ($self,@filter) = @_;
    return $self->_qc_filtered_lane_field('hierarchy_name',@filter);
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
    return $self->_qc_filtered_lane_field('name',@filter);
}

sub _qc_filtered_lane_field
{
    my ($self,$field,@filter) = @_;
    my $filterclause;
    if (@filter){
        # input validation
        my %allowed = map {$_ => 1} @{VRTrack::Core_obj::list_enum_vals($self,'lane','qc_status')};
        my @goodfilters = grep {$allowed{lc($_)}} @filter;

	$filterclause = 'and lane.qc_status in (';
	$filterclause .= join (",", map {"'$_'"} @goodfilters).')';
    }
    my @lane_names;
    my $sql =qq[select lane.$field
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


=head2 qc_filtered_lib_hnames

  Arg [1]    : [optional] list of qc_status filters
  Example    : my $all_libs = $track->qc_filtered_lib_hnames();
               my $pend_pass_libs = $track->qc_filtered_lib_hnames('pending','passed');
  Description: retrieves a (optionally filtered) list of all lib hierarchy names, ordered by project, sample, name
               This is a helper function for the qc web interface for speed.
  Returntype : arrayref

=cut

sub qc_filtered_lib_hnames {
    my ($self,@filter) = @_;
    my $filterclause;
    if (@filter){
        # input validation
        my %allowed = map {$_ => 1} @{VRTrack::Core_obj::list_enum_vals($self,'library','qc_status')};
        my @goodfilters = grep {$allowed{lc($_)}} @filter;

	$filterclause = 'and library.qc_status in (';
	$filterclause .= join (",", map {"'$_'"} @goodfilters).')';
    }
    my @lib_names;
    my $sql =qq[select library.hierarchy_name 
                from latest_project as project,
                    latest_sample as sample,
                    latest_library as library 
                where library.sample_id = sample.sample_id 
                      and sample.project_id = project.project_id 
                $filterclause 
                order by project.hierarchy_name, 
                        sample.name, 
                        library.hierarchy_name;];
    my $sth = $self->{_dbh}->prepare($sql);
	
    my $tmpname;
    if ($sth->execute()){
        $sth->bind_columns ( \$tmpname );
        push @lib_names, $tmpname while $sth->fetchrow_arrayref;
    }
    else{
        die(sprintf('Cannot retrieve projects: %s', $DBI::errstr));
    }
	
    return \@lib_names;
}

=head2 processed_file_hnames

  Arg [1]    : list of flags and values
  Example    : my $all_files   = $track->processed_file_hnames();
               my $qc_files    = $track->processed_file_hnames('qc'=>1);
               my $no_qc_files = $track->processed_file_hnames('qc'=>0);
  Description: retrieves a (optionally filtered) list of all file hierarchy names, ordered by project, sample, library names.
               This is a helper function for the qc web interface for speed.
  Returntype : arrayref

=cut

sub processed_file_hnames {
    my ($self,@filter) = @_;
    if ( scalar @filter % 2 ) { croak "Expected list of keys and values.\n"; }
    my %flags = VRTrack::Core_obj->allowed_processed_flags();
    my @goodfilters;
    my $filterclause = '';
    while (my ($key,$value)=splice(@filter,0,2))
    {
        if ( !exists($flags{$key}) ) { croak qq[The flag "$key" not recognised.\n]; }
        push @goodfilters, ($value ? '' : '!') . qq[(file.processed & $flags{$key})];
    }
    if ( scalar @goodfilters )
    {
        $filterclause = ' AND ' . join(' AND ', @goodfilters);
    }
    my @file_names;
    my $sql =qq[select file.hierarchy_name 
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
        die(sprintf('Cannot retrieve processed file hnames: %s. The query was %s', $DBI::errstr, $sql));
    }

    return \@file_names;
}


=head2 individual_names 

  Arg [1]    : None
  Example    : my $name_list   = $track->individual_names();
  Description: retrieves a list of all individual names.
               This is a helper function for updating the vrtrack database
               from the warehouse.
  Returntype : arrayref

=cut

sub individual_names {
    my ($self) = @_;
    my @individual_names;
    my $sql =qq[select name
                from individual
                order by name
                ];
    my $sth = $self->{_dbh}->prepare($sql);

    my $tmpname;
    if ($sth->execute()){
        $sth->bind_columns ( \$tmpname );
        push @individual_names, $tmpname while $sth->fetchrow_arrayref;
    }
    else{
        die(sprintf('Cannot retrieve individual names: %s. The query was %s', $DBI::errstr, $sql));
    }

    return \@individual_names;
}

=head2 lane_info

 Title   : lane_info
 Usage   : my %info = $obj->lane_info('lane_name');
 Function: Get information about a lane from the database.
 Returns : hash of information, with keys:
           hierarchy_path => string,
           study          => string, (the true project code)
           project        => string, (may not be the true project code)
           sample         => string,
           individual     => string,
           individual_acc => string,
           individual_coverage => float, (the coverage of this lane's individual)
           population     => string,
           technology     => string, (aka platform, the way DCC puts it, eg.
                                      'ILLUMINA' instead of 'SLX')
           seq_tech       => string, (aka platform, the way Sanger puts it, eg.
                                      'SLX' instead of 'ILLUMINA')
           library        => string, (the hierarchy name, which is most likely
                                      similar to the true original library name)
           library_raw    => string, (the name stored in the database, which may
                                      be a uniquified version of the original
                                      library name)
           library_true   => string, (an attempt at getting the true original
                                      library name, as it was before it was
                                      munged in various ways to create library
                                      and library_raw)
           lane           => string, (aka read group)
           centre         => string, (the sequencing centre name)
           species        => string, (may be undef)
           insert_size    => int, (can be undef if this lane is single-ended)
           withdrawn      => boolean,
           imported       => boolean,
           mapped         => boolean,
           vrlane         => VRTrack::Lane object
           (returns undef if lane name isn't in the database)
 Args    : lane name (read group) OR a VRTrack::Lane object.

           optionally, pre_swap => 1 to get info applicable to the lane in its
           state immediately prior to the last time is_processed('swapped', 1)
           was called on it.

           optionally, get_coverage => 1 to calculate (can be very slow!)
           individual_coverage. To configure this, supply the optional args
           understood by individual_coverage()

=cut

sub lane_info {
    my ($vrtrack, $lane, %args) = @_;
    
    my $hist = VRTrack::History->new();
    my $orig_time_travel = $hist->time_travel;
    
    my ($rg, $vrlane);
    if (ref($lane) && $lane->isa('VRTrack::Lane')) {
        $vrlane = $lane;
        $rg = $vrlane->hierarchy_name;
        $lane = $rg;
    }
    else {
        $vrlane = VRTrack::Lane->new_by_hierarchy_name($vrtrack, $lane);
        $rg = $lane;
    }
    
    return unless ($rg && $vrlane);
    
    my $datetime = 'latest';
    if ($args{pre_swap}) {
        $datetime = $hist->was_processed($vrlane, 'swapped');
    }
    # make sure we've got a lane of the correct time period
    $hist->time_travel($datetime);
    $vrlane = VRTrack::Lane->new_by_hierarchy_name($vrtrack, $lane) || confess("Could not get a vrlane with name $lane prior to $datetime");
    
    my %info = (lane => $rg, vrlane => $vrlane);
    
    $info{hierarchy_path} = $vrtrack->hierarchy_path_of_lane($vrlane);
    $info{withdrawn} = $vrlane->is_withdrawn;
    $info{imported} = $vrlane->is_processed('import');
    $info{mapped} = $vrlane->is_processed('mapped');
    
    my %objs = $vrtrack->lane_hierarchy_objects($vrlane);
    
    $info{insert_size} = $objs{library}->insert_size;
    $info{library} = $objs{library}->hierarchy_name || confess("library hierarchy_name wasn't known for $rg");
    my $lib_name = $objs{library}->name || confess("library name wasn't known for $rg");
    $info{library_raw} = $lib_name;
    ($lib_name) = split(/\|/, $lib_name);
    $info{library_true} = $lib_name;
    $info{centre} = $objs{centre}->name || confess("sequencing centre wasn't known for $rg");
    my $seq_tech = $objs{platform}->name || confess("sequencing platform wasn't known for $rg");
    $info{seq_tech} = $seq_tech;
    if ($seq_tech =~ /illumina|slx/i) {
        $info{technology} = 'ILLUMINA';
    }
    elsif ($seq_tech =~ /solid/i) {
        $info{technology} = 'ABI_SOLID';
    }
    elsif ($seq_tech =~ /454/) {
        $info{technology} = 'LS454';
    }
    $info{sample} = $objs{sample}->name || confess("sample name wasn't known for $rg");
    $info{individual} = $objs{individual}->name || confess("individual name wasn't known for $rg");
    $info{species} =  $objs{species}->name if $objs{species};#|| $self->throw("species name wasn't known for $rg");
    $info{individual_acc} = $objs{individual}->acc; # || $self->throw("sample accession wasn't known for $rg");
    if ($args{get_coverage}) {
        $info{individual_coverage} = $vrtrack->hierarchy_coverage(individual => [$info{individual}],
                                                                  $args{genome_size} ? (genome_size => $args{genome_size}) : (),
                                                                  $args{gt_confirmed} ? (gt_confirmed => $args{gt_confirmed}) : (),
                                                                  $args{qc_passed} ? (qc_passed => $args{qc_passed}) : (),
                                                                  $args{mapped} ? (mapped => $args{mapped}) : ());
    }
    $info{population} = $objs{population}->name;
    $info{project} = $objs{project}->name;
    $info{study} = $objs{study} ? $objs{study}->acc : $info{project};
    
    $hist->time_travel($orig_time_travel);
    
    return %info;
}

=head2 lane_hierarchy_objects

 Title   : lane_hierarchy_objects
 Usage   : my %objects = $obj->lane_hierarchy_objects($lane);
 Function: Get all the parent objects of a lane, from the library up to the
           project.
 Returns : hash with these key and value pairs:
           study => VRTrack::Study object
           project => VRTrack::Project object
           sample => VRTrack::Sample object
           individual => VRTrack::Individual object
           population => VRTrack::Population object
           platform => VRTrack::Seq_tech object
           centre => VRTrack::Seq_centre object
           library => VRTrack::Library object
           species => VRTrack::Species object
 Args    : VRTrack::Lane object

=cut

sub lane_hierarchy_objects {
    my ($vrtrack, $vrlane) = @_;
    
    my $lib = VRTrack::Library->new($vrtrack, $vrlane->library_id);
    my $sc = $lib->seq_centre;
    my $st = $lib->seq_tech;
    my $sample = VRTrack::Sample->new($vrtrack, $lib->sample_id);
    my $individual = $sample->individual;
    my $species = $individual->species;
    my $pop = $individual->population;
    my $project_obj = VRTrack::Project->new($vrtrack, $sample->project_id);
    my $study_obj = VRTrack::Study->new($vrtrack, $project_obj->study_id) if $project_obj->study_id;
    
    return (study => $study_obj,
            project => $project_obj,
            sample => $sample,
            individual => $individual,
            population => $pop,
            platform => $st,
            centre => $sc,
            library => $lib,
            species => $species);
}

=head2 hierarchy_coverage

 Title   : hierarchy_coverage
 Usage   : my $coverage = $obj->hierarchy_coverage(sample => ['NA19239'],
                                                   genome_size => 3e9);
 Function: Discover the sequencing coverage calculated over certain lanes.
 Returns : float
 Args    : At least one hierarchy level as a key, and an array ref of names
           as values, eg. sample => ['NA19239'], platform => ['SLX', '454'].
           Valid key levels are project, sample, individual, population,
           platform, centre and library. (With no options at all, coverage will
           be calculated over all lanes in the database)
           -OR-
           A special mode can be activated by supplying a single lane name with
           the key lane, and a desired level with the level key, eg.:
           lane => 'lane_name', level => 'individual'. This would calculate the
           coverage of all the lanes that belong to the individual that the
           supplied lane belonged to.

           plus optional hash:
           genome_size => int (total genome size in bp; default 3e9)
           gt_confirmed => boolean (only consider genotype confirmed lanes;
                                    default false)
           qc_passed => boolean (only consider qc passed lanes; default false)
           mapped => boolean (coverage of mapped bases; default false: coverage
                              of total bases)

=cut

sub hierarchy_coverage {
    my ($vrtrack, %args) = @_;
    my $genome_size = delete $args{genome_size} || 3e9;
    my $gt = delete $args{gt_confirmed} ? 1 : 0;
    my $qc = delete $args{qc_passed} ? 1 : 0;
    my $mapped = delete $args{mapped} ? 1 : 0;
    
    if (exists $args{lane} || exists $args{level}) {
        my $lane = delete $args{lane};
        my $level = delete $args{level};
        confess("Both lane and level options must be supplied if either of them are") unless $lane && $level;
        
        my @levels = qw(project sample individual population platform centre library);
        foreach my $valid_level (@levels) {
            confess("'$valid_level' option is mutually exclusive of lane&level") if exists $args{$valid_level};
        }
        my %levels = map { $_ => 1 } @levels;
        confess("Supplied level '$level' wasn't valid") unless exists $levels{$level};
        
        my $vrlane = VRTrack::Lane->new_by_name($vrtrack, $lane) || confess("Could not get a lane from the db with name '$lane'");
        
        my %objs = $vrtrack->lane_hierarchy_objects($vrlane);
        confess("Could not get the $level of lane $lane") unless defined $objs{$level};
        
        $args{$level} = [$objs{$level}->name];
    }
    
    my @store = ($gt, $qc, $mapped);
    while (my ($key, $val) = each %args) {
        unless (ref($val)) {
            push(@store, $val);
        }
        else {
            if (ref($val) eq 'ARRAY') {
                push(@store, @{$val});
            }
            elsif (ref($val) eq 'HASH') {
                while (my ($sub_key, $sub_val) = each %{$val}) {
                    push(@store, $sub_val);
                }
            }
        }
    }
    my $store = join(",", sort @store);
    
    unless (defined $vrtrack->{_cover_bases}->{$store}) {
        my @lanes = $vrtrack->get_lanes(%args);
        @lanes || return 0;
        my $bps = 0;
        
        # sum raw bases for all the qc passed, gt confirmed and not withdrawn
        # lanes
        foreach my $lane (@lanes) {
            next if $lane->is_withdrawn;
            if ($gt) {
                next unless ($lane->genotype_status && $lane->genotype_status eq 'confirmed');
            }
            if ($qc) {
                next unless ($lane->qc_status && $lane->qc_status eq 'passed');
            }
            
            my $bp = $lane->raw_bases || 0;
            
            if ($mapped) {
                my $mapstats = $lane->latest_mapping;
                
                if ($mapstats && $mapstats->raw_bases){
                    if ($mapstats->genotype_ratio) {
                        # this is a QC mapped lane, so we make a projection
                        $bps += $bp * ($mapstats->rmdup_bases_mapped / $mapstats->raw_bases);
                    }
                    else {
                        # this is a fully mapped lane, so we know the real answer
                        $bps += $mapstats->bases_mapped;
                    }
                }
                else {
                    $bps += $bp * 0.9; # not sure what else to do here?
                }
            }
            else {
                $bps += $bp;
            }
        }
        
        $vrtrack->{_cover_bases}->{$store} = $bps;
    }
    
    return sprintf('%.2f', $vrtrack->{_cover_bases}->{$store} / $genome_size);
}

=head2 get_lanes

 Title   : get_lanes
 Usage   : my @lanes = $obj->get_lanes(sample => ['NA19239']);
 Function: Get all the lanes under certain parts of the hierarchy, excluding
           withdrawn lanes.
 Returns : list of VRTrack::Lane objects
 Args    : At least one hierarchy level as a key, and an array ref of names
           as values, eg. sample => ['NA19239'], platform => ['SLX', '454'].
           Valid key levels are project, sample, individual, population,
           platform, centre, library and species. (With no options at all, all 
           active lanes in the database will be returned)
           Alternatively to supplying hierarchy level keys and array refs of
           allowed values, you can supply *_regex keys with regex string values
           to select all members of that hierarchy level that match the regex,
           eg. project_regex => 'low_coverage' to limit to projects with
           "low_coverage" in the name. _regex only applies to project, sample
           and library.

           By default it won't return withdrawn lanes; change that:
           return_withdrawn => bool

=cut

sub get_lanes {
    my ($vrtrack, %args) = @_;
    
    my @good_lanes;
    foreach my $project (@{$vrtrack->projects}) {
        my $ok = 1;
        if (defined $args{project}) {
            $ok = 0;
            foreach my $name (@{$args{project}}) {
                if ($name eq $project->name || $name eq $project->hierarchy_name || ($project->study && $name eq $project->study->acc)) {
                    $ok = 1;
                    last;
                }
            }
        }
        $ok || next;
        if (defined $args{project_regex}) {
            $project->name =~ /$args{project_regex}/ || next;
        }
        
        foreach my $sample (@{$project->samples}) {
            my $ok = 1;
            if (defined ($args{sample})) {
                $ok = 0;
                foreach my $name (@{$args{sample}}) {
                    if ($name eq $sample->name) {
                        $ok = 1;
                        last;
                    }
                }
            }
            $ok || next;
            if (defined $args{sample_regex}) {
                $sample->name =~ /$args{sample_regex}/ || next;
            }
            
            my %objs;
            $objs{individual} = $sample->individual;
            $objs{population} = $objs{individual}->population;
            $objs{species}    = $objs{individual}->species;
            
            my ($oks, $limits) = (0, 0);
            foreach my $limit (qw(individual population species)) {
                if (defined $args{$limit}) {
                    $limits++;
                    if ($limit eq 'species' && !(defined $objs{'species'})) {
                        confess('species not defined for sample '.$sample->name);
                        last;
                    }
                    my $ok = 0;
                    foreach my $name (@{$args{$limit}}) {
                        if ($name eq $objs{$limit}->name || ($objs{$limit}->can('hierarchy_name') && $name eq $objs{$limit}->hierarchy_name)) {
                            $ok = 1;
                            last;
                        }
                    }
                    $oks += $ok;
                }
            }
            next unless $oks == $limits;
            
            foreach my $library (@{$sample->libraries}) {
                my $ok = 1;
                if (defined ($args{library})) {
                    $ok = 0;
                    foreach my $name (@{$args{library}}) {
                        if ($name eq $library->name || $name eq $library->hierarchy_name) {
                            $ok = 1;
                            last;
                        }
                    }
                }
                $ok || next;
                if (defined $args{library_regex}) {
                    $library->name =~ /$args{library_regex}/ || next;
                }
                
                my %objs;
                $objs{centre} = $library->seq_centre;
                $objs{platform} = $library->seq_tech;
                
                my ($oks, $limits) = (0, 0);
                foreach my $limit (qw(centre platform)) {
                    if (defined $args{$limit}) {
                        $limits++;
                        my $ok = 0;
                        foreach my $name (@{$args{$limit}}) {
                            if ($name eq $objs{$limit}->name) {
                                $ok = 1;
                                last;
                            }
                        }
                        $oks += $ok;
                    }
                }
                next unless $oks == $limits;
                
                push(@good_lanes, @{$library->lanes});
            }
        }
    }
    
    if ($args{return_withdrawn}) {
        return @good_lanes;
    }
    else {
        # filter out withdrawn lanes
        my @active;
        foreach my $lane (@good_lanes) {
            next if $lane->is_withdrawn;
            push(@active, $lane);
        }
        return @active;
    }
}

=head2 transaction

  Arg [1]    : Code ref
  Arg [2]    : optional hash ref: { read => [], write => []} where the array
               ref values contain table names to lock for reading/writing
  Example    : my $worked = $vrtrack->transaction(sub { $lane->update; }, { write => ['lane'] });
  Description: Run code safely in a transaction, with automatic retries in the
               case of deadlocks. If the transaction fails for some other reason,
               the error message can be found in $vrtrack->{transaction_error}.
  Returntype : Boolean

=cut

sub transaction {
    my ($self, $code, $locks) = @_;
    
    my $dbh = $self->{_dbh};
    
    # we wanted to use begin_work() to handle turning off and on AutoCommit, but
    # to be compatible with transaction_start() et al. we'll have to use the
    # same mechanisms
    my $autocommit = $dbh->{AutoCommit};
    $dbh->{AutoCommit} = 0;
    $self->{transaction}++;
    
    # Raise Errors if there are any problems, which we will catch in evals
    my $raiseerror = $dbh->{RaiseError};
    $dbh->{RaiseError} = 1;
    
    #eval {
    #    $dbh->begin_work;
    #};
    #if ($@) {
    #    my $err = $@;
    #    if ($err =~ /Already in a transaction/) {
    #        # we're already in a transaction, so just run the $code without
    #        # doing a commit
    #        $dbh->{RaiseError} = 1;
    #        $dbh->{AutoCommit} = 0;
    #        &$code;
    #        $self->{transaction}--;
    #        return 1;
    #    }
    #    else {
    #        confess $err;
    #    }
    #}
    if ($self->{transaction} > 1) {
        &$code;
        $self->{transaction}--;
        return 1;
    } 
    
    # turn off warnings that may be generated when we call &$code
    my $sig_warn = $SIG{'__WARN__'};
    $SIG{'__WARN__'} = sub { };
    
    # try to run the $code and commit, repeating if deadlock found, die with
    # stack trace for other issues
    my $success = 0;
    delete $self->{transaction_error};
    while (1) {
        eval {
            # make extra sure RaiseError and AutoCommit are set as necessary,
            # incase something nested changes these
            $dbh->{RaiseError} = 1;
            $dbh->{AutoCommit} = 0;
            &$code;
            $dbh->{RaiseError} = 1;
            $dbh->{AutoCommit} = 0;
            $dbh->commit;
        };
        if ($@) {
            my $err = $@;
            eval { $dbh->rollback };
            if ($err =~ /Deadlock found/) {
                sleep(2);
            }
            else {
                chomp($err);
                $self->{transaction_error} = "Transaction failed, rolled back. Error was: $err";
                last;
            }
        }
        else {
            $success = 1;
            last;
        }
    }
    
    $dbh->{AutoCommit} = $autocommit;
    $dbh->{RaiseError} = $raiseerror;
    $SIG{'__WARN__'} = $sig_warn;
    $self->{transaction}--;
    
    return $success;
}

=head2 transaction_start

  Arg [1]    : None
  Example    : $vrtrack->transaction_start();
  Description: 
  Returntype : none

=cut

sub transaction_start {
    my ($self) = @_;

    my $dbh = $self->{_dbh};

    $self->{transaction}++;                         # Increase the counter
    if ( $self->{transaction}>1 ) { return; }       # If already inside a transaction, we are done.

    $self->{_AutoCommit} = $dbh->{AutoCommit};      # Remember the previous state
    $dbh->{AutoCommit} = 0;                         # Start the transaction

    return;
}


=head2 transaction_commit

  Arg [1]    : None
  Example    : $vrtrack->transaction_commit();
  Description: 
  Returntype : none

=cut

sub transaction_commit {
    my ($self) = @_;

    die "transaction_commit: no active transaction\n" unless $self->{transaction}>0;

    $self->{transaction}--;
    if ( $self->{transaction} ) { return; }         # If inside a nested transactions, don't commit yet.

    $self->{_dbh}->commit;
    $self->{_dbh}->{AutoCommit} = $self->{_AutoCommit};

    return;
}


=head2 transaction_rollback

  Arg [1]    : None
  Example    : $vrtrack->transaction_rollback();
  Description: 
  Returntype : none

=cut

sub transaction_rollback {
    my ($self) = @_;

    die "transaction_commit: no active transaction\n" unless $self->{transaction}>0;

    $self->{transaction}--;

    # roll back within eval to prevent rollback
    # failure from terminating the script
    eval { $self->{_dbh}->rollback; };

    $self->{_dbh}->{AutoCommit} = $self->{_AutoCommit};

    if ( $self->{transaction} ) { die "Transaction failed\n"; }     # If inside a nested transaction, return the control higher

    return;
}


=head2 vrtrack

  Arg [1]    : None
  Example    : my $vrtrack = $obj->vrtrack();
  Description: Get vrtrack.  This is self, as we _are_ a VRTrack object.  This call is just to provide consistency to getting new objects through the api by my $sub = VRTrack::Whatever($parent->vrtrack, $id);
  Returntype : VRTrack::VRTrack

=cut

sub vrtrack {
    my ($self) = @_;
    return $self;
}

=head2 database_params

  Arg [1]    : None
  Example    : my $parms = $obj->database_params();
  Description: Get the database parameters that were supplied to new() to create
               this instance of VRTrack::VRTrack. 
  Returntype : hash ref

=cut

sub database_params {
    my $self = shift;
    return $self->{_db_params};
}


=head2 assembly_names

  Arg [1]    : None
  Example    : my @assemblies = $vrtrack->assembly_names();
  Description: Returns a reference to an array of the names in the assembly table sorted by name
  Returntype : reference to array of strings

=cut

sub assembly_names {
    my $self = shift;
    return $self->_list_names('assembly');
}

=head2 species_names

  Arg [1]    : None
  Example    : my @species = $vrtrack->species_names();
  Description: Returns a reference to an array of the names in the species table sorted by name
  Returntype : reference to array of strings

=cut

sub species_names {
    my $self = shift;
    return $self->_list_names('species');
}

=head2 mapper_names

  Arg [1]    : None
  Example    : my @mappers = $vrtrack->mapper_names();
  Description: Returns a reference to an array of the distinct names in the species table
  Returntype : reference to array of strings

=cut

sub mapper_names {
    my $self = shift;
    return $self->_list_names('mapper');
}

# Returns reference to an array of names contained in the name column of a table
# Tables are restricted to tables contained in hash %permitted_table within the function.
sub _list_names {
    my ($self,$table) = @_;

    # List of permitted tables
    my %permitted_table = (assembly => 0, species => 0, mapper => 1);

    unless(exists $permitted_table{$table}){croak qq[The listing names for '$table' not permitted];}
    my $distinct = $permitted_table{$table} ? 'distinct' : '';

    my @names;
    my $sql =qq[select $distinct $table.name 
                from $table
                order by $table.name;];

    my $sth = $self->{_dbh}->prepare($sql);
	
    my $tmpname;
    if ($sth->execute()){
        $sth->bind_columns ( \$tmpname );
        push @names, $tmpname while $sth->fetchrow_arrayref;
    }
    else{
        die(sprintf('Cannot retrieve $table names: %s', $DBI::errstr));
    }

    return \@names;
}


1;


__DATA__
--
-- Table structure for table `version`
--

DROP TABLE IF EXISTS `schema_version`;
CREATE TABLE `schema_version` (
  `schema_version` mediumint(8) unsigned NOT NULL,
  PRIMARY KEY  (`schema_version`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

insert into schema_version(schema_version) values (19);

--
-- Table structure for table `assembly`
--

DROP TABLE IF EXISTS `assembly`;
CREATE TABLE `assembly` (
  `assembly_id` smallint(5) unsigned NOT NULL auto_increment,
  `name` varchar(255) NOT NULL DEFAULT '',
  `reference_size` integer DEFAULT NULL,
  `taxon_id` mediumint(8) unsigned DEFAULT NULL,
  `translation_table` smallint(5) unsigned DEFAULT NULL,
  PRIMARY KEY  (`assembly_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

--
-- Table structure for table `exome_design`
--

DROP TABLE IF EXISTS `exome_design`;
CREATE TABLE `exome_design` (
  `exome_design_id` smallint(5) unsigned NOT NULL auto_increment,
  `name` varchar(255) NOT NULL DEFAULT '',
  `bait_bases` bigint(20) unsigned DEFAULT NULL,
  `target_bases` bigint(20) unsigned DEFAULT NULL,
  PRIMARY KEY (`exome_design_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

--
-- Table structure for table `note`
--

DROP TABLE IF EXISTS `note`;
CREATE TABLE `note` (
  `note_id` mediumint(8) unsigned NOT NULL auto_increment,
  `note` text DEFAULT NULL,
  PRIMARY KEY  (`note_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;


--
-- Table structure for table `file`
--

DROP TABLE IF EXISTS `file`;
CREATE TABLE `file` (
  `row_id` int unsigned NOT NULL auto_increment key,
  `file_id` mediumint(8) unsigned NOT NULL DEFAULT 0,
  `lane_id` mediumint(8) unsigned NOT NULL DEFAULT 0,
  `name` varchar(255) NOT NULL DEFAULT '',
  `hierarchy_name` varchar(255) DEFAULT NULL,
  `processed` int(10) DEFAULT 0,
  `type` tinyint(4) DEFAULT NULL,
  `readlen` smallint(5) unsigned DEFAULT NULL,
  `raw_reads` bigint(20) unsigned DEFAULT NULL,
  `raw_bases` bigint(20) unsigned DEFAULT NULL,
  `mean_q` float unsigned DEFAULT NULL,
  `md5` varchar(40) DEFAULT NULL,
  `note_id` mediumint(8) unsigned DEFAULT NULL,
  `changed` datetime NOT NULL DEFAULT '0000-00-00',
  `latest` tinyint(1) DEFAULT '0',
  KEY `file_id` (`file_id`),
  KEY `lane_id` (`lane_id`),
  KEY `hierarchy_name` (`hierarchy_name`),
  KEY `name` (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

--
-- Table structure for table `image`
--

DROP TABLE IF EXISTS `image`;
CREATE TABLE `image` (
  `image_id` mediumint(8) unsigned NOT NULL auto_increment,
  `mapstats_id` mediumint(8) unsigned NOT NULL DEFAULT 0,
  `name` varchar(40) NOT NULL DEFAULT '',
  `caption` varchar(40) DEFAULT NULL,
  `image` MEDIUMBLOB,
  PRIMARY KEY (`image_id`),
  KEY  `mapstats_id` (`mapstats_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

--
-- Table structure for table `lane`
--

DROP TABLE IF EXISTS `lane`;
CREATE TABLE `lane` (
  `row_id` int unsigned NOT NULL auto_increment key,
  `lane_id` mediumint(8) unsigned NOT NULL DEFAULT 0,
  `library_id` smallint(5) unsigned NOT NULL DEFAULT 0,
  `seq_request_id` mediumint(8) unsigned NOT NULL DEFAULT 0,
  `name` varchar(255) NOT NULL DEFAULT '',
  `hierarchy_name` varchar(255) NOT NULL DEFAULT '',
  `acc` varchar(40) DEFAULT NULL,
  `readlen` smallint(5) unsigned DEFAULT NULL,
  `paired` tinyint(1) DEFAULT NULL,
  `raw_reads` bigint(20) unsigned DEFAULT NULL,
  `raw_bases` bigint(20) unsigned DEFAULT NULL,
  `npg_qc_status` enum('pending','pass','fail','-') DEFAULT 'pending',
  `processed` int(10) DEFAULT 0,
  `auto_qc_status` enum('no_qc','passed','failed') DEFAULT 'no_qc',
  `qc_status` enum('no_qc','pending','passed','failed','gt_pending','investigate') DEFAULT 'no_qc',
  `gt_status` enum('unchecked','confirmed','wrong','unconfirmed','candidate','unknown','swapped') DEFAULT 'unchecked',
  `submission_id` smallint(5) unsigned DEFAULT NULL,
  `withdrawn` tinyint(1) DEFAULT NULL,
  `note_id` mediumint(8) unsigned DEFAULT NULL,
  `changed` datetime NOT NULL DEFAULT '0000-00-00',
  `run_date` datetime DEFAULT NULL,
  `storage_path` varchar(255) DEFAULT NULL,
  `latest` tinyint(1) DEFAULT '0',
  KEY `lane_id` (`lane_id`),
  KEY `lanename` (`name`),
  KEY `library_id` (`library_id`),
  KEY `hierarchy_name` (`hierarchy_name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

--
-- Table structure for table `library`
--

DROP TABLE IF EXISTS `library`;
CREATE TABLE `library` (
  `row_id` int unsigned NOT NULL auto_increment key,
  `library_id` smallint(5) unsigned NOT NULL DEFAULT 0,
  `library_request_id` mediumint(8) unsigned NOT NULL DEFAULT 0,
  `sample_id` smallint(5) unsigned NOT NULL DEFAULT 0,
  `ssid` mediumint(8) unsigned DEFAULT NULL,
  `name` varchar(255) NOT NULL DEFAULT '',
  `hierarchy_name` varchar(255) NOT NULL DEFAULT '',
  `prep_status` enum('unknown','pending','started','passed','failed','cancelled','hold') DEFAULT 'unknown',
  `auto_qc_status` enum('no_qc','passed','failed') DEFAULT 'no_qc',
  `qc_status` enum('no_qc','pending','passed','failed') DEFAULT 'no_qc',
  `fragment_size_from` mediumint(8) unsigned DEFAULT NULL,
  `fragment_size_to` mediumint(8) unsigned DEFAULT NULL,
  `library_type_id` smallint(5) unsigned DEFAULT NULL,
  `library_tag` smallint(5) unsigned,
  `library_tag_group` smallint(5) unsigned,
  `library_tag_sequence` varchar(1024),
  `seq_centre_id` smallint(5) unsigned DEFAULT NULL,
  `seq_tech_id` smallint(5) unsigned DEFAULT NULL,
  `open` tinyint(1) DEFAULT '1',
  `note_id` mediumint(8) unsigned DEFAULT NULL,
  `changed` datetime NOT NULL DEFAULT '0000-00-00',
  `latest` tinyint(1) DEFAULT '0',
  KEY `ssid` (`ssid`),
  KEY `name` (`name`),
  KEY `hierarchy_name` (`hierarchy_name`),
  KEY `sample_id` (`sample_id`),
  KEY `library_id` (`library_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

--
-- Table structure for table `multiplex_pool`
--
DROP TABLE IF EXISTS `multiplex_pool`;
CREATE TABLE `multiplex_pool` (
  `multiplex_pool_id` mediumint(8) unsigned NOT NULL auto_increment key,
  `ssid` mediumint(8) unsigned DEFAULT NULL,
  `name` varchar(255) NOT NULL DEFAULT '',
  `note_id` mediumint(8) unsigned DEFAULT NULL,
  KEY `multiplex_pool_id` (`multiplex_pool_id`),
  KEY `name` (`name`),
  UNIQUE KEY `ssid` (`ssid`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

--
-- Table structure for table `library_multiplex_pool`
--
DROP TABLE IF EXISTS `library_multiplex_pool`;
CREATE TABLE `library_multiplex_pool` (
  `library_multiplex_pool_id` mediumint(8) unsigned NOT NULL auto_increment key,
  `multiplex_pool_id` smallint(5) unsigned NOT NULL DEFAULT 0,
  `library_id` smallint(5) unsigned NOT NULL DEFAULT 0,
  KEY `library_multiplex_pool_id` (`library_multiplex_pool_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;



--
-- Table structure for table `library_request`
--
DROP TABLE IF EXISTS `library_request`;
CREATE TABLE `library_request` (
  `row_id` int unsigned NOT NULL auto_increment key,
  `library_request_id` mediumint(8) unsigned NOT NULL DEFAULT 0,
  `sample_id` smallint(5) unsigned NOT NULL DEFAULT 0,
  `ssid` mediumint(8) unsigned DEFAULT NULL,
  `prep_status` enum('unknown','pending','started','passed','failed','cancelled','hold') DEFAULT 'unknown',
  `note_id` mediumint(8) unsigned DEFAULT NULL,
  `changed` datetime NOT NULL DEFAULT '0000-00-00',
  `latest` tinyint(1) DEFAULT '0',
  KEY `library_request_id` (`library_request_id`),
  KEY `ssid` (`ssid`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

--
-- Table structure for table `seq_request`
--
DROP TABLE IF EXISTS `seq_request`;
CREATE TABLE `seq_request` (
  `row_id` int unsigned NOT NULL auto_increment key,
  `seq_request_id` mediumint(8) unsigned NOT NULL DEFAULT 0,
  `library_id` smallint(5) unsigned,
  `multiplex_pool_id` smallint(5) unsigned,
  `ssid` mediumint(8) unsigned DEFAULT NULL,
  `seq_type` enum('Single ended sequencing','Paired end sequencing','HiSeq Paired end sequencing','MiSeq sequencing','Single ended hi seq sequencing') DEFAULT 'Single ended sequencing',
  `seq_status` enum('unknown','pending','started','passed','failed','cancelled','hold') DEFAULT 'unknown',
  `note_id` mediumint(8) unsigned DEFAULT NULL,
  `changed` datetime NOT NULL DEFAULT '0000-00-00',
  `latest` tinyint(1) DEFAULT '0',
  KEY `seq_request_id` (`seq_request_id`),
  KEY `ssid` (`ssid`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;


--
-- Table structure for table `library_type`
--

DROP TABLE IF EXISTS `library_type`;
CREATE TABLE `library_type` (
  `library_type_id` smallint(5) unsigned NOT NULL auto_increment,
  `name` varchar(40) NOT NULL DEFAULT '',
  PRIMARY KEY  (`library_type_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

--
-- Table structure for table `mapper`
--

DROP TABLE IF EXISTS `mapper`;
CREATE TABLE `mapper` (
  `mapper_id` smallint(5) unsigned NOT NULL auto_increment,
  `name` varchar(40) NOT NULL DEFAULT '',
  `version` varchar(40) NOT NULL DEFAULT 0,
  PRIMARY KEY  (`mapper_id`),
  UNIQUE KEY `name_v` (`name`, `version`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

--
-- Table structure for table `mapstats`
--

DROP TABLE IF EXISTS `mapstats`;
CREATE TABLE `mapstats` (
  `row_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `mapstats_id` mediumint(8) unsigned NOT NULL DEFAULT 0,
  `lane_id` mediumint(8) unsigned NOT NULL DEFAULT 0,
  `mapper_id` smallint(5) unsigned DEFAULT NULL,
  `assembly_id` smallint(5) unsigned DEFAULT NULL,
  `raw_reads` bigint(20) unsigned DEFAULT NULL,
  `raw_bases` bigint(20) unsigned DEFAULT NULL,
  `clip_bases` bigint(20) unsigned DEFAULT NULL,
  `reads_mapped` bigint(20) unsigned DEFAULT NULL,
  `reads_paired` bigint(20) unsigned DEFAULT NULL,
  `bases_mapped` bigint(20) unsigned DEFAULT NULL,
  `rmdup_reads_mapped` bigint(20) unsigned DEFAULT NULL,
  `rmdup_bases_mapped` bigint(20) unsigned DEFAULT NULL,
  `adapter_reads` bigint(20) unsigned DEFAULT NULL,
  `error_rate` float unsigned DEFAULT NULL,
  `mean_insert` float unsigned DEFAULT NULL,
  `sd_insert` float unsigned DEFAULT NULL,
  `gt_expected` varchar(40) DEFAULT NULL,
  `gt_found` varchar(40) DEFAULT NULL,
  `gt_ratio` float unsigned DEFAULT NULL,
  `note_id` mediumint(8) unsigned DEFAULT NULL,
  `changed` datetime NOT NULL DEFAULT '0000-00-00',
  `latest` tinyint(1) DEFAULT '0',
  `bait_near_bases_mapped` bigint(20) unsigned DEFAULT NULL,
  `target_near_bases_mapped` bigint(20) unsigned DEFAULT NULL,
  `bait_bases_mapped` bigint(20) unsigned DEFAULT NULL,
  `mean_bait_coverage` float unsigned DEFAULT NULL,
  `bait_coverage_sd` float unsigned DEFAULT NULL,
  `off_bait_bases` bigint(20) unsigned DEFAULT NULL,
  `reads_on_bait` bigint(20) unsigned DEFAULT NULL,
  `reads_on_bait_near` bigint(20) unsigned DEFAULT NULL,
  `reads_on_target` bigint(20) unsigned DEFAULT NULL,
  `reads_on_target_near` bigint(20) unsigned DEFAULT NULL,
  `target_bases_mapped` bigint(20) unsigned DEFAULT NULL,
  `mean_target_coverage` float unsigned DEFAULT NULL,
  `target_coverage_sd` float unsigned DEFAULT NULL,
  `target_bases_1X` float unsigned DEFAULT NULL,
  `target_bases_2X` float unsigned DEFAULT NULL,
  `target_bases_5X` float unsigned DEFAULT NULL,
  `target_bases_10X` float unsigned DEFAULT NULL,
  `target_bases_20X` float unsigned DEFAULT NULL,
  `target_bases_50X` float unsigned DEFAULT NULL,
  `target_bases_100X` float unsigned DEFAULT NULL,
  `exome_design_id` smallint(5) unsigned DEFAULT NULL,
  `percentage_reads_with_transposon` float unsigned DEFAULT NULL,
  `is_qc` tinyint(1) DEFAULT '0',
  `prefix` varchar(40) DEFAULT '_',
  PRIMARY KEY (`row_id`),
  KEY `mapstats_id` (`mapstats_id`),
  KEY `lane_id` (`lane_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

--
-- Table structure for table `population`
--

DROP TABLE IF EXISTS `population`;
CREATE TABLE `population` (
  `population_id` smallint(5) unsigned NOT NULL auto_increment,
  `name` varchar(40) NOT NULL DEFAULT '',
  PRIMARY KEY  (`population_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

--
-- Table structure for table `species`
--

DROP TABLE IF EXISTS `species`;
CREATE TABLE `species` (
  `species_id` smallint(5) unsigned NOT NULL auto_increment,
  `name` varchar(255) NOT NULL,
  `taxon_id` mediumint(8) unsigned NOT NULL,
  PRIMARY KEY  (`species_id`),
  KEY `name` (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

--
-- Table structure for table `individual`
--

DROP TABLE IF EXISTS `individual`;
CREATE TABLE `individual` (
  `individual_id` smallint(5) unsigned NOT NULL auto_increment,
  `name` varchar(255) NOT NULL DEFAULT '',
  `hierarchy_name` varchar(255) NOT NULL DEFAULT '',
  `alias` varchar(40) NOT NULL DEFAULT '',
  `sex` enum('M','F','unknown') DEFAULT 'unknown',
  `acc` varchar(40) DEFAULT NULL,
  `species_id` smallint(5) unsigned DEFAULT NULL,
  `population_id` smallint(5) unsigned DEFAULT NULL,
  PRIMARY KEY  (`individual_id`),
  UNIQUE KEY `name` (`name`),
  UNIQUE KEY `hierarchy_name` (`hierarchy_name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

--
-- Table structure for table `project`
--

DROP TABLE IF EXISTS `project`;
CREATE TABLE `project` (
  `row_id` int unsigned NOT NULL auto_increment key,
  `project_id` smallint(5) unsigned NOT NULL DEFAULT 0,
  `ssid` mediumint(8) unsigned DEFAULT NULL,
  `name` varchar(255) NOT NULL DEFAULT '',
  `hierarchy_name` varchar(255) NOT NULL DEFAULT '',
  `study_id` smallint(5) DEFAULT NULL,
  `note_id` mediumint(8) unsigned DEFAULT NULL,
  `changed` datetime NOT NULL DEFAULT '0000-00-00',
  `latest` tinyint(1) DEFAULT '0',
  KEY `project_id` (`project_id`),
  KEY `ssid` (`ssid`),
  KEY `latest` (`latest`),
  KEY `hierarchy_name` (`hierarchy_name`),
  KEY `name` (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;


--
-- Table structure for table `study`
--

DROP TABLE IF EXISTS `study`;
CREATE TABLE `study` (
`study_id` smallint(5) unsigned NOT NULL auto_increment,
`name` varchar(40) NOT NULL DEFAULT '',
`acc` varchar(40) DEFAULT NULL,
`ssid` mediumint(8) unsigned DEFAULT NULL,
`note_id` mediumint(8) unsigned DEFAULT NULL,
PRIMARY KEY  (`study_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

--
-- Table structure for table `allocation`
--

DROP TABLE IF EXISTS `allocation`;
CREATE TABLE `allocation` (
`study_id` smallint(5) unsigned DEFAULT NULL,
`individual_id` smallint(5) unsigned DEFAULT NULL,
`seq_centre_id` smallint(5) unsigned DEFAULT NULL,
PRIMARY KEY  (`study_id`,`individual_id`,`seq_centre_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;


--
-- Table structure for table `sample`
--

DROP TABLE IF EXISTS `sample`;
CREATE TABLE `sample` (
  `row_id` int unsigned NOT NULL auto_increment key,
  `sample_id` smallint(5) unsigned NOT NULL DEFAULT 0,
  `project_id` smallint(5) unsigned NOT NULL DEFAULT 0,
  `ssid` mediumint(8) unsigned DEFAULT NULL,
  `name` varchar(40) NOT NULL DEFAULT '',
  `hierarchy_name` varchar(40) NOT NULL DEFAULT '',
  `individual_id` smallint(5) unsigned DEFAULT NULL,
  `note_id` mediumint(8) unsigned DEFAULT NULL,
  `changed` datetime NOT NULL DEFAULT '0000-00-00',
  `latest` tinyint(1) DEFAULT '0',
  KEY  (`sample_id`),
  KEY `ssid` (`ssid`),
  KEY `latest` (`latest`),
  KEY `project_id` (`project_id`),
  KEY `name` (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

--
-- Table structure for table `seq_centre`
--

DROP TABLE IF EXISTS `seq_centre`;
CREATE TABLE `seq_centre` (
  `seq_centre_id` smallint(5) unsigned NOT NULL auto_increment,
  `name` varchar(40) NOT NULL DEFAULT '',
  PRIMARY KEY  (`seq_centre_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

--
-- Table structure for table `seq_tech`
--

DROP TABLE IF EXISTS `seq_tech`;
CREATE TABLE `seq_tech` (
  `seq_tech_id` smallint(5) unsigned NOT NULL auto_increment,
  `name` varchar(40) NOT NULL DEFAULT '',
  PRIMARY KEY  (`seq_tech_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

--
-- Table structure for table `submission`
--

DROP TABLE IF EXISTS `submission`;
CREATE TABLE `submission` (
  `submission_id` smallint(5) unsigned NOT NULL auto_increment,
  `date` datetime NOT NULL DEFAULT '0000-00-00',
  `name` varchar(40) NOT NULL DEFAULT '',
  `acc` varchar(40) DEFAULT NULL,
  PRIMARY KEY  (`submission_id`),
  UNIQUE KEY `acc` (`acc`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;


--
-- Views
--

DROP VIEW if EXISTS `latest_project`;
create view latest_project as select * from project where latest=true;
DROP VIEW if EXISTS `latest_sample`;
create view latest_sample as select * from sample where latest=true;
DROP VIEW if EXISTS `latest_library`;
create view latest_library as select * from library where latest=true;
DROP VIEW if EXISTS `latest_library_request`;
create view latest_library_request as select * from library_request where latest=true;
DROP VIEW if EXISTS `latest_seq_request`;
create view latest_seq_request as select * from seq_request where latest=true;
DROP VIEW if EXISTS `latest_lane`;
create view latest_lane as select * from lane where latest=true;
DROP VIEW if EXISTS `latest_file`;
create view latest_file as select * from file where latest=true;
DROP VIEW if EXISTS `latest_mapstats`;
create view latest_mapstats as select * from mapstats where latest=true;
