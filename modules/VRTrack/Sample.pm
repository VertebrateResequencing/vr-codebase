package VRTrack::Sample;

=head1 NAME

VRTrack::Sample - Sequence Tracking Sample object

=head1 SYNOPSIS
    my $samp = VRTrack::Sample->new($vrtrack, $sample_id);

    #get arrayref of library objects in a sample
    my $libs = $sample->libraries();

    my $id = $sample->id();
    my $name = $sample->name();

=head1 DESCRIPTION

An object describing the tracked properties of a sample.

=head1 AUTHOR

jws@sanger.ac.uk (author)

=head1 METHODS

=cut

use strict;
use warnings;
use Carp qw(cluck confess);
use VRTrack::Library;
use VRTrack::Library_request;
use VRTrack::Individual;
use VRTrack::Allocations;

use base qw(VRTrack::Core_obj
            VRTrack::Hierarchy_obj
	    VRTrack::SequenceScape_obj);


=head2 fields_dispatch

  Arg [1]    : none
  Example    : my $fieldsref = $sample->fields_dispatch();
  Description: Returns hashref dispatch table keyed on database field
               Used internally for new and update methods
  Returntype : hashref

=cut

sub fields_dispatch {
    my $self = shift;
    
    my %fields = %{$self->SUPER::fields_dispatch()};
    %fields = (%fields,
               sample_id      => sub { $self->id(@_)},
               project_id     => sub { $self->project_id(@_)},
               ssid           => sub { $self->ssid(@_)},
               individual_id  => sub { $self->individual_id(@_)},
	       name           => sub { $self->name(@_)},
	       hierarchy_name => sub { $self->hierarchy_name(@_)});

    return \%fields;
}


###############################################################################
# Class methods
###############################################################################


=head2 new_by_name_project

  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : sample name
  Arg [3]    : project id
  Example    : my $sample = VRTrack::Sample->new_by_name_project($vrtrack, $name, $project_id)
  Description: Class method. Returns latest Sample object by name and
               project_id. If no such name is in the database, returns undef
  Returntype : VRTrack::Sample object

=cut

sub new_by_name_project {
    my ($class, $vrtrack, $name, $project_id) = @_;
    confess "Need to call with a vrtrack handle, name, project_id" unless ($vrtrack && $name && $project_id);
    if ( $vrtrack->isa('DBI::db') ) { confess "The interface has changed, expected vrtrack reference.\n"; }
    my $dbh = $vrtrack->{_dbh};
    my $history_sql = $class->_history_sql;
    my $sql = qq[select sample_id from sample where name=? and project_id = ? $history_sql];
    my $sth = $dbh->prepare($sql);
    
    my $id;
    if ($sth->execute($name, $project_id)){
        my $data = $sth->fetchrow_hashref;
        unless ($data){
            return undef;
        }
        $id = $data->{'sample_id'};
    }
    else{
        confess(sprintf('Cannot retrieve sample by $name, $project: %s', $DBI::errstr));
    }
    
    return $class->new($vrtrack, $id);
}


=head2 new_by_ssid

  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : sample sequencescape id
  Example    : my $sample = VRTrack::Sample->new_by_ssid($vrtrack, $ssid);
  Description: Class method. Returns latest Sample object by ssid.  If no such ssid is in the database, returns undef
  Returntype : VRTrack::Sample object

=cut


=head2 create
  
  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : name
  Arg [3]    : project id (optional)
  Example    : my $file = VRTrack::Sample->create($vrtrack, $name, $project_id)
  Description: Class method.  Creates new Sample object in the database.
               Overrides Core_obj method to allow allow creating samples with the 
               same name, but different project ids  .
  Returntype : VRTrack::Sample object
   
=cut

sub create {
    my ($class, $vrtrack, $name, $pid) = @_;
    confess "Need to call with a vrtrack handle" unless $vrtrack;
    confess "The interface has changed, expected vrtrack reference." if $vrtrack->isa('DBI::db');
    
    my $dbh = $vrtrack->{_dbh};
    my $table = $class->_class_to_table;
    
    # prevent adding an object with an existing name, if name supplied. In case of mapstats, the name is void
    if ($name && $class->is_name_in_database($vrtrack, $name, $name, $pid)){
        confess "Already a $table entry with value $name";
    }
    
    my $next_id;
    my $success = $vrtrack->transaction(sub {
	# insert a fake record to obtain a unique id (row_id)
	my $query = qq[INSERT INTO $table SET ${table}_id=0];
	my $sth   = $dbh->prepare($query) or confess qq[The query "$query" failed: $!];
	my $rv    = $sth->execute or confess qq[The query "$query" failed: $!];
	
	# now update the inserted record
	$next_id = $dbh->last_insert_id(undef, undef, $table, 'row_id') or confess "No last_insert_id? $!";
	
	if ($name) {
	    my $hierarchy_name;
	    
	    my $fieldsref = $class->fields_dispatch();
	    if ( exists($fieldsref->{hierarchy_name}) )
	    {
		$hierarchy_name = $name;
		$hierarchy_name =~ s/\W+/_/g;
	    }
	    
	    $name = qq[name='$name' ];
	    if ($hierarchy_name) {
		$name .= qq[, hierarchy_name='$hierarchy_name' ];
	    }
	}
	
	$query = qq[UPDATE $table SET ${table}_id=$next_id];
	if ($name){
	    $query .= qq[, $name ];     # add name, hierarchy_name clause
	}
	
	$query .= qq[, changed=now(), latest=true WHERE row_id=$next_id];
	$sth   = $dbh->prepare($query) or confess qq[The query "$query" failed: $!];
	$sth->execute or confess qq[The query "$query" failed: $!];
    });
    
    unless ($success) {
	confess $vrtrack->{transaction_error};
    }
    
    return $class->new($vrtrack, $next_id);
}

=head2 is_name_in_database

  Arg [1]    : sample name
  Arg [2]    : hierarchy name
  Arg [3]    : project id (optional)
  Example    : if(VRTrack::Sample->is_name_in_database($vrtrack, $name, $hname, $project_id)
  Description: Class method. Checks to see if a name or hierarchy name is already used in the sample table.
               Overrides Core_obj method.
  Returntype : boolean

=cut

sub is_name_in_database {
    my ($class, $vrtrack, $name, $hname, $pid) = @_;
    confess "Need to call with a vrtrack handle, name, hierarchy name" unless ($vrtrack && $name && $hname);
    if ($vrtrack->isa('DBI::db')) {
        confess "The interface has changed, expected vrtrack reference.\n";
    }
    
    my $table = $class->_class_to_table;
    
    my $dbh = $vrtrack->{_dbh};
    my $sql = qq[select ${table}_id from $table where latest=true and (name = ? or hierarchy_name = ?)];
    if ($pid) {
        $sql .= qq[ and project_id = $pid];
    }
    my $sth = $dbh->prepare($sql);
    
    my $already_used = 0;
    if ($sth->execute($name, $hname)) {
        my $data = $sth->fetchrow_hashref;
        if ($data) {
            $already_used = 1;
        }
    }
    else {
        confess "Cannot retrieve $table by $name: ".$DBI::errstr;
    }
    
    return $already_used;
}


###############################################################################
# Object methods
###############################################################################

=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $samp->id();
               $samp->id('104');
  Description: Get/Set for ID of a sample
  Returntype : Internal ID integer

=cut


=head2 project_id

  Arg [1]    : project_id (optional)
  Example    : my $project_id = $samp->project_id();
               $samp->project_id('104');
  Description: Get/Set for ID of a sample
  Returntype : SequenceScape ID (usu. integer)

=cut

sub project_id {
    my $self = shift;
    return $self->_get_set('project_id', 'number', @_);
}


=head2 hierarchy_name

  Arg [1]    : directory name (optional)
  Example    : my $hname = $sample->hierarchy_name();
  Description: Get/set sample hierarchy name.  This is the directory name
               (without path) that the sample will be named in a file hierarchy.
  Returntype : string

=cut

sub hierarchy_name {
    my $self = shift;
    return $self->_get_set('hierarchy_name', 'string', @_);
}


=head2 name

  Arg [1]    : name (optional)
  Example    : my $name = $samp->name();
               $samp->name('104');
  Description: Get/Set for sample name
  Returntype : string

=cut

sub name {
    # we can't be a Named_obj since we don't allow new_by_name(), so have to
    # implement this ourselves
    my $self = shift;
    return $self->_get_set('name', 'string', @_);
}


=head2 ssid

  Arg [1]    : ssid (optional)
  Example    : my $ssid = $samp->ssid();
               $samp->ssid(104);
  Description: Get/Set for sample SequenceScape ID
  Returntype : SequenceScape ID integer

=cut


=head2 individual_id

  Arg [1]    : individual_id (optional)
  Example    : my $individual_id = $samp->individual_id();
               $samp->individual_id(123);
  Description: Get/Set for sample internal individual_id
  Returntype : integer

=cut

sub individual_id {
    my $self = shift;
    return $self->_get_set('individual_id', 'number', @_);
}


=head2 individual

  Arg [1]    : individual name (optional)
  Example    : my $individual = $samp->individual();
               $samp->individual('NA19820');
  Description: Get/Set for sample individual.  Lazy-loads individual object from $self->individual_id.  If a individual name is supplied, then individual_id is set to the corresponding individual in the database.  If no such individual exists, returns undef.  Use add_individual to add a individual in this case.
  Returntype : VRTrack::Individual object

=cut

sub individual {
    my $self = shift;
    return $self->_get_set_child_object('get_individual_by_name', 'VRTrack::Individual', @_);
}


=head2 add_individual

  Arg [1]    : individual name
  Example    : my $ind = $samp->add_individual('NA19820');
  Description: create a new individual, and if successful, return the object
  Returntype : VRTrack::Library object

=cut

sub add_individual {
    my $self = shift;
    return $self->_create_child_object('get_individual_by_name', 'VRTrack::Individual', @_);
}


=head2 get_individual_by_name

  Arg [1]    : individual_name
  Example    : my $ind = $samp->get_individual_by_name('NA19820');
  Description: Retrieve a VRTrack::Individual object by name
  Returntype : VRTrack::Individual object

=cut

sub get_individual_by_name {
    my ($self,$name) = @_;
    return VRTrack::Individual->new_by_name($self->{vrtrack}, $name);
}


=head2 libraries

  Arg [1]    : None
  Example    : my $libraries = $sample->libraries();
  Description: Returns a ref to an array of the sample objects that are associated with this sample.
  Returntype : ref to array of VRTrack::Sample objects

=cut

sub libraries {
    my $self = shift;
    return $self->_get_child_objects('VRTrack::Library');
}


=head2 library_ids

  Arg [1]    : None
  Example    : my $library_ids = $sample->library_ids();
  Description: Returns a ref to an array of the library IDs that are associated with this sample
  Returntype : ref to array of integer library IDs

=cut

sub library_ids {
    my $self = shift;
    return $self->_get_child_ids('VRTrack::Library');
}


=head2 add_library

  Arg [1]    : library name
  Example    : my $newlib = $samp->add_library('NOD_500_SLX_1');
  Description: create a new library, and if successful, return the object
  Returntype : VRTrack::Library object

=cut

sub add_library {
    my $self = shift;
    return $self->_add_child_object('new_by_name', 'VRTrack::Library', @_);
}


=head2 get_library_by_id

  Arg [1]    : library internal id
  Example    : my $library = $sam->get_library_by_id(1930);
  Description: retrieve library object by internal id
  Returntype : VRTrack::Library object

=cut

sub get_library_by_id {
    my $self = shift;
    return $self->_get_child_by_field_value('libraries', 'id', @_);
}


=head2 get_library_by_ssid

  Arg [1]    : library sequencescape id
  Example    : my $library = $sam->get_library_by_ssid(1930);
  Description: retrieve library object by sequencescape id
  Returntype : VRTrack::Library object

=cut

sub get_library_by_ssid {
    my $self = shift;
    return $self->_get_child_by_field_value('libraries', 'ssid', @_);
}


=head2 get_library_by_name

  Arg [1]    : library name
  Example    : my $library = $sam->get_library_by_name('My library');
  Description: retrieve library object by name
  Returntype : VRTrack::Library object

=cut

sub get_library_by_name {
    my $self = shift;
    return $self->_get_child_by_field_value('libraries', 'name', @_);
}


=head2 library requests

  Arg [1]    : None
  Example    : my $library_requests = $sample->library_requests();
  Description: Returns a ref to an array of the Library_Reqest objects that are associated with this sample.
  Returntype : ref to array of VRTrack::Library_Request objects

=cut

sub library_requests {
    my $self = shift;
    return $self->_get_child_objects('VRTrack::Library_request');
}


=head2 library_request_ids

  Arg [1]    : None
  Example    : my $library_request_ids = $sample->library_request_ids();
  Description: Returns a ref to an array of the library request IDs that are associated with this sample
  Returntype : ref to array of integer library request IDs

=cut

sub library_request_ids {
    my $self = shift;
    return $self->_get_child_ids('VRTrack::Library_request');
}


=head2 add_library_request

  Arg [1]    : library request id
  Example    : my $newlibrequest = $samp->add_library_request('1234');
  Description: create a new library request , and if successful, return the object
  Returntype : VRTrack::Library_Request object

=cut

sub add_library_request {
    my ($self,$ssid) = @_;
    return $self->_add_child_object('new_by_ssid','VRTrack::Library_request',$ssid);
}


=head2 get_library_request_by_id

  Arg [1]    : library_request internal id
  Example    : my $library_request = $sam->get_library_request_by_id(1930);
  Description: retrieve library_request  object by internal id
  Returntype : VRTrack::Library_Request object

=cut

sub get_library_request_by_id {
    my $self = shift;
    return $self->_get_child_by_field_value('library_requests', 'id', @_);
}


=head2 get_library_request_by_ssid

  Arg [1]    : library_request sequencescape id
  Example    : my $library_request = $sam->get_library_request_by_ssid(1930);
  Description: retrieve library request object by sequencescape id
  Returntype : VRTrack::Library_Request object

=cut

sub get_library_request_by_ssid {
    my $self = shift;
    return $self->_get_child_by_field_value('library_requests', 'ssid', @_);
}

=head2 get_allocated_seq_centres

  Arg [1]    : None
  Example    : my @centres = $sam->get_allocated_seq_centres();
  Description: retrieve sequencing centres allocated to this sample (in this project)
  Returntype : arrayref of VRTrack::Seq_centre objects

=cut

sub get_allocated_seq_centres {
    my ($self) = @_;
    unless ($self->{'seq_centres'}){
        my $allocs = VRTrack::Allocations->new($self->{vrtrack});
	eval "require VRTrack::Project;"; # (we avoid using this since Project uses us)
        my $project = VRTrack::Project->new($self->{vrtrack},$self->project_id);
        unless ($allocs && $project){
            confess "Can't retrieve Allocations and Project";
        }
        my $centres = $allocs->get_centres_for_study_ind($project->study->id,$self->individual->id);
        $self->{'seq_centres'} = $centres;
    }
    return $self->{'seq_centres'};
}


=head2 is_sanger_sample

  Arg [1]    : None
  Example    : my $is_ours = $sam->is_sanger_sample();
  Description: checks that this sample in this project has been allocated to the Sanger (SC).
  Returntype : Boolean

=cut

sub is_sanger_sample {
    my ($self) = @_;
    my $is_ours = 0;
    if (grep {$_->name eq 'SC'} @{$self->get_allocated_seq_centres}){
        $is_ours = 1;
    }
    return $is_ours;
}


=head2 changed

  Arg [1]    : changed (optional)
  Example    : my $changed = $sample->changed();
               $sample->changed('20080810123000');
  Description: Get/Set for sample changed
  Returntype : string

=cut


=head2 descendants

  Arg [1]    : none
  Example    : my $desc_objs = $obj->descendants();
  Description: Returns a ref to an array of all objects that are descendants of this object
  Returntype : arrayref of objects

=cut

sub _get_child_methods {
    return qw(library_requests libraries);
}

1;
