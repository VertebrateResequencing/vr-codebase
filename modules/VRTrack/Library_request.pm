package VRTrack::Library_request; 

=head1 NAME

VRTrack::Library_Request - Tracking Library_Request object

=head1 SYNOPSIS
    my $librequest= VRTrack::Library_Request->new($vrtrack, $librequest_id);

    my $id = $librequest->id();
    my $status = $librequest->status();

=head1 DESCRIPTION

An object describing the tracked properties of a library_request.

=head1 AUTHOR

rn2@sanger.ac.uk (author)

=head1 METHODS

=cut

use strict;
use warnings;
use Carp qw(cluck confess);
use VRTrack::Sample;
use VRTrack::Study;

use base qw(VRTrack::Core_obj
	    VRTrack::SequenceScape_obj);


=head2 fields_dispatch

  Arg [1]    : none
  Example    : my $fieldsref = $req->fields_dispatch();
  Description: Returns hashref dispatch table keyed on database field
               Used internally for new and update methods
  Returntype : hashref

=cut

sub fields_dispatch {
    my $self = shift;
    
    my %fields = %{$self->SUPER::fields_dispatch()};
    %fields = (%fields,
               library_request_id  => sub { $self->id(@_)},
	       sample_id       => sub {	$self->sample_id(@_)},
	       ssid            => sub { $self->ssid(@_)},
	       prep_status     => sub { $self->prep_status(@_)});
    return \%fields;
}


=head2 new_by_ssid

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : request sequencescape id
  Example    : my $request = VRTrack::Request->new_by_ssid($vrtrack, $ssid);
  Description: Class method. Returns latest Request object by ssid.  If no such ssid is in the database, returns undef
  Returntype : VRTrack::Request object

=cut


=head2 create
  
  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : id
  Example    : my $library_request = VRTrack::Library_request->create($vrtrack, $id)
  Description: Class method.  Creates new Library_request object in the database.
  Returntype : VRTrack::Library_request object
   
=cut


###############################################################################
# Class methods
###############################################################################


=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $librequest->id();
	       $librequest->id(104);
  Description: Get/Set for internal db ID of a library_request
  Returntype : integer

=cut


=head2 sample_id

  Arg [1]    : sample_id (optional)
  Example    : my $sample_id = $librequest->sample_id();
               $librequest->sample_id(104);
  Description: Get/Set for ID of a library
  Returntype : Internal ID integer

=cut

sub sample_id {
    my $self = shift;
    return $self->_get_set('sample_id', 'number', @_);
}


=head2 ssid

  Arg [1]    : ssid (optional)
  Example    : my $ssid = $librequest->ssid();
               $librequest->ssid('104');
  Description: Get/Set for SequenceScape ID of a request
  Returntype : integer

=cut


=head2 prep_status

  Arg [1]    :  prep_status (optional)
  Example    : my $prep_status = $librequest->prep_status();
	       $librequest->prep_status('104');
  Description: Get/Set for library request  preparation status from SequenceScape
  Returntype : string

=cut

sub prep_status {
    my $self = shift;
    $self->_check_status_value('prep_status', @_);
    return $self->_get_set('prep_status', 'string', @_);
}


=head2 changed

  Arg [1]    : changed (optional)
  Example    : my $changed = $librequest->changed();
	       $librequest->changed('104');
  Description: Get/Set for request changed
  Returntype : string

=cut

=head2 libraries

  Arg [1]    : None
  Example    : my $libraries = $library_request->libraries();
  Description: Returns a ref to an array of the library objects that are associated with this library_request.
  Returntype : ref to array of VRTrack::Library_request objects

=cut

sub libraries {
    my $self = shift;
    return $self->_get_child_objects('VRTrack::Library');
}


=head2 library_ids

  Arg [1]    : None
  Example    : my $library_ids = $library_request->library_ids();
  Description: Returns a ref to an array of the library IDs that are associated with this library_request
  Returntype : ref to array of integer library IDs

=cut

sub library_ids {
    my $self = shift;
    return $self->_get_child_ids('VRTrack::Library');
}


=head2 add_library

  Arg [1]    : library name
  Example    : my $newlib = $library_request->add_library('NOD_500_SLX_1');
  Description: create a new library, and if successful, return the object
  Returntype : VRTrack::Library object

=cut

sub add_library {
    my $self = shift;
    my $lib = $self->_add_child_object('new_by_name', 'VRTrack::Library', @_);
    if ($lib){
        $lib->sample_id($self->sample_id);
        $lib->update;
    }
    return $lib;
}

=head2 get_library_by_id

  Arg [1]    : library internal id
  Example    : my $library = $librequest->get_library_by_id(1930);
  Description: retrieve library object by internal id
  Returntype : VRTrack::Library object

=cut

sub get_library_by_id {
    my $self = shift;
    return $self->_get_child_by_field_value('libraries', 'id', @_);
}


=head2 get_library_by_ssid

  Arg [1]    : library sequencescape id
  Example    : my $library = $librequest->get_library_by_ssid(1930);
  Description: retrieve library object by sequencescape id
  Returntype : VRTrack::Library object

=cut

sub get_library_by_ssid {
    my $self = shift;
    return $self->_get_child_by_field_value('libraries', 'ssid', @_);
}


=head2 get_library_by_name

  Arg [1]    : library name
  Example    : my $library = $librequest->get_library_by_name('My library');
  Description: retrieve library object by name
  Returntype : VRTrack::Library object

=cut

sub get_library_by_name {
    my $self = shift;
    return $self->_get_child_by_field_value('libraries', 'name', @_);
}

=head2 descendants

  Arg [1]    : none
  Example    : my $desc_objs = $obj->descendants();
  Description: Returns a ref to an array of all objects that are descendants of this object
  Returntype : arrayref of objects

=cut

sub _get_child_methods {
    return qw(libraries);
}

1;
