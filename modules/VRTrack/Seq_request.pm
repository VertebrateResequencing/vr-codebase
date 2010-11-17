package VRTrack::Seq_request; 

=head1 NAME

VRTrack::Seq_Request - Tracking Seq_Request object

=head1 SYNOPSIS
    my $seqrequest= VRTrack::Seq_request->new($vrtrack, $seqrequest_id);

    my $id = $seqrequest->id();
    my $status = $seqrequest->status();

=head1 DESCRIPTION

An object describing the tracked properties of a request.

=head1 AUTHOR

rn2@sanger.ac.uk (author)

=head1 METHODS

=cut

use strict;
use warnings;
use Carp qw(cluck confess);

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
               seq_request_id  => sub { $self->id(@_)},
	       library_id       => sub {$self->library_id(@_)},
	       multiplex_pool_id       => sub {$self->multiplex_pool_id(@_)},
	       ssid            => sub { $self->ssid(@_)},
	       seq_type     => sub { $self->seq_type(@_)},
	       seq_status     => sub { $self->seq_status(@_)});
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
  Example    : my $seq_request = VRTrack::Seq_request->create($vrtrack, $id)
  Description: Class method.  Creates new Seq_request object in the database.
  Returntype : VRTrack::Seq_request object
   
=cut


###############################################################################
# Class methods
###############################################################################


=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $seq_request->id();
	       $seq_request->id(104);
  Description: Get/Set for internal db ID of a request
  Returntype : integer

=cut


=head2 library_id

  Arg [1]    : library_id (optional)
  Example    : my $library_id = $seq_request->library_id();
               $seqrequest->library_id(104);
  Description: Get/Set for ID of a library
  Returntype : Internal ID integer

=cut

sub library_id {
    my $self = shift;
    return $self->_get_set('library_id', 'number', @_);
}


=head2 multiplex_pool_id

  Arg [1]    : multiplex_pool_id (optional)
  Example    : my $multiplex_pool_id = $Multiplex_Seq_request->multiplex_pool_id();
               $seqrequest->multiplex_pool_id(104);
  Description: Get/Set for ID of a library
  Returntype : Internal ID integer

=cut

sub multiplex_pool_id {
    my $self = shift;
    return $self->_get_set('multiplex_pool_id', 'number', @_);
}


=head2 ssid

  Arg [1]    : ssid (optional)
  Example    : my $ssid = $seq_request->ssid();
               $seq_request->ssid('104');
  Description: Get/Set for SequenceScape ID of a request
  Returntype : integer

=cut

=head2 seq_status

  Arg [1]    :  seq_status (optional)
  Example    : my $seq_status = $seqrequest->seq_status();
	       $seqrequest->seq_status('104');
  Description: Get/Set for seq request  sequencing status from SequenceScape
  Returntype : string

=cut

sub seq_status {
    my $self = shift;
    $self->_check_status_value('seq_status', @_);
    return $self->_get_set('seq_status', 'string', @_);
}

=head2 seq_type

  Arg [1]    :  seq_type (optional)
  Example    : my $seq_type = $seqrequest->seq_type();
	       $seqrequest->seq_type('Single ended sequencing');
  Description: Get/Set for sequencing  type from SequenceScape
  Returntype : string

=cut

sub seq_type {
    my $self = shift;
    $self->_check_status_value('seq_type', @_);
    return $self->_get_set('seq_type', 'string', @_);
}


=head2 changed

  Arg [1]    : changed (optional)
  Example    : my $changed = $seq_request->changed();
	       $seq_request->changed('2010-01-20');
  Description: Get/Set for request changed
  Returntype : string

=cut

=head2 lanes

  Arg [1]    : None
  Example    : my $lanes = $seq_request->lanes();
  Description: Returns a ref to an array of the lane objects that are associated with this seq_request.
  Returntype : ref to array of VRTrack::Seq_request objects

=cut

sub lanes{
    my $self = shift;
    return $self->_get_child_objects('VRTrack::Lane');
}


=head2 lane_ids

  Arg [1]    : None
  Example    : my $lane_ids = $seq_request->lane_ids();
  Description: Returns a ref to an array of the lane IDs that are associated with this sequencing request
  Returntype : ref to array of integer lane IDs

=cut

sub lane_ids {
    my $self = shift;
    return $self->_get_child_ids('VRTrack::Lane');
}

=head2 add_lane

  Arg [1]    : lane name
  Example    : my $newlane = $seq_request->add_lane('2631_3');
  Description: create a new lane, and if successful, return the object
  Returntype : VRTrack::Lane object

=cut

sub add_lane {
    my $self = shift;
    my $lane = $self->_add_child_object('new_by_name', 'VRTrack::Lane', @_);
    if ($lane && $self->library_id){ # non-multiplexed seq_req
        $lane->library_id($self->library_id);
        $lane->update;
    }
    return $lane;
}



=head2 get_lane_by_id

  Arg [1]    : lane internal id
  Example    : my $lane = $seq_request->get_lane_by_id(1930);
  Description: retrieve lane object by internal id
  Returntype : VRTrack::Lane object

=cut

sub get_lane_by_id {
    my $self = shift;
    return $self->_get_child_by_field_value('lanes', 'id', @_);
}


=head2 get_lane_by_ssid

  Arg [1]    : lane sequencescape id
  Example    : my $lane = $seq_request->get_lane_by_ssid(1930);
  Description: retrieve lane object by sequencescape id
  Returntype : VRTrack::Lane object

=cut

sub get_lane_by_ssid {
    my $self = shift;
    return $self->_get_child_by_field_value('lane', 'ssid', @_);
}


=head2 get_lane_by_name

  Arg [1]    : lane name
  Example    : my $lane = $seq_request->get_lane_by_name('My lane');
  Description: retrieve lane object by name
  Returntype : VRTrack::Lane object

=cut

sub get_lane_by_name {
    my $self = shift;
    return $self->_get_child_by_field_value('lanes', 'name', @_);
}


=head2 descendants

  Arg [1]    : none
  Example    : my $desc_objs = $obj->descendants();
  Description: Returns a ref to an array of all objects that are descendants of this object
  Returntype : arrayref of objects

=cut

sub _get_child_methods {
    return qw(lanes);
}



1;
