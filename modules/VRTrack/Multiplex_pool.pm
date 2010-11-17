package VRTrack::Multiplex_pool; 

=head1 NAME

VRTrack::Multiplex_pool - Sequence Tracking Multiplex_library object

=head1 SYNOPSIS
    my $multiplex_lib = VRTrack::Multiplex_library->new_by_ssid($vrtrack,$ssid);

=head1 DESCRIPTION

An object describing the tracked properties of a library.

=head1 AUTHOR

rn2@sanger.ac.uk (author)

=head1 METHODS

=cut

use strict;
use warnings;
use Carp qw(cluck confess);
use VRTrack::Lane;
use VRTrack::Seq_request;
use VRTrack::Library_type;
use VRTrack::Seq_centre;
use VRTrack::Seq_tech;

use base qw(VRTrack::SequenceScape_obj
            VRTrack::Named_obj);


=head2 fields_dispatch

  Arg [1]    : none
  Example    : my $fieldsref = $lib->fields_dispatch();
  Description: Returns hashref dispatch table keyed on database field
               Used internally for new and update methods
  Returntype : hashref

=cut

sub fields_dispatch {
    my $self = shift;
    
    my %fields = %{$self->SUPER::fields_dispatch()};
    %fields = (%fields, 
               multiplex_pool_id    => sub { $self->id(@_)},
               ssid                 => sub { $self->ssid(@_)},
               name                 => sub { $self->name(@_)},
               );
    return \%fields;
}


###############################################################################
# Class methods
###############################################################################

=head2 new_by_ssid

  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : library sequencescape id
  Example    : my $multiplex_library = VRTrack::Multiplex_library->new_by_ssid($vrtrack, $ssid);
  Description: Class method. Returns latest Library object by ssid.  If no such ssid is in the database, returns undef
  Returntype : VVRTrack::Multiplex_library object

=cut

=head2 create
  
  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : ssid
  Example    : my $file = VRTrack::Multiplex_library->create($vrtrack, $ssid)
  Description: Class method.  Creates new Multiplex Library object in the database.
  Returntype : VRTrack::Multiplex_library object
   
=cut

sub create{
    my ($class, $vrtrack, $ssid) = @_;
    my $obj = $class->new_by_ssid($vrtrack,$ssid);
    if ($obj) {
        cluck "$class ($ssid) is already present in the database\n";
        return;
    }

    return $class->SUPER::create($vrtrack,'ssid'=>$ssid);
}

###############################################################################
# Object methods
###############################################################################

=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $multiplex_lib->id();
               $multiplex_lib->id(104);
  Description: Get/Set for internal db ID of a multiplex_library
  Returntype : integer

=cut


=head2 ssid

  Arg [1]    : ssid (optional)
  Example    : my $ssid = $multiplex_lib->ssid();
               $multiplex_lib->ssid('104');
  Description: Get/Set for SequenceScape ID of a multiplex_library
  Returntype : integer

=cut

=head2 library_multiplex_pools

  Arg [1]    : None
  Example    : my $library_multiplex_pools = $multiplex->library_multiplex_pools();
  Description: Returns a ref to an array of the library_multiplex_pool objects that are associated with this multiplex pool.
  Returntype : ref to array of VRTrack::Library_Multiplex_pool objects

=cut

sub library_multiplex_pools {
    my $self = shift;
    return $self->_get_child_objects('VRTrack::Library_Multiplex_pool');
}


=head2 library_multiplex_pools

  Arg [1]    : None
  Example    : my $library_multiplex_pool_ids = $multiplex->library_multiplex_pool_ids();
  Description: Returns a ref to an array of the library_multiplex_pool IDs that are associated with this multiplex library
  Returntype : ref to array of integer library_multiplex_pool IDs

=cut

sub library_multiplex_pool_ids {
    my $self = shift;
    return $self->_get_child_ids('VRTrack::Library_Multiplex_pool');
}


=head2 sequencing requests

  Arg [1]    : None
  Example    : my $seq_requests = $multiplex_lib->seq_requests();
  Description: Returns a ref to an array of the Seq_Request objects that are associated with this multiplex_library.
  Returntype : ref to array of VRTrack::Seq_Request objects

=cut

sub seq_requests {
    my $self = shift;
    return $self->_get_child_objects('VRTrack::Seq_request');
}


=head2 seq_request_ids

  Arg [1]    : None
  Example    : my $multiplex_Seq_Request_ids = $multiplex_lib->seq_Request_ids();
  Description: Returns a ref to an array of the Seq_Request IDs that are associated with this multiplex_library
  Returntype : ref to array of integer Seq_Request IDs

=cut

sub seq_request_ids {
    my $self = shift;
    return $self->_get_child_ids('VRTrack::Seq_request');
}


=head2 add_seq_request

  Arg [1]    : sequencing request id
  Example    : my $newseqrequest = $multiplex_lib->add_seq_request('12121212');
  Description: create a new Seq_Request , and if successful, return the object
  Returntype : VRTrack::Seq_Request object

=cut

sub add_seq_request {
    my ($self,$ssid) = @_;
    return $self->_add_child_object('new_by_ssid','VRTrack::Seq_request',$ssid);
}


=head2 get_seq_request_by_id

  Arg [1]    : Seq_Request internal id
  Example    : my $lseq_request = $multiplex_lib->get_seq_request_by_id(1930);
  Description: retrieve Seq_Request  object by internal id
  Returntype : VRTrack::Seq_Request object

=cut

sub get_seq_request_by_id {
    my $self = shift;
    return $self->_get_child_by_field_value('seq_requests', 'id', @_);
}


=head2 get_seq_request_by_ssid

  Arg [1]    : Seq_Request sequencescape id
  Example    : my $seq_request = $multiplex_lib->get_seq_request_by_ssid(1930);
  Description: retrieve Seq_Request object by sequencescape id
  Returntype : VRTrack::Seq_Request object

=cut

sub get_seq_request_by_ssid {
    my $self = shift;
    return $self->_get_child_by_field_value('seq_requests', 'ssid', @_);
}


=head2 changed

  Arg [1]    : timestamp (optional)
  Example    : my $changed = $multiplex_lib->changed();
               $multiplex_lib->changed('20080810123000');
  Description: Get/Set for multiplexlibrary changed
  Returntype : string

=cut

1;
