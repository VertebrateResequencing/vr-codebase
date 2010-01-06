package VRTrack::Request; 

=head1 NAME

VRTrack::Request - Sequence Tracking Request object

=head1 SYNOPSIS
    my $request= VRTrack::Request->new($vrtrack, $request_id);

    my $id = $request->id();
    my $status = $request->status();

=head1 DESCRIPTION

An object describing the tracked properties of a request.

=head1 CONTACT

jws@sanger.ac.uk (author)

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
               request_id      => sub { $self->id(@_)},
               library_id      => sub { $self->library_id(@_)},
               ssid            => sub { $self->ssid(@_)},
               seq_status      => sub { $self->seq_status(@_)},
               name            => sub { $self->name(@_)});
    
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

  Arg [1]    : none
  Example    : n/a
  Description: create does not apply to Request objects, since they can't be
               created?!
  Returntype : n/a

=cut

sub create {
    return undef;
}


###############################################################################
# Class methods
###############################################################################


=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $request->id();
	       $request->id(104);
  Description: Get/Set for internal db ID of a request
  Returntype : integer

=cut


=head2 library_id

  Arg [1]    : library_id (optional)
  Example    : my $library_id = $request->library_id();
	       $request->library_id('104');
  Description: Get/Set for internal library ID of a request
  Returntype : integer

=cut

sub library_id {
    my $self = shift;
    return $self->_get_set('library_id', 'number', @_);
}


=head2 ssid

  Arg [1]    : ssid (optional)
  Example    : my $ssid = $lib->ssid();
               $lib->ssid('104');
  Description: Get/Set for SequenceScape ID of a request
  Returntype : integer

=cut


=head2 seq_status

  Arg [1]    : seq_status (optional)
  Example    : my $seq_status = $request->seq_status();
	       $request->seq_status('104');
  Description: Get/Set for request sequencing status from SequenceScape
  Returntype : string

=cut

sub seq_status {
    my $self = shift;
    $self->_check_status_value('seq_status', @_);
    return $self->_get_set('seq_status', 'string', @_);
}


=head2 changed

  Arg [1]    : changed (optional)
  Example    : my $changed = $request->changed();
	       $request->changed('104');
  Description: Get/Set for request changed
  Returntype : string

=cut

1;
