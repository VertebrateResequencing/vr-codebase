package VRTrack::Request; 
# author: jws
=head1 NAME

VRTrack::Request - Sequence Tracking Request object

=head1 SYNOPSIS
    my $request= VRTrack::Request->new($vrtrack, $request_id);

    my $id = $request->id();
    my $status = $request->status();

=head1 DESCRIPTION

An object describing the tracked properties of a request.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;
no warnings 'uninitialized';
use VRTrack::Core_obj;
our @ISA = qw(VRTrack::Core_obj);

=head2 fields_dispatch

  Arg [1]    : none
  Example    : my $fieldsref = $req->fields_dispatch();
  Description: Returns hashref dispatch table keyed on database field
               Used internally for new and update methods
  Returntype : hashref

=cut

sub fields_dispatch {
    my $self = shift;
    my %fields = ( 
                'request_id'        => sub { $self->id(@_)},
                'library_id'        => sub { $self->library_id(@_)},
                'ssid'              => sub { $self->ssid(@_)},
                'seq_status'        => sub { $self->seq_status(@_)},
                'note_id'           => sub { $self->note_id(@_)},
                'changed'           => sub { $self->changed(@_)},
                'latest'            => sub { $self->is_latest(@_)},
                );

    return \%fields;
}


=head2 new_by_ssid

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : request sequencescape id
  Example    : my $request = VRTrack::Request->new_by_ssid($vrtrack, $ssid);
  Description: Class method. Returns latest Request object by ssid.  If no such ssid is in the database, returns undef
  Returntype : VRTrack::Request object

=cut

sub new_by_ssid {
    my ($class,$vrtrack, $ssid) = @_;
    die "Need to call with a vrtrack handle, ssid" unless ($vrtrack && $ssid);
    return $class->new_by_field_value($vrtrack, 'ssid',$ssid);
}


###############################################################################
# Class methods
###############################################################################

=head2 dirty

  Arg [1]    : boolean for dirty status
  Example    : $request->dirty(1);
  Description: Get/Set for request properties having been altered.
  Returntype : boolean

=cut

sub dirty {
    my ($self,$dirty) = @_;
    if (defined $dirty){
	$self->{_dirty} = $dirty ? 1 : 0;
    }
    return $self->{_dirty};
}



=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $request->id();
	       $request->id(104);
  Description: Get/Set for internal db ID of a request
  Returntype : integer

=cut

sub id {
    my ($self,$id) = @_;
    if (defined $id and $id != $self->{'id'}){
	$self->{'id'} = $id;
	$self->dirty(1);
    }
    return $self->{'id'};
}


=head2 library_id

  Arg [1]    : library_id (optional)
  Example    : my $library_id = $request->library_id();
	       $request->library_id('104');
  Description: Get/Set for internal library ID of a request
  Returntype : integer

=cut

sub library_id {
    my ($self,$library_id) = @_;
    if (defined $library_id and $library_id != $self->{'library_id'}){
	$self->{'library_id'} = $library_id;
	$self->dirty(1);
    }
    return $self->{'library_id'};
}


=head2 ssid

  Arg [1]    : ssid (optional)
  Example    : my $ssid = $lib->ssid();
               $lib->ssid('104');
  Description: Get/Set for SequenceScape ID of a request
  Returntype : integer

=cut

sub ssid {
    my ($self,$ssid) = @_;
    if (defined $ssid and $ssid != $self->{'ssid'}){
        $self->{'ssid'} = $ssid;
        $self->dirty(1);
    }
    return $self->{'ssid'};
}


=head2 seq_status

  Arg [1]    : seq_status (optional)
  Example    : my $seq_status = $request->seq_status();
	       $request->seq_status('104');
  Description: Get/Set for request sequencing status from SequenceScape
  Returntype : string

=cut

sub seq_status {
    my ($self,$seq_status) = @_;
    if (defined $seq_status and $seq_status ne $self->{'seq_status'}){
        my %allowed = map {$_ => 1} @{$self->list_enum_vals('request','seq_status')};
        unless ($allowed{lc($seq_status)}){
            die "'$seq_status' is not a defined seq_status";
        }
	$self->{'seq_status'} = $seq_status;
	$self->dirty(1);
    }
    return $self->{'seq_status'};
}


=head2 changed

  Arg [1]    : changed (optional)
  Example    : my $changed = $request->changed();
	       $request->changed('104');
  Description: Get/Set for request changed
  Returntype : string

=cut

sub changed {
    my ($self,$changed) = @_;
    if (defined $changed and $changed ne $self->{'changed'}){
	$self->{'changed'} = $changed;
	$self->dirty(1);
    }
    return $self->{'changed'};
}


1;
