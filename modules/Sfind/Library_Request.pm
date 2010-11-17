package Sfind::Library_Request; 
=head1 NAME

Sfind::Library_Request - Sequence Tracking Library_Request object

=head1 SYNOPSIS
    my $librequest= Sfind::Library_Request->new($dbh, $request_id);

    my $id = $librequest->id();
    my $status = $librequest->status();

=head1 DESCRIPTION

An object describing the tracked properties of a library request.

=head1 CONTACT

rn2@sanger.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;
no warnings 'uninitialized';

=head2 new

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : library request id
  Example    : my $librequest= Sfind::Sfind->new($dbh, $id)
  Description: Returns Library_Request object by lib_request_id
  Returntype : Sfind::Library_Request object

=cut

sub new {
    my ($class,$dbh, $id) = @_;
    die "Need to call with a db handle and id" unless ($dbh && $id);
    my $self = {};
    bless ($self, $class);
    $self->{_dbh} = $dbh;

    my $sql = qq[select type, state, created_at, read_length from requests_new where request_id = ?];
    my $sth = $self->{_dbh}->prepare($sql);

    $sth->execute($id);
    # note, while-ing over this array would give all attempts
    my $data = $sth->fetchall_arrayref()->[0];
    unless ($data){
	return undef;
    }
    $self->id($id);
    $self->type($data->[0]);
    $self->created($data->[2]);
    $self->read_len($data->[3]);

    # cancelled requests are currently in the database with no status
    if ($data->[1]){
        $self->status($data->[1]);
    }
    else {
        $self->status('cancelled');
    }

    return $self;
}


=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $request->id();
	       $request->id('104');
  Description: Get/Set for ID of a request
  Returntype : SequenceScape ID (usu. integer)

=cut

sub id {
    my ($self,$id) = @_;
    if (defined $id and $id ne $self->{'id'}){
	$self->{'id'} = $id;
    }
    return $self->{'id'};
}


=head2 created

  Arg [1]    : created (optional)
  Example    : my $created = $request->created();
	       $request->created('2008-10-24 11:40:59');
  Description: Get/Set for created timestamp
  Returntype : timestamp string

=cut

sub created {
    my ($self,$created) = @_;
    if (defined $created and $created ne $self->{'created'}){
	$self->{'created'} = $created;
    }
    return $self->{'created'};
}


=head2 type

  Arg [1]    : type (optional)
  Example    : my $type = $request->type();
	       $request->type('Paired end sequencing');
  Description: Get/Set for type of request
  Returntype : sequencescape request type string

=cut

sub type {
    my ($self,$type) = @_;
    if (defined $type and $type ne $self->{'type'}){
	$self->{'type'} = $type;
    }
    return $self->{'type'};
}

=head2 read_len

  Arg [1]    : read_len (optional)
  Example    : my $read_len = $request->read_len();
	       $request->read_len(54);
  Description: Get/Set for request read_len
  Returntype : integer

=cut

sub read_len {
    my ($self,$read_len) = @_;
    if (defined $read_len and $read_len ne $self->{'read_len'}){
	$self->{'read_len'} = $read_len;
    }
    return $self->{'read_len'};
}


=head2 status

  Arg [1]    : status (optional)
  Example    : my $status = $request->status();
	       $request->status('pending');
  Description: Get/Set for request status
  Returntype : string

=cut

sub status {
    my ($self,$status) = @_;
    if (defined $status and $status ne $self->{'status'}){
	$self->{'status'} = $status;
    }
    return $self->{'status'};
}

=head2 libraries

  Arg [1]    : None
  Example    : my $libraries = $library_request->libraries();
  Description: Returns a ref to an array of the library objects that are associated with this library_request.
  Returntype : ref to array of Sfind::Library objects

=cut

sub libraries {
    my ($self) = @_;

    unless ($self->{'libraries'}){
	my @libraries;
    	foreach my $id (@{$self->library_ids()}){
	    my $obj = Sfind::Library->new($self->{_dbh},$id);
	    push @libraries, $obj; 
	}
	@libraries = sort {$a <=> $b} @libraries;
	$self->{'libraries'} = \@libraries;
    }

    return $self->{'libraries'};
}


=head2 library_ids

  Arg [1]    : None
  Example    : my $library_ids = $library_request->library_ids();
  Description: Returns a ref to an array of the library IDs that are associated with this library_request
  Returntype : ref to array of integer library IDs

=cut

sub library_ids {
    my ($self) = @_;


    # check if this is a multiplex library creation request
    # return the asset_ids as as library_ids
  
    # select all library_tube asset ids associated with this request
    my $sql= qq[select target_asset_id from requests_new where request_id=?];
          
    unless ($self->{'library_ids'}){
	my @libraries;
	my $sth = $self->{_dbh}->prepare($sql);

	$sth->execute($self->id);
	foreach(@{$sth->fetchall_arrayref()}){
	    push @libraries, $_->[0] if $_->[0];
	}
	$self->{'library_ids'} = \@libraries;
    }
    
    return $self->{'library_ids'};
}

1;
