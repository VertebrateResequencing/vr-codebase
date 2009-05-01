package Sfind::Request; 
=head1 NAME

Sfind::Request - Sequence Tracking Request object

=head1 SYNOPSIS
    my $request= Sfind::Request->new($dbh, $request_id);

    my $id = $request->id();
    my $status = $request->status();

=head1 DESCRIPTION

An object describing the tracked properties of a request.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut

use Sfind::Lane;

=head2 new

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : request id
  Example    : my $request= Sfind::Sfind->new($dbh, $id)
  Description: Returns Request object by request_id
  Returntype : Sfind::Request object

=cut

sub new {
    my ($class,$dbh, $id) = @_;
    die "Need to call with a db handle and id" unless ($dbh && $id);
    my $self = {};
    bless ($self, $class);
    $self->{_dbh} = $dbh;

    my $sql = qq[select item_id, item_name, type, state, created_at, read_length from requests where request_id = ?];
    my $sth = $self->{_dbh}->prepare($sql);

    $sth->execute($id);
    # note, while-ing over this array would give all attempts
    my $data = $sth->fetchall_arrayref()->[0];
    unless ($data){
	return undef;
    }
    $self->id($id);
    $self->library_id($data->[0]);
    $self->library_name($data->[1]);
    $self->type($data->[2]);
    $self->status($data->[3]);
    $self->created($data->[4]);
    $self->read_len($data->[5]);
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


=head2 library_id

  Arg [1]    : library_id (optional)
  Example    : my $library_id = $request->library_id();
	       $request->library_id('104');
  Description: Get/Set for library ID of a request
  Returntype : SequenceScape ID integer

=cut

sub library_id {
    my ($self,$library_id) = @_;
    if (defined $library_id and $library_id ne $self->{'library_id'}){
	$self->{'library_id'} = $library_id;
    }
    return $self->{'library_id'};
}


=head2 library_name

  Arg [1]    : library_name (optional)
  Example    : my $library_name = $request->library_name();
	       $request->library_name('foo');
  Description: Get/Set for library name of request
  Returntype : SequenceScape name

=cut

sub library_name {
    my ($self,$library_name) = @_;
    if (defined $library_name and $library_name ne $self->{'library_name'}){
	$self->{'library_name'} = $library_name;
    }
    return $self->{'library_name'};
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



=head2 lanes

  Arg [1]    : None
  Example    : my $lanes = $request->lanes();
  Description: Returns a ref to an array of the file objects that are associated with this request.
  Returntype : ref to array of Sfind::File objects

=cut

sub lanes {
    my ($self) = @_;
    unless ($self->{'lanes'}){
	my @lanes;
    	foreach my $id (@{$self->lane_ids()}){
	    my $obj = Sfind::File->new($self->{_dbh},$id);
	    push @lanes, $obj;
	}
	$self->{'lanes'} = \@lanes;
    }

    return $self->{'lanes'};
}


=head2 lane_ids

  Arg [1]    : None
  Example    : my $lane_ids = $request->lane_ids();
  Description: Returns a ref to an array of the file names that are associated with this request
  Returntype : ref to array of file names

=cut

sub lane_ids {
    my ($self) = @_;
    unless ($self->{'lane_ids'}){
	my $sql = qq[select distinct(name) from file where request_id=?];
	my @files;
	my $sth = $self->{_dbh}->prepare($sql);

	if ($sth->execute($self->id)){
	    foreach(@{$sth->fetchall_arrayref()}){
		push @files, $_->[0];
	    }
	}
	else{
	    die(sprintf('Cannot retrieve files: %s', $DBI::errstr));
	}

	$self->{'lane_ids'} = \@files;
    }
 
    return $self->{'lane_ids'};
}


1;
