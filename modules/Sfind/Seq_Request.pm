package Sfind::Seq_Request; 
=head1 NAME

Sfind::Request - Sequence Tracking Request object

=head1 SYNOPSIS
    my $seqrequest= Sfind::Seq_Request->new($dbh, $request_id);

    my $id = $seqrequest->id();
    my $status = $seqrequest->status();

=head1 DESCRIPTION

An object describing the tracked properties of a sequencing request.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;
no warnings 'uninitialized';
use Sfind::Lane;

=head2 new

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : request id
  Example    : my $seqrequest= Sfind::Sfind->new($dbh, $id)
  Description: Returns Seq_Request object by request_id
  Returntype : Sfind::Seq_Request object

=cut

sub new {
    my ($class,$dbh, $id) = @_;
    die "Need to call with a db handle and id" unless ($dbh && $id);
    my $self = {};
    bless ($self, $class);
    $self->{_dbh} = $dbh;

    my $sql = qq[select requests_new.asset_id,assets.name,requests_new.type,requests_new.state,requests_new.created_at,requests_new.read_length from 
		    requests_new join assets on (requests_new.asset_id=assets.asset_id)
		    where requests_new.request_id = ?];
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
    $self->created($data->[4]);
    $self->read_len($data->[5]);

    # cancelled requests are currently in the database with no status
    if ($data->[3]){
        $self->status($data->[3]);
    }
    else {
        $self->status('cancelled');
    }

    return $self;
}


=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $seqrequest->id();
	       $seqrequest->id('104');
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
  Example    : my $created = $seqrequest->created();
	       $seqrequest->created('2008-10-24 11:40:59');
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
  Example    : my $type = $seqrequest->type();
	       $seqrequest->type('Paired end sequencing');
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
  Example    : my $library_id = $seqrequest->library_id();
	       $seqrequest->library_id('104');
  Description: Get/Set for library ID of a request
		This is the multiplex_tube_asset_id if the library is a multiplex library
		else
		it is the library tube asset id
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
  Example    : my $library_name = $seqrequest->library_name();
	       $seqrequest->library_name('foo');
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
  Example    : my $read_len = $seqrequest->read_len();
	       $seqrequest->read_len(54);
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
  Example    : my $status = $seqrequest->status();
	       $seqrequest->status('pending');
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
  Example    : my $lanes = $seqrequest->lanes();
  Description: Returns a ref to an array of the file objects that are associated with this request.
  Returntype : ref to array of Sfind::File objects

=cut

sub lanes {
    my ($self) = @_;
    unless ($self->{'lanes'}){
	my @lanes;
    	foreach my $id (@{$self->lane_ids()}){
	    my $obj = Sfind::Lane->new($self->{_dbh},$id);
	    push @lanes, $obj;
	}
	$self->{'lanes'} = \@lanes;
    }

    return $self->{'lanes'};
}


=head2 lane_ids

  Arg [1]    : None
  Example    : my $lane_ids = $seqrequest->lane_ids();
  Description: Returns a ref to an array of the file names that are associated with this request
  Returntype : ref to array of file names

=cut

sub lane_ids {
    my ($self) = @_;
    unless ($self->{'lane_ids'}){
	my $sql = qq[select id_npg_information from npg_information n, library l where l.request_id=? and l.batch_id =n.batch_id and l.position = n.position and (n.id_run_pair=0 or n.id_run_pair is null);];
	my @lanes;
	my $sth = $self->{_dbh}->prepare($sql);

	$sth->execute($self->id);
	foreach(@{$sth->fetchall_arrayref()}){
	    push @lanes, $_->[0];
	}
	@lanes = sort {$a <=> $b} @lanes;

	$self->{'lane_ids'} = \@lanes;
    }
 
    return $self->{'lane_ids'};
}

=head2 sequenced_bases

  Arg [1]    : none
  Example    : my $tot_bp = $lib->sequenced_bases();
  Description: the total number of sequenced bases on this library.
		This is the sum of the bases from the fastq files associated
		with this library in NPG.
  Returntype : integer

=cut

sub sequenced_bases {
    my ($self, $id) = @_;
    unless ($self->{'seq_bases'}){
	$self->{'seq_bases'} = 0;
	foreach my $lane(@{$self->lanes}){
	    $self->{'seq_bases'} += $lane->basepairs;
	}
    }
    return $self->{'seq_bases'};
}


1;
