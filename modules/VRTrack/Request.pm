package VRTrack::Request; 
=head1 NAME

VRTrack::Request - Sequence Tracking Request object

=head1 SYNOPSIS
    my $request= VRTrack::Request->new($dbh, $request_id);

    my $id = $request->id();
    my $status = $request->status();

=head1 DESCRIPTION

An object describing the tracked properties of a request.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut

use constant DBI_DUPLICATE => '1062';
use VRTrack::Utils;

=head2 new

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : request id
  Example    : my $request= VRTrack::Request->new($dbh, $id)
  Description: Returns Request object by request_id
  Returntype : VRTrack::Request object

=cut

sub new {
    my ($class,$dbh, $id) = @_;
    die "Need to call with a db handle and id" unless ($dbh && $id);
    my $self = {};
    bless ($self, $class);
    $self->{_dbh} = $dbh;

    my $sql = qq[select request_id, library_id, ssid, seq_status, changed, latest from request where request_id = ? and latest = true];
    my $sth = $self->{_dbh}->prepare($sql);

    if ($sth->execute($id)){
        my $data = $sth->fetchrow_hashref;
	unless ($data){
	    return undef;
	}
	$self->id($data->{'request_id'});
	$self->library_id($data->{'library_id'});
	$self->ssid($data->{'ssid'});
	$self->seq_status($data->{'seq_status'});
        $self->changed($data->{'changed'});
	$self->dirty(0);    # unset the dirty flag
    }
    else{
	die(sprintf('Cannot retrieve request: %s', $DBI::errstr));
    }

    return $self;
}


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
        my %allowed = map {$_ => 1} @{VRTrack::Utils::list_enum_vals($self->{_dbh},'request','seq_status')};
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


=head2 update

  Arg [1]    : None
  Example    : $request->update();
  Description: Update a request whose properties you have changed.  If properties haven't changed (i.e. dirty flag is unset) do nothing.  
	       Changes the changed datestamp to now() on the mysql server (i.e. you don't have to set changed yourself, and indeed if you do, it will be overridden).
  Returntype : 1 if successful, otherwise undef.

=cut

sub update {
    die "NOT YET IMPLEMENTED\n";
    my ($self) = @_;
    my $success = undef;
    if ($self->dirty){
	my $dbh = $self->{_dbh};
	my $save_re = $dbh->{RaiseError};
	my $save_pe = $dbh->{PrintError};
	my $save_ac = $dbh->{AutoCommit};
	$dbh->{RaiseError} = 1; # raise exception if an error occurs
	$dbh->{PrintError} = 0; # don't print an error message
	$dbh->{AutoCommit} = 0; # disable auto-commit

	eval {
	    # Need to unset 'latest' flag on current latest request and add
	    # the new request details with the latest flag set
	    my $updsql = qq[UPDATE request SET latest=false WHERE request_id = ? and latest=true];
	    
	    my $addsql = qq[INSERT INTO request (request_id, library_id, read_len, seq_status, seq_status, changed, latest) 
			    VALUES (?,?,?,?,?,now(),true)];
	    $dbh->do ($updsql, undef,$self->id);
	    $dbh->do ($addsql, undef,$self->id, $self->library_id, $self->read_len, $self->seq_status, $self->seq_status );
	    $dbh->commit ( );
	};

	if ($@) {
	    warn "Transaction failed, rolling back. Error was:\n$@\n";
	    # roll back within eval to prevent rollback
	    # failure from terminating the script
	    eval { $dbh->rollback ( ); };
	}
	else {
	    $success = 1;
	}

	# restore attributes to original state
	$dbh->{AutoCommit} = $save_ac;
	$dbh->{PrintError} = $save_pe;
	$dbh->{RaiseError} = $save_re;

    }

    return $success;
}


1;
