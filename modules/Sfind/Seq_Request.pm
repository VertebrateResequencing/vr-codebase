package Sfind::Seq_Request; 
=head1 NAME

Sfind::Request - Sequence Tracking Request object

=head1 SYNOPSIS
    my $seqrequest= Sfind::Seq_Request->new({dbh => $dbh, id => $req_id});

    my $id = $seqrequest->id();
    my $status = $seqrequest->status();

=head1 DESCRIPTION

An object describing the tracked properties of a sequencing request.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS


=head2 new

  Arg [1]    : hashref: dbh => database handle to seqtracking database
                        id  => sequencescape request id
  Example    : my $req = Sfind::Seq_Request->new({dbh=>$dbh, id=>$id};)
  Description: Returns Seq_Request object by request_id
  Returntype : Sfind::Seq_Request object


=head2 id

  Arg [1]    : None
  Example    : my $id = $seqrequest->id();
  Description: Returns ID of a request
  Returntype : SequenceScape ID (usu. integer)


=head2 created

  Arg [1]    : None
  Example    : my $created = $seqrequest->created();
  Description: Returns created timestamp
  Returntype : timestamp string


=head2 type

  Arg [1]    : None
  Example    : my $type = $seqrequest->type();
  Description: Returns type of request
  Returntype : sequencescape request type string


=head2 library_id

  Arg [1]    : None
  Example    : my $library_id = $seqrequest->library_id();
  Description: Returns library ID of a request
                This is the multiplex_tube_asset_id if the library is a
                multiplex library else it is the library tube asset id
  Returntype : SequenceScape ID integer


=head2 library_name

  Arg [1]    : None
  Example    : my $library_name = $seqrequest->library_name();
  Description: Returns library name of request
  Returntype : SequenceScape name


=head2 read_len

  Arg [1]    : None
  Example    : my $read_len = $seqrequest->read_len();
	       $seqrequest->read_len(54);
  Description: Returns request read_len
  Returntype : integer


=head2 status

  Arg [1]    : None
  Example    : my $status = $seqrequest->status();
  Description: Returns request status
  Returntype : string


=head2 lane_ids

  Arg [1]    : None
  Example    : my $lane_ids = $seqrequest->lane_ids();
  Description: Returns a ref to an array of the file names that are associated with this request
  Returntype : ref to array of file names

=head2 lanes

  Arg [1]    : None
  Example    : my $lanes = $seqrequest->lanes();
  Description: Returns a ref to an array of the file objects that are associated with this request.
  Returntype : ref to array of Sfind::File objects

=head2 sequenced_bases

  Arg [1]    : None
  Example    : my $tot_bp = $lib->sequenced_bases();
  Description: the total number of sequenced bases on this library.
		This is the sum of the bases from the fastq files associated
		with this library in NPG.
  Returntype : integer

=cut

use Moose;
use namespace::autoclean;
use Sfind::Types qw(MysqlDateTime);
use Sfind::Lane;

has '_dbh'  => (
    is          => 'ro',
    isa         => 'DBI::db',
    required    => 1,
    init_arg    => 'dbh',
);

has 'id'    => (
    is          => 'ro',
    isa         => 'Int',
    required    => 1,
);

has 'uuid'    => (
    is          => 'ro',
    isa         => 'Str',
    required    => 1,
);

has 'type' => (
    is          => 'ro',
    isa         => 'Str',
    init_arg    => 'request_type',
);

has 'library_id'    => (
    is          => 'ro',
    isa         => 'Int',
    init_arg    => 'source_asset_internal_id',
);

has 'library_name'    => (
    is          => 'ro',
    isa         => 'Str',
    init_arg    => 'source_asset_name',
);

has 'read_len'    => (
    is          => 'ro',
    isa         => 'Int',
    init_arg    => 'read_length',
);

has 'status'    => (
    is          => 'ro',
    isa         => 'Str',
    init_arg    => 'state',
);

has 'lane_ids'=> (
    is          => 'ro',
    isa         => 'ArrayRef[Int]',
    lazy        => 1,
    builder     => '_get_lane_ids',
);

has 'lanes'=> (
    is          => 'ro',
    isa         => 'ArrayRef[Sfind::Lane]',
    lazy        => 1,
    builder     => '_get_lanes',
);

has 'created' => (
    is          => 'ro',
    isa         => MysqlDateTime,
    coerce      => 1,   # accept mysql dates
);

# Populate the parameters from the database
around BUILDARGS => sub {
    my $orig  = shift;
    my $class = shift;
    
    my $argref = $class->$orig(@_);

    die "Need to call with a seq_request id" unless $argref->{id};
    my $sql = qq[select * from requests where internal_id = ? and is_current=1];
    my $id_ref = $argref->{dbh}->selectrow_hashref($sql, undef, ($argref->{id}));
    if ($id_ref){
        foreach my $field(keys %$id_ref){
            $argref->{$field} = $id_ref->{$field};
        }
    };
    return $argref;
};


###############################################################################
# BUILDERS
###############################################################################


sub _get_lanes {
    my ($self) = @_;
    my @lanes;
    foreach my $id (@{$self->lane_ids()}){
        my $obj = Sfind::Lane->new({dbh=>$self->{_dbh},id=>$id});
        push @lanes, $obj;
    }

    return \@lanes
}


sub _get_lane_ids {
    my ($self) = @_;
    my $sql = qq[select id_npg_information 
                from npg_information 
                where request_id = ? and
                (id_run_pair=0 or id_run_pair is null);];
    my @lane_ids;
    my $sth = $self->{_dbh}->prepare($sql);

    $sth->execute($self->id);
    foreach(@{$sth->fetchall_arrayref()}){
        push @lane_ids, $_->[0];
    }
    @lane_ids = sort {$a <=> $b} @lane_ids;
 
    return \@lane_ids;
}


###############################################################################
# Additional Methods
###############################################################################

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
