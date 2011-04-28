package Sfind::Library; 
=head1 NAME

Sfind::Library - Sequence Tracking Library object

=head1 SYNOPSIS
    my $lib = Sfind::Library->new({dbh => $dbh, id => $library_id});

    #get arrayref of requests on a library
    my $reqs = $library->requests();

    my $id = $library->id();
    my $name = $library->name();


=head1 DESCRIPTION

An object describing the tracked properties of a library.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=head2 new

  Arg [1]    : hashref: dbh => database handle to seqtracking database
                        id  => librarytube asset id 
  Example    : my $lib = Sfind::Library->new({dbh=>$dbh, id=>$id};)
  Description: Returns Library object by library_id
  Returntype : Sfind::Library object


=head2 seq_requests

  Arg [1]    : None
  Example    : my $seq_requests = $library->seq_requests();
  Description: Returns a ref to an array of the seq_request objects that are associated with this library.
  Returntype : ref to array of Sfind::Seq_Request objects


=head2 seq_request_ids

  Arg [1]    : None
  Example    : my $seq_request_ids = $library->seq_request_ids();
  Description: Returns a ref to an array of the seq_request IDs that are associated with this library, including those on multiplexes that this library is in
  Returntype : ref to array of integer seq request IDs


=head2 id

  Arg [1]    : None
  Example    : my $id = $lib->id();
  Description: Retrieve ID of a library
  Returntype : Integer asset_id


=head2 sample_id

  Arg [1]    : None
  Example    : my $sample_id = $lib->sample_id();
  Description: Retrieve sample ID of a library
  Returntype : Integer sample id
  

=head2 name

  Arg [1]    : None
  Example    : my $name = $lib->name();
  Description: Retrieve library name
  Returntype : string


=head2 prep_status

  Arg [1]    : None
  Example    : my $prep_status = $lib->prep_status();
  Description: Retrieve library state
  Returntype : string


=head2 get_seq_request_by_id

  Arg [1]    : None
  Example    : my $seqrequest = $lib->get_seq_request_by_id(7447);
  Description: retrieve seq_request object by sequencescape id
  Returntype : Sfind::Seq_Request object


=head2 get_lane_by_id

  Arg [1]    : None
  Example    : my $lane = $lib->get_lane_by_id(7447);
  Description: retrieve NPG lane object by id
  Returntype : Sfind::Lane object


=head2 fragment_size

  Arg [1]    : None
  Example    : my ($frag_from, $frag_to) = @{$lib->fragment_size};
  Description: Fetches requested library fragment size, as a pair of from and to sizes
  Returntype : arrayref of [from, to] sizes or undefs


=head2 fragment_size_from

  Arg [1]    : none
  Example    : my $frag_from = $request->fragment_size_from();
  Description: Retrieve fragment size from on request
  Returntype : integer


=head2 fragment_size_to

  Arg [1]    : none
  Example    : my $frag_to = $request->fragment_size_to();
  Description: Retrieve fragment size to on request
  Returntype : integer


=head2 type

  Arg [1]    : None
  Example    : my $lib_type = $lib->type();
  Description: returns type of library, if in database.  e.g. 'Standard', 'No PCR'
  Returntype : string or undef


=head2 multiplex_pool_asset_ids

  Arg [1]    : None
  Description: Retrieve multiplex_tube_asset_ids
  Returntype : reference to a array


=head2 is_tagged

  Arg [1]    : None
  Example    : my $tag = $lib->is_tagged();
  Description: Retrieve whether Library is tagged.  The tag is a short sequence tag added to each Library molecule so that it can be sequenced in a multiplex and the resulting mix resolved by the tag sequence.
  Returntype : boolean


=head2 tag_id

  Arg [1]    : None
  Example    : my $tag_id = $lib->tag_id();
  Description: Retrieve tag ID of a library
  Returntype : int


=head2 tag_group_id

  Arg [1]    : None
  Example    : my $tag_group_id = $lib->tag_group_id();
  Description: Retrieve tag group ID.  Can have the same tag id in multiple tag groups.
  Returntype : int


=head2 tag_sequence

  Arg [1]    : None
  Example    : my $tag_sequence = $lib->tag_sequence();
  Description: Retrieve tag sequence.
  Returntype : string

=cut

use Moose;
use namespace::autoclean;
use Sfind::Types qw(MysqlDateTime);
use Sfind::Seq_Request;
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

has 'name'  => (
    is          => 'ro',
    isa         => 'Str',
);

has 'sample_id'    => (
    is          => 'ro',
    isa         => 'Int',
    init_arg    => 'sample_internal_id',
);

has 'prep_status'  => (
    is          => 'ro',
    isa         => 'Maybe[Str]',
    # should this be the state of the library_request instead?
    init_arg    => 'state',
);

has 'fragment_size_from'  => (
    is          => 'ro',
    isa         => 'Maybe[Int]',
    init_arg    => 'fragment_size_required_from',
);

has 'fragment_size_to'  => (
    is          => 'ro',
    isa         => 'Maybe[Int]',
    init_arg    => 'fragment_size_required_to',
);

has 'fragment_size'  => (
    is          => 'ro',
    isa         => 'ArrayRef[Int]',
    lazy        => 1,
    default => sub {
        my $self = shift;
        if ($self->fragment_size_from() && $self->fragment_size_to()) {
        	return [$self->fragment_size_from(), $self->fragment_size_to()];
        }
        return [0,0];	
    },
);

has 'type'  => (
    is          => 'ro',
    isa         => 'Str',
    init_arg    => 'library_type',
);

has 'is_tagged'  => (
    is          => 'ro',
    isa         => 'Bool',
    lazy        => 1,
    builder     => '_get_is_tagged',
);

has 'tag_id'  => (
    is          => 'ro',
    isa         => 'Maybe[Int]',
    init_arg    => 'tag_internal_id',
);

has 'tag_group_id'  => (
    is          => 'ro',
    isa         => 'Maybe[Int]',
    init_arg    => 'tag_group_internal_id',
);

has 'tag_sequence'  => (
    is          => 'ro',
    isa         => 'Maybe[Str]',
    init_arg    => 'expected_sequence',
);

has 'multiplex_pool_asset_ids'=> (
    is          => 'ro',
    isa         => 'ArrayRef[Int]',
    lazy        => 1,
    builder     => '_get_mplex_pool_ids',
);

has 'seq_request_ids'=> (
    is          => 'ro',
    isa         => 'ArrayRef[Int]',
    lazy        => 1,
    builder     => '_get_seq_req_ids',
);

has 'seq_requests'=> (
    is          => 'ro',
    isa         => 'ArrayRef[Sfind::Seq_Request]',
    lazy        => 1,
    builder     => '_get_seq_reqs',
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

    die "Need to call with a librarytube asset id" unless $argref->{id};
    my $sql = qq[select * from library_tubes where internal_id = ? and is_current=1];
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
sub _get_seq_reqs {
    my ($self) = @_;
    my @seq_requests;
    foreach my $id (@{$self->seq_request_ids()}){
        my $obj = $self->get_seq_request_by_id($id);
        push @seq_requests, $obj;
    }

    return \@seq_requests;
}


sub _get_seq_req_ids {
    my ($self) = @_;
    # There are three cases:
    # 1 non mplex lib, requests are on this lib as source_asset
    # 2 mplex lib, requests are on the mplexes this lib is in
    # 3 cherrypick mplex, handled in Well_Library::multiplex_pool_asset_ids

        
    # seq requests can be on this library or a multiplex that this library
    # is in.
    my @asset_ids = ($self->id);

    # case 2 first, then do case 1 & 2 at once
    if ($self->is_tagged){
        push @asset_ids,@{$self->multiplex_pool_asset_ids};
    }
        
    my $sql = qq[select distinct internal_id 
            from 
            requests where source_asset_internal_id=? 
            and request_type like '%sequencing'
            and is_current = 1
            ];

    my @seq_requests;
    my $sth = $self->_dbh->prepare($sql);
    foreach my $asset_id(@asset_ids){
        $sth->execute($asset_id);
        foreach(@{$sth->fetchall_arrayref()}){
            push @seq_requests, $_->[0];
        }
    }
    @seq_requests = sort {$a <=> $b} @seq_requests;

    return \@seq_requests;
}

sub _get_mplex_pool_ids{
    my ($self) = @_;
    my @mplex_ids;

    my $sql = qq[select descendant_internal_id as mplex_id
                from asset_links 
                where ancestor_type="library_tubes" 
                and ancestor_internal_id=?
                and descendant_type="multiplexed_library_tubes"
                and is_current=1];

    my $sth = $self->_dbh->prepare($sql);
    $sth->execute($self->id);
    foreach(@{$sth->fetchall_arrayref()}){
        push @mplex_ids, $_->[0];
    }

    @mplex_ids = sort {$a <=> $b} @mplex_ids;

    return \@mplex_ids;
}


sub _get_is_tagged {
    my ($self) = @_;
    my $tag = $self->tag_id ? 1 : 0;
    return $tag;
}

###############################################################################
# Additional Methods
###############################################################################

sub get_seq_request_by_id {
    my ($self, $id) = @_;
    my $obj = Sfind::Seq_Request->new({dbh=>$self->{_dbh},id=>$id});
    return $obj;
}


sub get_lane_by_id {
    die "NOT IMPLEMENTED\n";
    my ($self, $id) = @_;
    my $obj = Sfind::Lane->new($self->{_dbh},$id);
    return $obj;
}


1;
