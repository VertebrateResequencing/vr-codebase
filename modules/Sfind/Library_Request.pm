package Sfind::Library_Request; 
=head1 NAME

Sfind::Library_Request - Sequence Tracking Library_Request object

=head1 SYNOPSIS
    my $librequest= Sfind::Library_Request->new({dbh => $dbh, id => $request_id});

    my $id = $librequest->id();
    my $status = $librequest->status();

=head1 DESCRIPTION

An object describing the tracked properties of a library request.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=head2 new

  Arg [1]    : hashref: dbh => database handle to seqtracking database
                        id  => library_request_id
  Example    : my $librequest= Sfind::Library_Request->new({dbh=>$dbh, id=>$id};)
  Description: Returns Library_Request object by lib_request_id
  Returntype : Sfind::Library_Request object

=head2 id

  Arg [1]    : none
  Example    : my $id = $request->id();
  Description: Retrieve ID of a request
  Returntype : Integer

=head2 type

  Arg [1]    : none
  Example    : my $type = $request->type();
  Description: Retrieve type of request, e.g "Paired end sequencing"
  Returntype : String

=head2 status

  Arg [1]    : none
  Example    : my $status = $request->status();
  Description: Retrieve request status
  Returntype : string


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


=head2 created

  Arg [1]    : none
  Example    : my $libreq_when = $request->created();
  Description: Retrieves request creation datetime
  Returntype : DateTime object


=head2 libraries

  Arg [1]    : none
  Example    : my $libraries = $library_request->libraries();
  Description: Returns a ref to an array of the library objects that are associated with this library_request.
  Returntype : ref to array of Sfind::Library objects


=head2 library_ids

  Arg [1]    : none
  Example    : my $library_ids = $library_request->library_ids();
  Description: Returns a ref to an array of the library IDs that are associated with this library_request
  Returntype : ref to array of integer library IDs

=cut

use Moose;
use namespace::autoclean;
use Sfind::Types qw(MysqlDateTime);
use Sfind::Library;
use Sfind::Well_Library;

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

has 'type'  => (
    is          => 'ro',
    isa         => 'Str',
    init_arg    => 'request_type',
);

has 'fragment_size_from'  => (
    is          => 'ro',
    isa         => 'Maybe[Int]',
);

has 'fragment_size_to'  => (
    is          => 'ro',
    isa         => 'Maybe[Int]',
);

has 'status'  => (
    is          => 'ro',
    isa         => 'Str',
    init_arg    => 'state',
);

has 'library_ids'=> (
    is          => 'ro',
    isa         => 'ArrayRef[Int]',
    lazy        => 1,
    builder     => '_get_library_ids',
);

has 'libraries'=> (
    is          => 'ro',
    isa         => 'ArrayRef[Sfind::Library]',
    lazy        => 1,
    builder     => '_get_libraries',
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

    die "Need to call with a library_request id" unless $argref->{id};
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

sub _get_libraries { 
    my ($self) = @_;
    my @libraries;
    foreach my $id (@{$self->library_ids()}){
        my $obj;
        if ($self->type eq "Pulldown Multiplex Library Preparation"){
            $obj = Sfind::Well_Library->new({dbh=>$self->{_dbh},id=>$id});
        }
        else {
            $obj = Sfind::Library->new({dbh=>$self->{_dbh},id=>$id});
        }
        push @libraries, $obj; 
    }

    return \@libraries;
}


sub _get_library_ids {
    my ($self) = @_;

    # If this request is "normal" then the library ids will be library_tube
    # asset ids If this is a cherrypick pulldown request, then the library is
    # in a well on a plate so the warehouse queries are different and the
    # library_ids aren't for library_tubes but for wells
    my @lib_ids;

    if ($self->type eq "Pulldown Multiplex Library Preparation"){
        my $sql= qq[select descendant_internal_id from asset_links, requests where requests.source_asset_internal_id = asset_links.ancestor_internal_id and descendant_type="wells" and requests.internal_id=? and requests.is_current=1 and asset_links.is_current=1];
        my $sth = $self->{_dbh}->prepare($sql);

        $sth->execute($self->id);
        foreach(@{$sth->fetchall_arrayref()}){
            if ($_->[0]){
                push @lib_ids, $_->[0];
            }
        }
    }
    else {
        # the target of the request should be either a standard library tube
        # for a non-multiplex request, or the indexed library tube that will be
        # pooled for a multiplexed request

        my $sql= qq[select target_asset_internal_id, target_asset_type from requests where internal_id=? and is_current=1];
          
        my $sth = $self->{_dbh}->prepare($sql);

        $sth->execute($self->id);
        foreach(@{$sth->fetchall_arrayref()}){
            if ($_->[0]){
                die "Unexpected target type ".$_->[1] unless $_->[1] eq 'library_tubes';
                push @lib_ids, $_->[0];
            }
        }
    }
    
    return \@lib_ids;
}


__PACKAGE__->meta->make_immutable;
1;

