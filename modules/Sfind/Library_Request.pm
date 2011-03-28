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
    isa         => 'Int',
);

has 'fragment_size_to'  => (
    is          => 'ro',
    isa         => 'Int',
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
        my $obj = Sfind::Library->new({dbh=>$self->{_dbh},id=>$id});
        push @libraries, $obj; 
    }
    @libraries = sort {$a <=> $b} @libraries;

    return \@libraries;
}


sub _get_library_ids {
    my ($self) = @_;

    # check if this is a multiplex library creation request
    # return the asset_ids as as library_ids
    die "NOT IMPLEMENTED";
  
    # select all library_tube asset ids associated with this request
    my $sql= qq[select target_asset_id from requests_new where request_id=?];
          
    my @libraries;
    my $sth = $self->{_dbh}->prepare($sql);

    $sth->execute($self->id);
    foreach(@{$sth->fetchall_arrayref()}){
        push @libraries, $_->[0] if $_->[0];
    }
    
    return \@libraries;
}


__PACKAGE__->meta->make_immutable;
1;

1;
