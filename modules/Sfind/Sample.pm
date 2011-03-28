package Sfind::Sample;
=head1 NAME

Sfind::Sample - Sequence Tracking Sample object

=head1 SYNOPSIS
    my $samp = Sfind::Sample->new({dbh => $dbh, 
                                    id => $sample_id,
                                study_id=> $study_id});

    #get arrayref of library objects in a sample
    my $libs = $sample->libraries();
    
    my $id = $sample->id();
    my $name = $sample->name();

=head1 DESCRIPTION

An object describing the tracked properties of a sample.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=head2 new

  Arg [1]    : hashref: dbh => database handle to seqtracking database
                        id  => sample id
                        study_id  => study id
  Example    : my $sample = Sfind::Sample->new({dbh=>$dbh, id=>$id};)
  Description: Returns Sample object by sample id
  Returntype : Sfind::Sample object


=head2 id

  Arg [1]    : none
  Example    : my $id = $samp->id();
  Description: Retrieve ID of a sample
  Returntype : Integer


=head2 study_id

  Arg [1]    : none
  Example    : my $study_id = $samp->study_id();
  Description: Retrieve study ID of a sample
  Returntype : Integer


=head2 name

  Arg [1]    : none
  Example    : my $name = $samp->name();
  Description: Retrieve sample name
  Returntype : string

 
=head2 uuid

  Arg [1]    : none
  Example    : my $uuid = $samp->uuid();
  Description: Retrieve sample uuid
  Returntype : string


=head2 description

  Arg [1]    : none
  Example    : my $description = $samp->description();
  Description: Retrieve sample description
  Returntype : string


=head2 organism

  Arg [1]    : none
  Example    : my $organism = $samp->organism();
  Description: Retrieve sample organism
  Returntype : string


=head2 common_name

  Arg [1]    : none
  Example    : my $common_name = $samp->common_name();
  Description: Retrieve sample common_name
  Returntype : string


=head2 taxon_id

  Arg [1]    : none
  Example    : my $taxon_id = $samp->taxon_id();
  Description: Retrieve sample taxon_id
  Returntype : int


=head2 accession

  Arg [1]    : none
  Example    : my $accession = $samp->accession();
  Description: Retrieve sample accession
  Returntype : string


=head2 gender

  Arg [1]    : none
  Example    : my $gender = $samp->gender();
  Description: Retrieve sample gender
  Returntype : string


=head2 sanger_sample_id

  Arg [1]    : none
  Example    : my $sanger_sample_id = $samp->sanger_sample_id();
  Description: Retrieve sample sanger_sample_id
  Returntype : string


=head2 supplier_name

  Arg [1]    : none
  Example    : my $supplier_name = $samp->supplier_name();
  Description: Retrieve sample supplier_name
  Returntype : string


=head2 created

  Arg [1]    : none
  Example    : my $sample_when = $sample->created();
  Description: Retrieves sample creation datetime
  Returntype : DateTime object


=head2 library_requests

  Arg [1]    : none
  Example    : my $librequests = $sample->library_requests();
  Description: Returns a ref to an array of the library_request objects that are associated with this sample.
  Returntype : ref to array of Sfind::Library_Request objects


=head2 library_request_ids

  Arg [1]    : none
  Example    : my $librequestids = $sample->library_request_ids();
  Description: Returns an arrayref of the library_request ids on this sample
  Returntype : ref to array of Sfind::Library_Request IDs

=cut

use Moose;
use Sfind::Types qw(MysqlDateTime);
use namespace::autoclean;
use Sfind::Library_Request;

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

has 'study_id'    => (
    is          => 'ro',
    isa         => 'Int',
    required    => 1,
);

has 'name'  => (
    is          => 'ro',
    isa         => 'Str',
    required    => 1,
);

has 'uuid'  => (
    is          => 'ro',
    isa         => 'Str',
    required    => 1,
);

has 'description'  => (
    is          => 'ro',
    isa         => 'Maybe[Str]',
);

has 'organism'  => (
    is          => 'ro',
    isa         => 'Maybe[Str]',
);

has 'common_name'	=> (
    is          => 'ro',
    isa         => 'Maybe[Str]',
);

has 'taxon_id'	=> (
    is          => 'ro',
    isa         => 'Maybe[Int]',
);

has 'accession'  => (
    is          => 'ro',
    isa         => 'Maybe[Str]',
    init_arg    => 'accession_number',
);

has 'gender'	=> (
    is          => 'ro',
    isa         => 'Maybe[Str]',
);

has 'sanger_sample_id'	=> (
    is          => 'ro',
    isa         => 'Maybe[Str]',
);

has 'supplier_name'	=> (
    is          => 'ro',
    isa         => 'Maybe[Str]',
);


has 'created' => (
    is          => 'ro',
    isa         => MysqlDateTime,
    coerce      => 1,   # accept mysql dates
);
# Add these when fields appear from Andrew Page.
#has 'strain'=> (
#    is          => 'ro',
#    isa         => 'Maybe[Str]',
#);
#has 'public_name'=> (
#    is          => 'ro',
#    isa         => 'Maybe[Str]',
#);

has 'library_request_ids'=> (
    is          => 'ro',
    isa         => 'ArrayRef[Int]',
    lazy        => 1,
    builder     => '_get_library_request_ids',
);

has 'library_requests' => (
    is          => 'ro',
    isa         => 'ArrayRef[Sfind::Library_Request]',
    lazy        => 1,
    builder     => '_get_library_requests',
);

# Populate the parameters from the database
around BUILDARGS => sub {
    my $orig  = shift;
    my $class = shift;
    
    my $argref = $class->$orig(@_);

    die "Need to call with a sample id" unless $argref->{id};
    die "Need to call with a study id" unless $argref->{study_id};
    # note this retrieves samples limited to a particular study.
    my $sql = qq[select s.* from samples s, study_samples ss where ss.sample_internal_id = ? and ss.study_internal_id = ? and s.internal_id=ss.sample_internal_id and s.is_current=1 and ss.is_current=1];
    my $id_ref = $argref->{dbh}->selectrow_hashref($sql, undef, ($argref->{id},$argref->{study_id}));
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
sub _get_library_requests {
    my ($self) = @_;

    my @library_requests;
    foreach my $id (@{$self->library_request_ids()}){
        my $obj = Sfind::Library_Request->new({ dbh => $self->_dbh,
                                                id  => $id, 
                                                });
        push @library_requests, $obj; 
    }
    @library_requests = sort {$a <=> $b} @library_requests;

    return \@library_requests;
}
 
 
sub _get_library_request_ids {
    my ($self) = @_;
    my $sql = qq[select internal_id from requests where (request_type like '%library creation' or request_type like '%library preparation') and source_asset_sample_internal_id=? and study_internal_id=? and is_current=1];
    my @library_requests;
    my $sth = $self->{_dbh}->prepare($sql);
    $sth->execute($self->id, $self->study_id);
    foreach(@{$sth->fetchall_arrayref()}){
        push @library_requests, $_->[0] if $_->[0];
    }
    return \@library_requests;
}

###############################################################################
# Additional methods
###############################################################################



###############################################################################
# DEPRECATED CALLS
###############################################################################
sub get_organism_name {
    my ($self) = @_;
    warnings::warnif("deprecated",
    "get_organism_name is deprecated, use organism instead");
    return $self->organism();
}

sub get_strain_name {
    my ($self) = @_;
    warnings::warnif("deprecated",
    "get_strain_name is deprecated, use strain instead");
    return $self->strain();
}

sub get_accession {
    my ($self) = @_;
    warnings::warnif("deprecated",
    "get_accession is deprecated, use accession instead");
    return $self->accession();
}

sub get_public_name {
    my ($self) = @_;
    warnings::warnif("deprecated",
    "get_public_name is deprecated, use public_name instead");
    return $self->public_name();
}

__PACKAGE__->meta->make_immutable;
1;
