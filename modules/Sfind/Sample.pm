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

=cut

use Moose;
use namespace::autoclean;
use Sfind::Library;
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


has 'library_request_ids'=> (
    is          => 'ro',
    isa         => 'ArrayRef[Int]',
    lazy        => 1,
    builder     => '_get_library_request_ids',
);

has 'library_requests'=> (
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


 
=head2 library_requests

  Arg [1]    : library name from sequencescape
  Example    : my $librequests = $sample->library_requests();
  Description: Returns a ref to an array of the library_request objects that are associated with this sample.
  Returntype : ref to array of Sfind::Library_Request objects

=cut

sub _get_library_requests {
    my ($self) = @_;

    my @library_requests;
    foreach my $id (@{$self->library_request_ids()}){
        my $obj = Sfind::Library_Request->new($self->{_dbh},$id);
        push @library_requests, $obj; 
    }
    @library_requests = sort {$a <=> $b} @library_requests;

    return \@library_requests;
}
 
 
=head2 library_request_ids
 
   Arg [1]    : None
   Example    : my $library_request_ids = $sample->library_request_ids();
   Description: Returns a ref to an array of the library request IDs that are associated with this sample
   Returntype : ref to array of integer library request IDs
 
=cut
 
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


# =head2 id
# 
#   Arg [1]    : id (optional)
#   Example    : my $id = $samp->id();
# 	       $samp->id('104');
#   Description: Get/Set for ID of a sample
#   Returntype : SequenceScape ID (usu. integer)
# 
# =cut
# 
# sub id {
#     my ($self,$id) = @_;
#     if ($id){
# 	$self->{'id'} = $id;
#     }
#     return $self->{'id'};
# }
# 
# 
# =head2 study_id
# 
#   Arg [1]    : study_id (optional)
#   Example    : my $study_id = $samp->study_id();
# 	       $samp->study_id('104');
#   Description: Get/Set for study ID of a sample
#   Returntype : SequenceScape ID (usu. integer)
# 
# =cut
# 
# sub study_id {
#     my ($self,$study_id) = @_;
#     if ($study_id){
# 	$self->{'study_id'} = $study_id;
#     }
#     return $self->{'study_id'};
# }
# 
# 
# =head2 name
# 
#   Arg [1]    : name (optional)
#   Example    : my $name = $samp->name();
# 	       $samp->name('104');
#   Description: Get/Set for sample name
#   Returntype : string
# 
# =cut
# 
# sub name {
#     my ($self,$name) = @_;
#     if ($name){
# 	$self->{'name'} = $name;
#     }
#     return $self->{'name'};
# }
# 
# =head2 get_organism_name
# 
#   Arg [1]    : None
#   Example    : my $organism = $sample->get_organism_name();
#   Description: retrieve organism name from sample_common_name key for given sample ID
#   Returntype : string
# 
# =cut
# 
# sub get_organism_name {
#     my ($self) = @_;
#     my $sql = qq[select value from property_information where `key` like "sample_common_name" and property_information.obj_id=?];
#     my $org = $self->{_dbh}->selectrow_hashref($sql, undef, ($self->id));
#     unless ($org){
# 	warn "No organism ", $self->id,"\n";
# 	return undef;
#     }
#     my $orgname = $org->{value};
#     $orgname =~ s/^\s+//; #remove leading spaces
#     $orgname =~ s/\s+$//; #remove trailing spaces
#     return $orgname;
# }
# 
# 
# 
# 
# 
# =head2 get_strain_name
# 
#   Arg [1]    : None
#   Example    : my $strain = $sample->get_strain_name();
#   Description: retrieve strain information from given sample ID
#   Returntype : string
# 
# =cut
# 
# sub get_strain_name {
#     my ($self) = @_;
#     my $sql = qq[select value from property_information where `key` like "sample_strain_att" and property_information.obj_id=?];
#     my $strain = $self->{_dbh}->selectrow_hashref($sql, undef, ($self->id));
#     unless ($strain){
# 	warn "No strain ", $self->id,"\n";
# 	return undef;
#     }
#     my $strain_name = $strain->{value};
#     return $strain_name;
# }
# 
# =head2 get_accession
# 
#   Arg [1]    : None
#   Example    : my $acc = $sample->get_accession();
#   Description: retrieve EBI accession number from given sample ID
#   Returntype : string
# 
# =cut
# 
# sub get_accession {
#     my ($self) = @_;
#     my $sql = qq[select value from property_information where `key` like "sample_ebi_accession_number" and property_information.obj_id=?];
#     my $acc = $self->{_dbh}->selectrow_hashref($sql, undef, ($self->id));
#     unless ($acc){
# 	warn "No accession ", $self->id,"\n";
# 	return undef;
#     }
#     my $acc_name = $acc->{value};
#     return $acc_name;
# }
# 
# =head2 get_public_name
# 
#   Arg [1]    : None
#   Example    : my $pubname = $sample->get_public_name();
#   Description: retrieve sample public name from given sample ID
#   Returntype : string
# 
# =cut
# 
# sub get_public_name {
#     my ($self) = @_;
#     my $sql = qq[select value from property_information where `key` like "sample_public_name" and property_information.obj_id=?];
#     my $pub = $self->{_dbh}->selectrow_hashref($sql, undef, ($self->id));
#     unless ($pub){
# 	warn "No public name ", $self->id,"\n";
# 	return undef;
#     }
#     my $pub_name = $pub->{value};
#     return $pub_name;
# }
# 
__PACKAGE__->meta->make_immutable;
1;
