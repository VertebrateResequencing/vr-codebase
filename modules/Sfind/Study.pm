package Sfind::Study;
=head1 NAME

Sfind::Study - Sequence Tracking Study object

=head1 SYNOPSIS
    my $study = Sfind::Study->new({dbh => $dbh, id=>$study_id});

    #get arrayref of sample objects in a study
    my $samples = $study->samples();
    
    my $id = $study->id();
    my $name = $study->name();
    my $uuid = $study->uuid();

=head1 DESCRIPTION

An object describing the tracked properties of a study.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=head2 new

  Arg [1]    : hashref: dbh => database handle to seqtracking database
                        id  => study id
  Example    : my $study = Sfind::Study->new({dbh=>$dbh, id=>$id};)
  Description: Returns Study object by study id
  Returntype : Sfind::Study object


=head2 id

  Arg [1]    : none
  Example    : my $id = $study->id();
  Description: Return ID of a study
  Returntype : SequenceScape ID (usu. integer)


=head2 name

  Arg [1]    : none
  Example    : my $name = $study->name();
  Description: Return study name
  Returntype : string


=head2 uuid

  Arg [1]    : none
  Example    : my $uuid = $study->uuid();
  Description: Return study uuid
  Returntype : string


=head2 description

  Arg [1]    : none
  Example    : my $description = $study->description();
  Description: Return study description
  Returntype : string


=head2 abstract

  Arg [1]    : none
  Example    : my $abstract = $study->abstract();
  Description: Return study abstract
  Returntype : string


=head2 sponsor

  Arg [1]    : none
  Example    : my $sponsor = $study->sponsor();
  Description: Return study Faculty/SAC sponsor
  Returntype : string


=head2 accession

  Arg [1]    : none
  Example    : my $accession = $study->accession();
  Description: Return study SRA/ERA accession
  Returntype : string


=head2 study_type

  Arg [1]    : none
  Example    : my $study_type = $study->study_type();
  Description: Return study type
  Returntype : string


=head2 abbreviation

  Arg [1]    : none
  Example    : my $abbreviation = $study->abbreviation();
  Description: Return study abbreviation.
  Returntype : string


=head2 ref_genome

  Arg [1]    : none
  Example    : my $ref_genome = $study->ref_genome();
  Description: Return reference genome to align this study to
  Returntype : string


=head2 ethically_approved

  Arg [1]    : none
  Example    : my $ethically_approved = $study->ethically_approved();
  Description: Return ethical approval status
  Returntype : string


=head2 sample_ids

  Arg [1]    : none
  Example    : my $sample_ids = $study->sample_ids();
  Description: Returns a ref to an array of the sample IDs that are associated with this study
  Returntype : arrayref of sorted integer sample IDs


=head2 samples

  Arg [1]    : none
  Example    : my $samples = $study->samples();
  Description: Returns a ref to an array of the sample objects that are associated with this study
  Returntype : ref to array of Sfind::Sample objects


=head2 contains_human_dna

  Arg [1]    : none
  Example    : my $chd = $study->contains_human_dna();
  Description: Returns the value in the 'contains_human_dna' field for this study
  Returntype : bool


=head2 contaminated_with_human_dna

  Arg [1]    : none
  Example    : my $chd = $study->contaminated_with_human_dna();
  Description: Returns the value in the 'contaminated_human_dna' field for this study
  Returntype : bool


=head2 visibility

  Arg [1]    : none
  Example    : my $vis = $study->visibility();
  Description: Public archive visibility of this study
  Returntype : String


=head2 title

  Arg [1]    : none
  Example    : my $title = $study->title();
  Description: Public archive title of this study
  Returntype : String


=head2 sra_project_id

  Arg [1]    : none
  Example    : my $sra_id = $study->sra_project_id();
  Description: Public archive project ID for this study
  Returntype : String


=head2 organism_names

  Arg [1]    : none
  Example    : my $orghash = $study->organism_names();
  Description: Convenient method to get all the organism names associated with
                this study. The hash will also give the number of samples for 
                a particular organism (if needed later on)
  Returntype : hashref of strings


=head2 created

  Arg [1]    : none
  Example    : my $study_when = $study->created();
  Description: Retrieves study creation datetime
  Returntype : DateTime object


=head2 state

  Arg [1]    : none
  Example    : my $study_state = $study->state();
  Description: Retrieves study state
  Returntype : string

=cut

use Moose;
use Sfind::Types qw(MysqlDateTime YesNoBool);
use Sfind::Sample;
use namespace::autoclean;


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
    isa         => 'Str',
);

has 'abstract'  => (
    is          => 'ro',
    isa         => 'Maybe[Str]',
);

has 'sponsor'  => (
    is          => 'ro',
    isa         => 'Str',
    init_arg    => 'faculty_sponsor',
);

has 'accession'  => (
    is          => 'ro',
    isa         => 'Maybe[Str]',
    init_arg    => 'accession_number',
);

has 'abbreviation'  => (
    is          => 'ro',
    isa         => 'Str',
);

has 'study_type'    => (
    is          => 'ro',
    isa         => 'Str',
);

has 'ref_genome'=> (
    is          => 'ro',
    isa         => 'Maybe[Str]',
    init_arg    => 'reference_genome',
);

has 'ethically_approved'=> (
    is          => 'ro',
    isa         => 'Bool',
);

has 'created' => (
    is          => 'ro',
    isa         => MysqlDateTime,
    coerce      => 1,   # accept mysql dates
);

has 'contains_human_dna'=> (
    is          => 'ro',
    isa         => YesNoBool,
    coerce      => 1,   # accept yes/no
);

has 'contaminated_with_human_dna'=> (
    is          => 'ro',
    isa         => YesNoBool,
    coerce      => 1,   # accept yes/no
    init_arg    => 'contaminated_human_dna',
);

has 'visibility'=> (
    is          => 'ro',
    isa         => 'Str',
    init_arg    => 'study_visibility',
);

has 'title'=> (
    is          => 'ro',
    isa         => 'Maybe[Str]',
    init_arg    => 'study_title',
);

has 'ena_project_id'=> (
    is          => 'ro',
    isa         => 'Maybe[Str]',
);

has 'state'=> (
    is          => 'ro',
    isa         => 'Str',
);

has 'sample_ids'=> (
    is          => 'ro',
    isa         => 'ArrayRef[Int]',
    lazy        => 1,
    builder     => '_get_sample_ids',
);

has 'samples'=> (
    is          => 'ro',
    isa         => 'ArrayRef[Sfind::Sample]',
    lazy        => 1,
    builder     => '_get_samples',
);

has 'organism_names'=> (
    is          => 'ro',
    isa         => 'HashRef',
    lazy        => 1,
    builder     => '_get_organism_names',
);

# Populate the parameters from the database
around BUILDARGS => sub {
    my $orig  = shift;
    my $class = shift;
    
    my $argref = $class->$orig(@_);
    die "Need to call with a study id" unless $argref->{id};

    my $sql = qq[select * from current_studies where internal_id=? ];
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
# builder to retrieve samples
sub _get_samples {
    my ($self) = @_;
    my @samples;
    foreach my $id (@{$self->sample_ids()}){
        if($id && $id=~/^\d/){
            my $obj = $self->get_sample_by_id($id);
            push @samples, $obj if $obj;
        }
    }
    return \@samples;
}



# builder to retrieve sample ids
sub _get_sample_ids {
    my ($self) = @_;
    my $sql = qq[select distinct(ss.sample_internal_id) from current_study_samples ss, current_samples s where ss.study_internal_id=? and ss.sample_internal_id = s.internal_id order by ss.sample_internal_id];
    my @samples;
    my $sth = $self->_dbh->prepare($sql);

    $sth->execute($self->id);
    foreach(@{$sth->fetchall_arrayref()}){
        push @samples, $_->[0];
    }
     
    return \@samples;
}


# builder for orgs
sub _get_organism_names {
    my ($self) = @_;
    my %organisms;

    my $sql = qq[select common_name, count(*) from current_samples cs, current_study_samples css where css.study_internal_id=? and css.sample_internal_id = cs.internal_id group by cs.common_name];
    my $sth = $self->_dbh->prepare($sql);

    $sth->execute($self->id);
    foreach(@{$sth->fetchall_arrayref()}){
        next unless $_->[0];
        $organisms{$_->[0]} = $_->[1];
    }
    return \%organisms;
}


###############################################################################
# Additional methods
###############################################################################

=head2 get_sample_by_id

  Arg [1]    : none
  Example    : my $sample = $study->get_sample_by_id(1154);
  Description: retrieve sample object by sequencescape id
  Returntype : Sfind::Sample object

=cut

sub get_sample_by_id {
    my ($self, $id) = @_;
    my $obj = Sfind::Sample->new({      dbh => $self->_dbh,
                                        id  => $id, 
                                  study_id  => $self->id});
    return $obj;
}


=head2 get_sample_by_name

  Arg [1]    : none
  Example    : my $sample = $study->get_sample_by_name('NA12878')
  Description: retrieve sample object by sequencescape name
  Returntype : Sfind::Sample object

=cut

sub get_sample_by_name {
    my ($self, $name) = @_;
    my $sql = qq[select css.sample_internal_id from current_samples cs, current_study_samples css where cs.internal_id = css.sample_internal_id and cs.name=? and css.study_internal_id=?];
    my $id_ref = $self->_dbh->selectrow_hashref($sql, undef, ($name, $self->id));
    unless ($id_ref){
	warn "No sample with name $name\n";
	return undef;
    }

    my $id = $id_ref->{sample_internal_id};
    return $self->get_sample_by_id($id);
}


###############################################################################
# DEPRECATED CALLS
###############################################################################
sub get_accession {
    my ($self) = @_;
    warnings::warnif("deprecated",
    "get_accession is deprecated, use accession instead");
    return $self->accession();
}


sub get_SAC_sponsor {
    my ($self) = @_;
    warnings::warnif("deprecated",
    "get_SAC_sponsor is deprecated, use sponsor instead");
    return $self->sponsor();
}


sub get_study_description {
    my ($self) = @_;
    warnings::warnif("deprecated",
    "get_study_description is deprecated, use description instead");
    return $self->description();
}


sub get_study_type {
    my ($self) = @_;
    warnings::warnif("deprecated",
    "get_study_type is deprecated, use study_type instead");
    return $self->study_type();
}


sub get_SRA_study_id {
    my ($self) = @_;
    warnings::warnif("deprecated",
    "get_SRA_study_id is deprecated, use ena_project_id instead");
    return $self->sra_project_id;
}


sub get_study_title {
    my ($self) = @_;
    warnings::warnif("deprecated",
    "get_study_title is deprecated, use title instead");
    return $self->title;
}


sub get_study_abstract {
    my ($self) = @_;
    warnings::warnif("deprecated",
    "get_study_abstract is deprecated, use abstract instead");
    return $self->abstract;
}


sub get_study_visibility {
    my ($self) = @_;
    warnings::warnif("deprecated",
    "get_study_visibility is deprecated, use visibility instead");
    return $self->visibility;
}

sub get_organism_names {
    my ($self) = @_;
    warnings::warnif("deprecated",
    "get_organism_names is deprecated, use organism_names instead");
    return %{$self->organism_names};
}

__PACKAGE__->meta->make_immutable;
1;
