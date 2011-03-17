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

=cut

use Moose;
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
);

has 'description'  => (
    is          => 'ro',
    isa         => 'Str',
);

has 'abstract'  => (
    is          => 'ro',
    isa         => 'Str',
);

has 'sponsor'  => (
    is          => 'ro',
    isa         => 'Str',
    init_arg    => 'faculty_sponsor',
);

has 'accession'  => (
    is          => 'ro',
    isa         => 'Str',
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
    isa         => 'Str',
    init_arg    => 'reference_genome',
);

has 'ethically_approved'=> (
    is          => 'ro',
    isa         => 'Bool',
);

# Add these when fields appear from Andrew Page.
# Might need coercion of yes/no to bool? 
#has 'contains_human_dna'=> (
#    is          => 'ro',
#    isa         => 'Bool',
#);

#has 'contaminated_with_human_dna'=> (
#    is          => 'ro',
#    isa         => 'Bool',
#    init_arg    => 'contaminated_human_dna',
#);

#has 'visibility'=> (
#    is          => 'ro',
#    isa         => 'Str',
#    init_arg    => 'sra_study_hold',
#);

#has 'title'=> (
#    is          => 'ro',
#    isa         => 'Str',
#    init_arg    => 'study_title',
#);

#has 'sra_project_id'=> (
#    is          => 'ro',
#    isa         => 'Str',
#);

has 'sample_ids'=> (
    is          => 'ro',
    isa         => 'ArrayRef[Int]',
    lazy        => 1,
    builder     => '_get_sample_ids',
);

has 'samples'=> (
    is          => 'ro',
    isa         => 'ArrayRef[Sfind::Study]',
    lazy        => 1,
    builder     => '_get_samples',
);


# Populate the parameters from the database
around BUILDARGS => sub {
    my $orig  = shift;
    my $class = shift;
    
    my $argref = $class->$orig(@_);

    my $sql = qq[select * from studies where internal_id=? and is_current=1];
    my $id_ref = $argref->{dbh}->selectrow_hashref($sql, undef, ($argref->{id}));
    if ($id_ref){
        foreach my $field(keys %$id_ref){
            $argref->{$field} = $id_ref->{$field};
        }
    };
    return $argref;
};


=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $study->id();
	       $study->id('104');
  Description: Return ID of a study
  Returntype : SequenceScape ID (usu. integer)


=head2 name

  Arg [1]    : name (optional)
  Example    : my $name = $study->name();
	       $study->name('1000Genomes-A1-CEU');
  Description: Return study name
  Returntype : string


=head2 uuid

  Arg [1]    : uuid (optional)
  Example    : my $uuid = $study->uuid();
	       $study->uuid('60f676b2-dd7f-11df-a1d2-00144f206e2e');
  Description: Return study uuid
  Returntype : string


=head2 description

  Arg [1]    : description (optional)
  Example    : my $description = $study->description();
	       $study->description($big_long_description_text);
  Description: Return study description
  Returntype : string


=head2 abstract

  Arg [1]    : abstract (optional)
  Example    : my $abstract = $study->abstract();
	       $study->abstract($big_long_abstract_text);
  Description: Return study abstract
  Returntype : string


=head2 sponsor

  Arg [1]    : sponsor (optional)
  Example    : my $sponsor = $study->sponsor();
	       $study->sponsor('Richard Durbin');
  Description: Return study Faculty/SAC sponsor
  Returntype : string


=head2 accession

  Arg [1]    : accession (optional)
  Example    : my $accession = $study->accession();
	       $study->accession('SRP000542');
  Description: Return study SRA/ERA accession
  Returntype : string


=head2 study_type

  Arg [1]    : study_type (optional)
  Example    : my $study_type = $study->study_type();
	       $study->study_type('Whole Genome Sequencing');
  Description: Return study type
  Returntype : string


=head2 abbreviation

  Arg [1]    : abbreviation (optional)
  Example    : my $abbreviation = $study->abbreviation();
	       $study->abbreviation('UK10K_MUIR');
  Description: Return study abbreviation.
  Returntype : string


=head2 ref_genome

  Arg [1]    : ref_genome (optional)
  Example    : my $ref_genome = $study->ref_genome();
	       $study->ref_genome('Homo_sapiens (NCBI36)');
  Description: Return reference genome to align this study to
  Returntype : string


=head2 ethically_approved

  Arg [1]    : ethically_approved (optional)
  Example    : my $ethically_approved = $study->ethically_approved();
	       $study->ethically_approved(0);
  Description: Return ethical approval status
  Returntype : string


=head2 sample_ids

  Arg [1]    : None
  Example    : my $sample_ids = $study->sample_ids();
  Description: Returns a ref to an array of the sample IDs that are associated with this study
  Returntype : arrayref of sorted integer sample IDs


=head2 samples

  Arg [1]    : None
  Example    : my $samples = $study->samples();
  Description: Returns a ref to an array of the sample objects that are associated with this study
  Returntype : ref to array of Sfind::Sample objects


=head2 contains_human_dna

  Arg [1]    : None
  Example    : my $chd = $study->contains_human_dna();
  Description: Returns the value in the 'contains_human_dna' field for this study
  Returntype : bool


=head2 contaminated_with_human_dna

  Arg [1]    : None
  Example    : my $chd = $study->contaminated_with_human_dna();
  Description: Returns the value in the 'contaminated_human_dna' field for this study
  Returntype : bool


=head2 visibility

  Arg [1]    : None
  Example    : my $vis = $study->visibility();
  Description: Public archive visibility of this study
  Returntype : String


=head2 title

  Arg [1]    : None
  Example    : my $title = $study->title();
  Description: Public archive title of this study
  Returntype : String


=head2 sra_project_id

  Arg [1]    : None
  Example    : my $sra_id = $study->sra_project_id();
  Description: Public archive project ID for this study
  Returntype : String
=cut


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
    my $sql = qq[select distinct(sample_internal_id) from study_samples  where study_internal_id=? and is_current=1 order by sample_internal_id];
    my @samples;
    my $sth = $self->_dbh->prepare($sql);

    $sth->execute($self->id);
    foreach(@{$sth->fetchall_arrayref()}){
        push @samples, $_->[0];
    }
     
    return \@samples;
}


=head2 get_sample_by_id

  Arg [1]    : sample id from sequencescape
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

  Arg [1]    : sample name from sequencescape
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


=head2 get_accession

  Arg [1]    : None
  Example    : my $acc = $study->get_accession();
  Description: retrieve EBI accession.  Deprecated - use $study->accession() instead.
  Returntype : string

=cut

sub get_accession {
    my ($self) = @_;
    return $self->accession();
}


=head2 get_SAC_sponsor

  Arg [1]    : None
  Example    : my $acc = $study->get_SAC_sponsor();
  Description: retrieve SAC sponsor.  Deprecated - use $study->sponsor() instead.
  Returntype : string

=cut

sub get_SAC_sponsor {
    my ($self) = @_;
    return $self->sponsor();
}


=head2 get_study_description

  Arg [1]    : None
  Example    : my $desc = $study->get_study_description();
  Description: retrieve study description.  Deprecated - use $study->description instead.
  Returntype : string

=cut

sub get_study_description {
    my ($self) = @_;
    return $self->description();
}



=head2 get_study_type

  Arg [1]    : None
  Example    : my $type = $study->get_study_type();
  Description: Returns the study_type.  Deprecated: use $study->study_type instead
  Returntype : string

=cut

sub get_study_type {
    my ($self) = @_;
    return $self->study_type();
}

=head2 get_SRA_study_id

  Arg [1]    : None
  Example    : my $sra_ids = $study->get_SRA_study_id();
  Description: Returns the SRA id of this study
                Deprecated: use $self->sra_project_id instead
  Returntype : string

=cut

sub get_SRA_study_id {
    my ($self) = @_;
    return $self->sra_project_id;
}

=head2 get_study_title

  Arg [1]    : None
  Example    : my $title = $study->get_study_title();
  Description: Returns the title of this study
                Deprecated: use $study->title instead.
  Returntype : string

=cut

sub get_study_title {
    my ($self) = @_;
    return $self->title;
}

=head2 get_study_abstract

  Arg [1]    : None
  Example    : my $abstract = $study->get_study_abstract();
  Description: Returns the abstract of this study.  
                Deprecated: use $study->abstract instead.
  Returntype : string

=cut

sub get_study_abstract {
    my ($self) = @_;
    return $self->abstract;
}


=head2 get_study_visibility

  Arg [1]    : None
  Example    : my $visibility = $study->get_study_visibility();
  Description: Returns the sra visibility for this study
                Deprecated, use $study->visibility instead
  Returntype : string

=cut

sub get_study_visibility {
    my ($self) = @_;
    return $self->visibility;
}


=head2 get_organism_names

  Arg [1]    : None
  Example    : my %orgs = $study->get_organism_name();
  Description: Convenient method to get all the organism names associated with 
               this study, undef if no organism. The hash will also give the
               number of samples for a particular organism (if needed later on)
  Returntype : hash of strings

=cut

sub get_organism_names {
    my ($self) = @_;
    my %organisms;
    foreach my $id (@{$self->sample_ids()}){
	    if($id=~/^\d/){
	        my $sample = Sfind::Sample->new($self->_dbh,$id, $self->id);
	        my $org_name = $sample->get_organism_name();
		if($org_name and $org_name!~/^\s*$/){ #If defined and not empty string
		    if (exists $organisms{$org_name}){
	        	$organisms{$org_name} = $organisms{$org_name}++;
		    }else{
	        	$organisms{$org_name} = 1;	
		    }
		}

	    }
    }
    return %organisms;
}


__PACKAGE__->meta->make_immutable;
1;
