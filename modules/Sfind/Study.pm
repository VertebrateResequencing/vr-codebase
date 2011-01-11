package Sfind::Study;
=head1 NAME

Sfind::Study - Sequence Tracking Study object

=head1 SYNOPSIS
    my $study = Sfind::Study->new($dbh, $study_id);

    #get arrayref of sample objects in a study
    my $samples = $study->samples();
    
    my $id = $study->id();
    my $name = $study->name();

=head1 DESCRIPTION

An object describing the tracked properties of a study.

=head1 CONTACT

nds@sanger.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;
no warnings 'uninitialized';
use Sfind::Sample;
use DBI;


=head2 new

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : study id
  Example    : my $study = Sfind::Sfind->new($dbh, $id)
  Description: Returns Study object by study_id
  Returntype : Sfind::Study object

=cut

sub new {
    my ($class,$dbh, $id) = @_;    
    die "Need to call with a db handle and id" unless ($dbh && $id);
    my $self = {};
    bless ($self, $class);
    $self->{_dbh} = $dbh;
    my $sql = qq[select study_name from study_information where study_id=?];
    my $id_ref = $self->{_dbh}->selectrow_hashref($sql, undef, ($id));
    if ($id_ref){
	my $name = $id_ref->{study_name};
        $name =~ s/\s+$//;  # trim trailing whitespace
	$self->id($id);
	$self->name($name);
    }
    else {
	return undef;
    }

    return $self;
}



=head2 samples

  Arg [1]    : None
  Example    : my $samples = $study->samples();
  Description: Returns a ref to an array of the sample objects that are associated with this study
  Returntype : ref to array of Sfind::Sample objects

=cut

sub samples {
    my ($self) = @_;

    unless ($self->{'samples'}){
	my @samples;
    	foreach my $id (@{$self->sample_ids()}){

	    if($id && $id=~/^\d/){
	        my $obj = Sfind::Sample->new($self->{_dbh},$id, $self->id);
	        push @samples, $obj;
	    }
	}
	$self->{'samples'} = \@samples;
    }

    return $self->{'samples'};
}


=head2 sample_ids

  Arg [1]    : None
  Example    : my $sample_ids = $study->sample_ids();
  Description: Returns a ref to an array of the sample IDs that are associated with this study
  Returntype : ref to sorted array of integer sample IDs

=cut

sub sample_ids {
    my ($self) = @_;

    # jws 2011-01-05
    # Changed this subroutine
    # from:  get distinct sample_ids from requests table
    # to: get sample_ids from study_sample_reports table
    # this is cleaner, and gives the current set of samples associated
    # with a study, rather than any samples which may have been associated in
    # the past, and are presumably no longer wanted (e.g. sample moves)

    unless ($self->{'sample_ids'}){
	my $sql = qq[select distinct(sample_id) from study_sample_reports where study_id=? order by sample_id];
	my @samples;
	my $sth = $self->{_dbh}->prepare($sql);

	$sth->execute($self->id);
	foreach(@{$sth->fetchall_arrayref()}){
	    push @samples, $_->[0];
	}
	$self->{'sample_ids'} = \@samples;
    }
 
    return $self->{'sample_ids'};
}


=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $study->id();
	       $study->id('104');
  Description: Get/Set for ID of a study
  Returntype : SequenceScape ID (usu. integer)

=cut

sub id {
    my ($self,$id) = @_;
    if ($id){
	$self->{'id'} = $id;
    }
    return $self->{'id'};
}


=head2 name

  Arg [1]    : name (optional)
  Example    : my $name = $study->name();
	       $study->name('1000Genomes-A1-CEU');
  Description: Get/Set for study name
  Returntype : string

=cut

sub name {
    my ($self,$name) = @_;
    if ($name){
	$self->{'name'} = $name;
    }
    return $self->{'name'};
}


=head2 get_sample_by_id

  Arg [1]    : sample id from sequencescape
  Example    : my $sample = $study->get_sample_by_id(1154);
  Description: retrieve sample object by sequencescape id
  Returntype : Sfind::Sample object

=cut

sub get_sample_by_id {
    my ($self, $id) = @_;
    my $obj = Sfind::Sample->new($self->{_dbh},$id,$self->id);
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
    my $sql = qq[select distinct(sample_id) from requests where sample_name=? and study_id=?];
    my $id_ref = $self->{_dbh}->selectrow_hashref($sql, undef, ($name, $self->id));
    unless ($id_ref){
	warn "No sample with name $name\n";
	return undef;
    }

    my $id = $id_ref->{sample_id};
    return $self->get_sample_by_id($id);
}


=head2 get_accession

  Arg [1]    : None
  Example    : my $acc = $study->get_accession();
  Description: retrieve EBI accession number from given study ID
  Returntype : string

=cut

sub get_accession {
    my ($self) = @_;
    my $sql = qq[select value from study_information where `param` like 'study_ebi_accession_number' and study_information.study_id=?];
    my $acc = $self->{_dbh}->selectrow_hashref($sql, undef, ($self->id));
    unless ($acc){
	warn "No accession ", $self->id,"\n";
	return undef;
    }
    my $acc_name = $acc->{value};
    return $acc_name;
}

=head2 get_SAC_sponsor

  Arg [1]    : None
  Example    : my $acc = $study->get_SAC_sponsor();
  Description: retrieve SAC sponsor from given study ID
  Returntype : string

=cut

sub get_SAC_sponsor {
    my ($self) = @_;
    my $sql = qq[select value from study_information where `param` like "sac_sponsor" and study_information.study_id=?];
    my $acc = $self->{_dbh}->selectrow_hashref($sql, undef, ($self->id));
    unless ($acc){
	warn "No sponsor", $self->id,"\n";
	return undef;
    }
    my $acc_name = $acc->{value};
    return $acc_name;
}

=head2 get_study_description

  Arg [1]    : None
  Example    : my $desc = $study->get_study_description();
  Description: retrieve study description from given study ID
  Returntype : string

=cut

sub get_study_description {
    my ($self) = @_;
    my $sql = qq[select value from study_information where `param` like "study_description" and study_information.study_id=?];
    my $acc = $self->{_dbh}->selectrow_hashref($sql, undef, ($self->id));
    unless ($acc){
	warn "No description", $self->id,"\n";
	return undef;
    }
    my $acc_name = $acc->{value};
    return $acc_name;
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
	        my $sample = Sfind::Sample->new($self->{_dbh},$id, $self->id);
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

=head2 get_study_type

  Arg [1]    : None
  Example    : my $type = $study->get_study_type();
  Description: Returns the value for 'study_study_type' in the database for this study
  Returntype : string

=cut

sub get_study_type {
    my ($self) = @_;
    my $sql = qq[select value from study_information where `param` like "study_study_type" and study_information.study_id=?];
    my $type_ref = $self->{_dbh}->selectrow_hashref($sql, undef, ($self->id));
    unless ($type_ref){
	warn "No study type", $self->id,"\n";
	return undef;
    }
    my $type = $type_ref->{value};
    return $type;
}

=head2 get_SRA_study_id

  Arg [1]    : None
  Example    : my $sra_ids = $study->get_SRA_study_id();
  Description: Returns the SRA id of this study
  Returntype : string

=cut

sub get_SRA_study_id {
    my ($self) = @_;
    my $sql = qq[select value from study_information where `param` like "study_study_id" and study_information.study_id=?];
    my $id_ref = $self->{_dbh}->selectrow_hashref($sql, undef, ($self->id));
    unless ($id_ref){
	warn "No SRA id", $self->id,"\n";
	return undef;
    }
    my $id = $id_ref->{value};
    return $id;
}

=head2 get_study_title

  Arg [1]    : None
  Example    : my $title = $study->get_study_title();
  Description: Returns the title of this study
  Returntype : string

=cut

sub get_study_title {
    my ($self) = @_;
    my $sql = qq[select value from study_information where `param` like "study_study_title" and study_information.study_id=?];
    my $title_ref = $self->{_dbh}->selectrow_hashref($sql, undef, ($self->id));
    unless ($title_ref){
	
	return undef;
    }
    my $title = $title_ref->{value};
    return $title;
}

=head2 get_study_abstract

  Arg [1]    : None
  Example    : my $abstract = $study->get_study_abstract();
  Description: Returns the abstract of this study
  Returntype : string

=cut

sub get_study_abstract {
    my ($self) = @_;
    my $sql = qq[select value from study_information where `param` like "study_abstract" and study_information.study_id=?];
    my $abstract_ref = $self->{_dbh}->selectrow_hashref($sql, undef, ($self->id));
    unless ($abstract_ref){
	return undef;
    }
    my $abstract = $abstract_ref->{value};
    return $abstract;
}


=head2 get_study_visibility

  Arg [1]    : None
  Example    : my $visibility = $study->get_study_visibility();
  Description: Returns the sra visibility for this study
  Returntype : string

=cut

sub get_study_visibility {
    my ($self) = @_;
    my $sql = qq[select value from study_information where `param` like "study_sra_hold" and study_information.study_id=?];
    my $vis_ref = $self->{_dbh}->selectrow_hashref($sql, undef, ($self->id));
    unless ($vis_ref){
	return undef;
    }
    my $visibility = $vis_ref->{value};
    return $visibility;
}



=head2 contains_human_dna

  Arg [1]    : None
  Example    : my $chd = $study->contains_human_dna();
  Description: Returns the value in the 'contains_human_dna' field for this study
  Returntype : string

=cut

sub contains_human_dna {
    my ($self) = @_;
    my $sql = qq[select value from study_information where `param` like "contains_human_dna" and study_information.study_id=?];
    my $chd_ref = $self->{_dbh}->selectrow_hashref($sql, undef, ($self->id));
    unless ($chd_ref){
	return undef;
    }
    my $chd = $chd_ref->{value};
    return $chd;
}


=head2 contaminated_with_human_dna

  Arg [1]    : None
  Example    : my $chd = $study->contaminated_with_human_dna();
  Description: Returns the value in the 'contaminated_human_dna' field for this study
  Returntype : string

=cut

sub contaminated_with_human_dna {
    my ($self) = @_;
    my $sql = qq[select value from study_information where `param` like "contaminated_human_dna" and study_information.study_id=?];
    my $chd_ref = $self->{_dbh}->selectrow_hashref($sql, undef, ($self->id));
    unless ($chd_ref){
	return undef;
    }
    my $chd = $chd_ref->{value};
    return $chd;
}









1;
