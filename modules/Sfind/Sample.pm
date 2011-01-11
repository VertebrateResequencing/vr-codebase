package Sfind::Sample;
=head1 NAME

Sfind::Sample - Sequence Tracking Sample object

=head1 SYNOPSIS
    my $samp = Sfind::Sample->new($dbh, $sample_id);

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

use strict;
use warnings;
no warnings 'uninitialized';
use Sfind::Library;
use Sfind::Library_Request;

=head2 new

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : sample id
  Arg [3]    : study id
  Example    : my $samp = Sfind::Sfind->new($dbh, $id, $sid)
  Description: Returns Sample object by sample_id and study_id
  Returntype : Sfind::Sample object

=cut

sub new {
    my ($class,$dbh, $id, $study_id) = @_;
    die "Need to call with a db handle, id and study id" unless ($dbh && $id && $study_id);
    my $self = {};
    bless ($self, $class);
    $self->{_dbh} = $dbh;

    # jws 2011-01-05
    # Changed
    # from: get sample_name from requests_new table on this study
    # to: get sample_name from samples table

    my $sql = qq[select s.name as sample_name from study_sample_reports ssr, samples s where ssr.sample_id = ? and ssr.study_id = ? and ssr.sample_id = s.sample_id];

    my $id_ref = $self->{_dbh}->selectrow_hashref($sql, undef, ($id, $study_id));
    if ($id_ref){
	my $name = $id_ref->{sample_name};
        $name =~ s/\s+$//;  # trim trailing whitespace
	#warn "Sample name : $name\n";
	$self->id($id);
	$self->name($name);
	$self->study_id($study_id);
    }
    else {
	return undef;
    }
    return $self;
}


=head2 library_requests

  Arg [1]    : library name from sequencescape
  Example    : my $librequests = $sample->library_requests();
  Description: Returns a ref to an array of the library_request objects that are associated with this sample.
  Returntype : ref to array of Sfind::Library_Request objects

=cut

sub library_requests {
   my ($self) = @_;

    unless ($self->{'library_requests'}){
	my @library_requests;
    	foreach my $id (@{$self->library_request_ids()}){
	    my $obj = Sfind::Library_Request->new($self->{_dbh},$id);
	    push @library_requests, $obj; 
	}
	@library_requests = sort {$a <=> $b} @library_requests;
	$self->{'library_requests'} = \@library_requests;
    }

    return $self->{'library_requests'};
}


=head2 library_request_ids

  Arg [1]    : None
  Example    : my $library_request_ids = $sample->library_request_ids();
  Description: Returns a ref to an array of the library request IDs that are associated with this sample
  Returntype : ref to array of integer library request IDs

=cut

sub library_request_ids {
    my ($self) = @_;
    unless ($self->{'library_request_ids'}){
	my $sql = qq[select request_id as lib_request_id from requests_new where type like '%library creation%' and sample_id=? and study_id=?];
	my @library_requests;
	my $sth = $self->{_dbh}->prepare($sql);
	$sth->execute($self->id, $self->study_id);
	foreach(@{$sth->fetchall_arrayref()}){
	    push @library_requests, $_->[0] if $_->[0];
	}
	$self->{'library_request_ids'} = \@library_requests;
    }
    
    return $self->{'library_request_ids'};
}

=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $samp->id();
	       $samp->id('104');
  Description: Get/Set for ID of a sample
  Returntype : SequenceScape ID (usu. integer)

=cut

sub id {
    my ($self,$id) = @_;
    if ($id){
	$self->{'id'} = $id;
    }
    return $self->{'id'};
}


=head2 study_id

  Arg [1]    : study_id (optional)
  Example    : my $study_id = $samp->study_id();
	       $samp->study_id('104');
  Description: Get/Set for study ID of a sample
  Returntype : SequenceScape ID (usu. integer)

=cut

sub study_id {
    my ($self,$study_id) = @_;
    if ($study_id){
	$self->{'study_id'} = $study_id;
    }
    return $self->{'study_id'};
}


=head2 name

  Arg [1]    : name (optional)
  Example    : my $name = $samp->name();
	       $samp->name('104');
  Description: Get/Set for sample name
  Returntype : string

=cut

sub name {
    my ($self,$name) = @_;
    if ($name){
	$self->{'name'} = $name;
    }
    return $self->{'name'};
}

=head2 get_organism_name

  Arg [1]    : None
  Example    : my $organism = $sample->get_organism_name();
  Description: retrieve organism name from sample_common_name key for given sample ID
  Returntype : string

=cut

sub get_organism_name {
    my ($self) = @_;
    my $sql = qq[select value from property_information where `key` like "sample_common_name" and property_information.obj_id=?];
    my $org = $self->{_dbh}->selectrow_hashref($sql, undef, ($self->id));
    unless ($org){
	warn "No organism ", $self->id,"\n";
	return undef;
    }
    my $orgname = $org->{value};
    $orgname =~ s/^\s+//; #remove leading spaces
    $orgname =~ s/\s+$//; #remove trailing spaces
    return $orgname;
}





=head2 get_strain_name

  Arg [1]    : None
  Example    : my $strain = $sample->get_strain_name();
  Description: retrieve strain information from given sample ID
  Returntype : string

=cut

sub get_strain_name {
    my ($self) = @_;
    my $sql = qq[select value from property_information where `key` like "sample_strain_att" and property_information.obj_id=?];
    my $strain = $self->{_dbh}->selectrow_hashref($sql, undef, ($self->id));
    unless ($strain){
	warn "No strain ", $self->id,"\n";
	return undef;
    }
    my $strain_name = $strain->{value};
    return $strain_name;
}

=head2 get_accession

  Arg [1]    : None
  Example    : my $acc = $sample->get_accession();
  Description: retrieve EBI accession number from given sample ID
  Returntype : string

=cut

sub get_accession {
    my ($self) = @_;
    my $sql = qq[select value from property_information where `key` like "sample_ebi_accession_number" and property_information.obj_id=?];
    my $acc = $self->{_dbh}->selectrow_hashref($sql, undef, ($self->id));
    unless ($acc){
	warn "No accession ", $self->id,"\n";
	return undef;
    }
    my $acc_name = $acc->{value};
    return $acc_name;
}

=head2 get_public_name

  Arg [1]    : None
  Example    : my $pubname = $sample->get_public_name();
  Description: retrieve sample public name from given sample ID
  Returntype : string

=cut

sub get_public_name {
    my ($self) = @_;
    my $sql = qq[select value from property_information where `key` like "sample_public_name" and property_information.obj_id=?];
    my $pub = $self->{_dbh}->selectrow_hashref($sql, undef, ($self->id));
    unless ($pub){
	warn "No public name ", $self->id,"\n";
	return undef;
    }
    my $pub_name = $pub->{value};
    return $pub_name;
}
1;
