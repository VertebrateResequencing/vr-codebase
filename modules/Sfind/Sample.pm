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

=head2 new

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : sample id
  Arg [3]    : project id
  Example    : my $samp = Sfind::Sfind->new($dbh, $id, $pid)
  Description: Returns Sample object by sample_id and project_id
  Returntype : Sfind::Sample object

=cut

sub new {
    my ($class,$dbh, $id, $proj_id) = @_;
    die "Need to call with a db handle, id and project id" unless ($dbh && $id && $proj_id);
    my $self = {};
    bless ($self, $class);
    $self->{_dbh} = $dbh;

    my $sql = qq[select distinct(sample_name) from requests where sample_id=? and project_id = ?];
    my $id_ref = $self->{_dbh}->selectrow_hashref($sql, undef, ($id, $proj_id));
    if ($id_ref){
	my $name = $id_ref->{sample_name};
        $name =~ s/\s+$//;  # trim trailing whitespace
	#warn "Sample name : $name\n";
	$self->id($id);
	$self->name($name);
	$self->project_id($proj_id);
    }
    else {
	return undef;
    }
    return $self;
}



=head2 libraries

  Arg [1]    : None
  Example    : my $libraries = $sample->libraries();
  Description: Returns a ref to an array of the sample objects that are associated with this sample.
  Returntype : ref to array of Sfind::Sample objects

=cut

sub libraries {
    my ($self) = @_;

    unless ($self->{'libraries'}){
	my @libraries;
    	foreach my $id (@{$self->library_ids()}){
	    my $obj = Sfind::Library->new($self->{_dbh},$id);
	    push @libraries, $obj; 
	}
	@libraries = sort {$a <=> $b} @libraries;
	$self->{'libraries'} = \@libraries;
    }

    return $self->{'libraries'};
}


=head2 library_ids

  Arg [1]    : None
  Example    : my $library_ids = $sample->library_ids();
  Description: Returns a ref to an array of the library IDs that are associated with this sample
  Returntype : ref to array of integer library IDs

=cut

sub library_ids {
    my ($self) = @_;

    unless ($self->{'library_ids'}){
	my $sql = qq[select distinct(item_id) from requests where sample_id=? and project_id = ?];
	my @libraries;
	my $sth = $self->{_dbh}->prepare($sql);

	$sth->execute($self->id, $self->project_id);
	foreach(@{$sth->fetchall_arrayref()}){
	    push @libraries, $_->[0] if $_->[0];
	}
	$self->{'library_ids'} = \@libraries;
    }
    
    return $self->{'library_ids'};
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


=head2 project_id

  Arg [1]    : project_id (optional)
  Example    : my $project_id = $samp->project_id();
	       $samp->project_id('104');
  Description: Get/Set for Project ID of a sample
  Returntype : SequenceScape ID (usu. integer)

=cut

sub project_id {
    my ($self,$project_id) = @_;
    if ($project_id){
	$self->{'project_id'} = $project_id;
    }
    return $self->{'project_id'};
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


=head2 get_library_by_id

  Arg [1]    : library id from sequencescape
  Example    : my $library = $sam->get_library_by_id(1930);
  Description: retrieve library object by sequencescape id
  Returntype : Sfind::Library object

=cut

sub get_library_by_id {
    my ($self, $id) = @_;
    my $obj = Sfind::Library->new($self->{_dbh},$id);
    return $obj;
}


=head2 get_library_by_name

  Arg [1]    : library name from sequencescape
  Example    : my $library = $sam->get_library_by_name('NA18909-YRI-1');
  Description: retrieve library object by sequencescape name
  Returntype : Sfind::Library object

=cut

sub get_library_by_name {
    my ($self, $name) = @_;
    my $sql = qq[select distinct(item_id) from requests where item_name=? and sample_id = ? and project_id=?];
    my $id_ref = $self->{_dbh}->selectrow_hashref($sql, undef, ($name, $self->id, $self->project_id));
    unless ($id_ref){
	warn "No library with name $name\n";
	return undef;
    }

    my $id = $id_ref->{item_id};
    return $self->get_library_by_id($id);
}



=head2 get_organism_name

  Arg [1]    : None
  Example    : my $organism = $sample->get_organism_name();
  Description: retrieve organism name from given sample ID
  Returntype : string

=cut

sub get_organism_name {
    my ($self) = @_;
    my $sql = qq[select value from property_information where `key` like "organism" and property_information.obj_id=?];
    my $org = $self->{_dbh}->selectrow_hashref($sql, undef, ($self->id));
    unless ($org){
	warn "No organism ", $self->id,"\n";
	return undef;
    }
    my $orgname = $org->{value};
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
1;
