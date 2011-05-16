package Sfind::Sfind;

=head1 NAME

Sfind::Sfind - API for the v2 denormalised SeqScape & NPG tracking database

=head1 SYNOPSIS
    my $sfind = Sfind::Sfind->new();

    #get arrayref of studies being tracked for traversing hierarchy
    my $studies = $sfind->studies();

    #also provides accessors for arbitrary objects in hierarchy
    my $sample = $sfind->get_sample_by_name('NA12878');

=head1 DESCRIPTION

Retrieves data from the sequencing warehouse_two database.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut

use Moose;
use namespace::autoclean;
use DBI;
use Sfind::Study;


has '_dbh' => (
    is          => 'ro',
    isa         => 'DBI::db',
    required    => 1,
    lazy        => 1,
    builder     => '_db_connect',
    init_arg    => undef,
);


sub _db_connect { 
    my ($self) = @_; 
    my $dbh = DBI->connect("DBI:mysql:host=mcs7:port=3306;database=warehouse_two_production", "warehouse_ro",undef, {'RaiseError' => 1, 'PrintError'=>0});
    return $dbh;
} 



=head2 get_study_by_id

  Arg [1]    : study id from sequencescape
  Example    : my $study = $sfind->get_study_by_id(140);
  Description: retrieve study object by sequencescape id
  Returntype : Sfind::Study object

=cut

sub get_study_by_id {
    my ($self, $id) = @_;
    my $obj = Sfind::Study->new({dbh   => $self->_dbh,
                                 id     => $id});
    return $obj;
}


=head2 get_study_by_name

  Arg [1]    : study name in sequencescape
  Example    : my $study = $sfind->get_study_by_name('1000Genomes-A1-CEU');
  Description: retrieve study object by sequencescape name
  Returntype : Sfind::Study object

=cut


sub get_study_by_name {
    my ($self, $name) = @_;
    my $sql = qq[select internal_id from current_studies where name=? ];
    my $id_ref = $self->_dbh->selectrow_hashref($sql, undef, ($name));
    unless ($id_ref){
	warn "No study with name $name\n";
	return undef;
    }

    my $id = $id_ref->{internal_id};
    return $self->get_study_by_id($id);
}



=head2 study_names

  Arg [1]    : None
  Example    : my $study_names = $sfind->study_names();
  Description: Returns a ref to an array of the study names that are being tracked
  Returntype : ref to array of study name strings

=cut

sub study_names {
    my ($self) = @_;

    unless ($self->{'study_names'}){
	my $sql = qq[select distinct name from current_studies ];
	my @studies;
	my $sth = $self->_dbh->prepare($sql);

	$sth->execute();
	foreach(@{$sth->fetchall_arrayref()}){
	    push @studies, $_->[0];
	}
	$self->{'study_names'} = \@studies;
    }
    return $self->{'study_names'};
}


=head2 get_sample_by_name

  Arg [1]    : sample name
  Arg [2]    : study name.  Only required if sample is in > 1 study
  Example    : my $sample = $sfind->get_sample_by_name('NA12878');
  Description: retrieve sample object by sample name.  If study is not supplied and sample is in more than one study, sample object creation will fail.
  Returntype : Sfind::Sample object

=cut

sub get_sample_by_name {
    die "NOT YET IMPLEMENTED";
    my ($self, $name, $studyname) = @_;
    my ($sample_id, $study_id) = $self->get_sample_study_id_by_name($name,$studyname);
    my $obj;
    if ($sample_id && $study_id){
	$obj = $self->get_sample_by_id_study($sample_id, $study_id);
    }
    return $obj;
}


=head2 get_sample_study_id_by_name

  Arg [1]    : sample name
  Arg [2]    : study name.  Only required if sample is in > 1 study
  Example    : my ($sample_id,$study_id) = $sfind->get_sample_study_id_by_name('NA12878');
  Description: retrieve sample id by sample name.  If study is not supplied and sample is in more than one study, call will die.
  Returntype : sample_id, study_id

=cut

sub get_sample_study_id_by_name {
    die "NOT YET IMPLEMENTED";
    my ($self, $name, $studyname) = @_;
    my ($sample_id, $study_id);
    if($studyname){
	my $sql = qq[select distinct requests.sample_id, requests.study_id 
                     from requests 
                     where sample_name=? 
                     and study_id=(select distinct study_id from study_information where study_name=?)];
	my $id_ref = $self->_dbh->selectall_arrayref($sql, undef, ($name, $studyname));
	if (scalar @$id_ref == 0){
	    #warn "No sample with name $name in $studyname\n";
	}
	elsif (scalar @$id_ref > 1){
	    die "More than one study id for name $studyname\n";
	}
	else {
	    $sample_id = $id_ref->[0][0];
	    $study_id = $id_ref->[0][1];
	}
    }
    else { 
	# no study name supplied, so need to check if more than one	
	# study for this sample
	my $sql = qq[select distinct sample_id, study_id from requests where sample_name=?];
	my $id_ref = $self->_dbh->selectall_arrayref($sql, undef, ($name));
	if (scalar @$id_ref == 0){
	    # warn "No sample with name $name\n";
	}
	elsif (scalar @$id_ref > 1){
	    die "More than one study for sample $name\n";
	}
	else {
	    $sample_id = $id_ref->[0][0];
	    $study_id = $id_ref->[0][1];
	}
    }

    return ($sample_id, $study_id);
}


=head2 get_library_by_name

  Arg [1]    : library name
  Example    : my $library = $sfind->get_library_by_name('NA12878-1');
  Description: retrieve library object by library name.
  Returntype : Sfind::Library object

=cut

sub get_library_by_name {
    die "NOT YET IMPLEMENTED";
    my ($self, $name) = @_;
    my $sql = qq[select distinct item_id from requests where item_name=?];
    my $id_ref = $self->_dbh->selectrow_hashref($sql, undef, ($name));
    unless ($id_ref){
	warn "No library with name $name\n";
	return undef;
    }
    my $id = $id_ref->{item_id};
    return $self->get_library_by_id($id);
}



=head2 get_library_by_id

  Arg [1]    : library id from sequencescape
  Example    : my $library = $sfind->get_library_by_id(140);
  Description: retrieve library object by sequencescape id
  Returntype : Sfind::Library object

=cut

sub get_library_by_id {
    die "NOT YET IMPLEMENTED";
    my ($self, $id) = @_;
    my $obj = Sfind::Library->new($self->_dbh,$id);
    return $obj;
}


=head2 get_sample_by_id_study

  Arg [1]    : sample name
  Arg [2]    : study name
  Example    : my $sample = $sfind->get_sample_by_id_study(1000,104);
  Description: retrieve sample object by sample id.  If study is not supplied and sample is in more than one study, sample object creation will fail.
  Returntype : Sfind::Sample object

=cut

sub get_sample_by_id_study {
    die "NOT YET IMPLEMENTED";
    my ($self, $id, $studyid) = @_;
    my $obj = Sfind::Sample->new($self->_dbh,$id, $studyid);
    return $obj;
}

=head2 get_studies_by_organism

  Arg [1]    : organism name
  Example    : my @studies = $sfind->get_studies_by_organism('Babesia bovis');
  Description: retrieve all study names for an organism given the organism name
  Returntype : array of study name strings

=cut

sub get_studies_by_organism {
    my ($self, $org) = @_;
    $org ="%$org%";   
    my $sql = qq[select distinct current_studies.name from current_studies join current_study_samples on current_studies.internal_id = current_study_samples.study_internal_id join current_samples on current_study_samples.sample_internal_id= current_samples.internal_id where current_samples.common_name like ? ];

    my @studies;
    my $sth = $self->_dbh->prepare($sql);
    $sth->execute($org);
    foreach(@{$sth->fetchall_arrayref()}){
	push @studies, $_->[0] if $_->[0];
    }
    return @studies;
}

=head2 get_taxonid_for_organism

  Arg [1]    : organism name
  Example    : my $taxon_id = $sfind->get_taxonid_for_organism('Homo sapiens');
  Description: Get the taxon ID for a given organism name
  Returntype : taxon_id

=cut

sub get_taxonid_for_organism {
    my ($self, $org) = @_;
    
    my $sql = qq[select distinct taxon_id from current_samples 
            where common_name=? and taxon_id is not null limit 1;];

   	my $sth = $self->_dbh->prepare($sql);
   	$sth->execute($org);
   	            
	my ($taxon_id) = $sth->fetchrow_array();
	unless ($taxon_id){
		warn "No taxon_id found for $org\n";
		return undef;
    }

    # Some taxon ids in the warehouse have a .0 at the end which seems to be
    # incorrect We check and remove them here

    $taxon_id =~ s/(\d+)\.0/$1/g;
	
    return $taxon_id;
}

__PACKAGE__->meta->make_immutable;

1;
