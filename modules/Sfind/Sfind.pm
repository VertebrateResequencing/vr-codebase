package Sfind::Sfind;
=head1 NAME

Sfind::Sfind - API for the denormalised SeqScape & NPG tracking database

=head1 SYNOPSIS
    my $sfind = Sfind::Sfind->new();

    #get arrayref of projects being tracked for traversing hierarchy
    my $projects = $sfind->projects();

    #also provides accessors for arbitrary objects in hierarchy
    my $sample = $sfind->get_sample_by_name('NA12878');

=head1 DESCRIPTION

Retrieves data from the sequencing warehouse database.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;
no warnings 'uninitialized';
use DBI;
use Sfind::Project;


=head2 new

  Arg [1]    : None
  Example    : my $sfind = Sfind::Sfind->new()
  Description: Returns Sfind object if can connect to database
  Returntype : Sfind::Sfind object

=cut

sub new {
    my ($class) = @_;

    my $self = {};
    bless ($self, $class);
    #my $dbh = DBI->connect("DBI:mysql:host=mcs2a:port=3313;database=warehouse_staging", "warehouse_ro",undef, {'RaiseError' => 1, 'PrintError'=>0});
    #my $dbh = DBI->connect("DBI:mysql:host=psdp:port=3305;database=warehouse_production", "warehouse_ro",undef, {'RaiseError' => 1, 'PrintError'=>0});
    # changed feb 3rd 2010
    # my $dbh = DBI->connect("DBI:mysql:host=psdp:port=3306;database=warehouse_production", "warehouse_ro",undef, {'RaiseError' => 1, 'PrintError'=>0});
    # changed Mar 31 2010
    my $dbh = DBI->connect("DBI:mysql:host=mcs7:port=3306;database=warehouse_production", "warehouse_ro",undef, {'RaiseError' => 1, 'PrintError'=>0});
#    my $dbh = DBI->connect("DBI:mysql:host=mcs6:port=3321;database=warehouse_staging", "warehouse_ro",undef, {'RaiseError' => 1, 'PrintError'=>0});
    $self->{_dbh} = $dbh;

    return $self;
}


=head2 get_project_by_id

  Arg [1]    : project id from sequencescape
  Example    : my $project = $sfind->get_project_by_id(140);
  Description: retrieve project object by sequencescape id
  Returntype : Sfind::Project object

=cut

sub get_project_by_id {
    my ($self, $id) = @_;
    my $obj = Sfind::Project->new($self->{_dbh},$id);
    return $obj;
}


=head2 get_project_by_name

  Arg [1]    : project name in sequencescape
  Example    : my $project = $sfind->get_project_by_name('1000Genomes-A1-CEU');
  Description: retrieve project object by sequencescape name
  Returntype : Sfind::Project object

=cut

sub get_project_by_name {
    my ($self, $name) = @_;
    my $sql = qq[select project_id from project_information where project_name=?];
    my $id_ref = $self->{_dbh}->selectrow_hashref($sql, undef, ($name));
    unless ($id_ref){
	warn "No project with name $name\n";
	return undef;
    }

    my $id = $id_ref->{project_id};
    return $self->get_project_by_id($id);
}


=head2 project_names

  Arg [1]    : None
  Example    : my $project_names = $sfind->project_names();
  Description: Returns a ref to an array of the project names that are being tracked
  Returntype : ref to array of project name strings

=cut

sub project_names {
    my ($self) = @_;

    unless ($self->{'project_names'}){
	my $sql = qq[select distinct project_name from project_information];
	my @projects;
	my $sth = $self->{_dbh}->prepare($sql);

	$sth->execute();
	foreach(@{$sth->fetchall_arrayref()}){
	    push @projects, $_->[0];
	}
	$self->{'project_names'} = \@projects;
    }

    return $self->{'project_names'};
}


=head2 get_sample_by_name

  Arg [1]    : sample name
  Arg [2]    : project name.  Only required if sample is in > 1 project
  Example    : my $sample = $sfind->get_sample_by_name('NA12878');
  Description: retrieve sample object by sample name.  If project is not supplied and sample is in more than one project, sample object creation will fail.
  Returntype : Sfind::Sample object

=cut

sub get_sample_by_name {
    my ($self, $name, $projname) = @_;
    my ($sample_id, $project_id) = $self->get_sample_project_id_by_name($name,$projname);
    my $obj;
    if ($sample_id && $project_id){
	$obj = $self->get_sample_by_id_project($sample_id, $project_id);
    }
    return $obj;
}


=head2 get_sample_project_id_by_name

  Arg [1]    : sample name
  Arg [2]    : project name.  Only required if sample is in > 1 project
  Example    : my ($sample_id,$project_id) = $sfind->get_sample_project_id_by_name('NA12878');
  Description: retrieve sample id by sample name.  If project is not supplied and sample is in more than one project, call will die.
  Returntype : sample_id, project_id

=cut

sub get_sample_project_id_by_name {
    my ($self, $name, $projname) = @_;
    my ($sample_id, $project_id);
    if($projname){
	my $sql = qq[select distinct sample_id,project_id from requests where sample_name=? and project_name=?];
	my $id_ref = $self->{_dbh}->selectall_arrayref($sql, undef, ($name, $projname));
	if (scalar @$id_ref == 0){
	    # warn "No sample with name $name in $projname\n";
	}
	elsif (scalar @$id_ref > 1){
	    die "More than one project id for name $projname\n";
	}
	else {
	    $sample_id = $id_ref->[0][0];
	    $project_id = $id_ref->[0][1];
	}
    }
    else { 
	# no project name supplied, so need to check if more than one	
	# project for this sample
	my $sql = qq[select distinct sample_id,project_id from requests where sample_name=?];
	my $id_ref = $self->{_dbh}->selectall_arrayref($sql, undef, ($name));
	if (scalar @$id_ref == 0){
	    # warn "No sample with name $name\n";
	}
	elsif (scalar @$id_ref > 1){
	    die "More than one project for sample $name\n";
	}
	else {
	    $sample_id = $id_ref->[0][0];
	    $project_id = $id_ref->[0][1];
	}
    }

    return ($sample_id, $project_id);
}


=head2 get_library_by_name

  Arg [1]    : library name
  Example    : my $library = $sfind->get_library_by_name('NA12878-1');
  Description: retrieve library object by library name.
  Returntype : Sfind::Library object

=cut

sub get_library_by_name {
    my ($self, $name) = @_;
    my $sql = qq[select distinct item_id from requests where item_name=?];
    my $id_ref = $self->{_dbh}->selectrow_hashref($sql, undef, ($name));
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
    my ($self, $id) = @_;
    my $obj = Sfind::Library->new($self->{_dbh},$id);
    return $obj;
}


=head2 get_sample_by_id_project

  Arg [1]    : sample name
  Arg [2]    : project name
  Example    : my $sample = $sfind->get_sample_by_id_project(1000,104);
  Description: retrieve sample object by sample name.  If project is not supplied and sample is in more than one project, sample object creation will fail.
  Returntype : Sfind::Sample object

=cut

sub get_sample_by_id_project {
    my ($self, $id, $projid) = @_;
    my $obj = Sfind::Sample->new($self->{_dbh},$id, $projid);
    return $obj;
}

=head2 get_projects_by_organism

  Arg [1]    : organism name
  Example    : my @proj = $sfind->get_projects_by_organism('Babesia bovis');
  Description: retrieve all project names for an organism given the organism name
  Returntype : array of project name strings

=cut

sub get_projects_by_organism {
    my ($self, $org) = @_;
    $org ="%$org%";   
    my $sql = qq[select distinct project_name from project_information join project_sample_reports on project_information.project_id = project_sample_reports.project_id join property_information on project_sample_reports.sample_id= property_information.obj_id where `key` like "sample_common_name" and  property_information.value like ?];
#   my @proj = $self->{_dbh}->selectall_arrayref($sql, undef, $org);
    my @proj;
    my $sth = $self->{_dbh}->prepare($sql);
    $sth->execute($org);
    foreach(@{$sth->fetchall_arrayref()}){
	    push @proj, $_->[0] if $_->[0];
	}
    return @proj;
}

=head2 get_taxonid_for_organism

  Arg [1]    : organism name
  Example    : my $taxon_id = $sfind->get_taxonid_for_organism('Homo sapiens');
  Description: Get the taxon ID for a given organism name
  Returntype : taxon_id

=cut

sub get_taxonid_for_organism {
    my ($self, $org) = @_;
    
    my $sql = qq[select distinct(sample.value) 
    		     from property_information sample 
    		     join property_information organism using (obj_id) 
    		     where organism.`key`="sample_common_name" 
    		     and organism.value=? 
    		     and sample.`key`="sample_taxon_id" 
    		     and sample.obj_type="Sample"
    		     limit 1;]; #TODO: Is this query faster if we use a subquery instead?

   	my $sth = $self->{_dbh}->prepare($sql); 								
   	$sth->execute($org);
   	            
	my ($taxon_id) = $sth->fetchrow_array();
	unless ($taxon_id){
		warn "No taxon_id found for $org\n";
		return undef;
    }
#    Some taxon ids in the warehouse have a .0 at the end which appears to be incorrect
#    We check and remove them here
    $taxon_id =~ s/(\d+)\.0/$1/g;
	
    return $taxon_id;
}


1;
