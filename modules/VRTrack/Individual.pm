package VRTrack::Individual;
# author: jws
=head1 NAME

VRTrack::Individual - Sequence Tracking Individual object

=head1 SYNOPSIS
    my $ind = VRTrack::Individual->new($vrtrack, $individual_id);

    my $id      = $individual->id();
    my $name    = $individual->name();
    my $alias   = $individual->alias();
    my $sex     = $individual->sex();

=head1 DESCRIPTION

An object describing an individual, i.e. the entity that a sample is taken
from.  Individual are usually attached to a VRTrack::Sample by the
individual_id on the sample.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;
use Carp;
no warnings 'uninitialized';
use VRTrack::Species;
use VRTrack::Population;
use constant DBI_DUPLICATE => '1062';

###############################################################################
# Class methods
###############################################################################

=head2 new

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : individual id
  Example    : my $ind = VRTrack::Individual->new($vrtrack, $id)
  Description: Returns Individual object by individual_id
  Returntype : VRTrack::Individual object

=cut

sub new {
    my ($class,$vrtrack, $id) = @_;
    die "Need to call with a vrtrack handle and id" unless ($vrtrack && $id);
    if ( $vrtrack->isa('DBI::db') ) { croak "The interface has changed, expected vrtrack reference.\n"; }
    my $dbh = $vrtrack->{_dbh};
    my $self = {};
    bless ($self, $class);
    $self->{_dbh} = $dbh;
    $self->{vrtrack} = $vrtrack;

    my $sql = qq[select individual_id, name, hierarchy_name, alias, sex, acc, species_id, population_id from individual where individual_id = ?];
    my $sth = $self->{_dbh}->prepare($sql);

    if ($sth->execute($id)){
        my $data = $sth->fetchrow_hashref;
        unless ($data){
            return undef;
        }
        $self->id($data->{'individual_id'});
        $self->name($data->{'name'});
        $self->hierarchy_name($data->{'hierarchy_name'});
        $self->alias($data->{'alias'});
        $self->sex($data->{'sex'});
        $self->acc($data->{'acc'});
        $self->species_id($data->{'species_id'});
        $self->population_id($data->{'population_id'});
	$self->dirty(0); # unset the dirty flag
    }
    else{
        die(sprintf('Cannot retrieve individual: %s', $DBI::errstr));
    }

    return $self;
}


=head2 new_by_name

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : individual name
  Example    : my $ind = VRTrack::Individual->new($vrtrack, $name)
  Description: Class method. Returns Individual object by name.  If no such name is in the database, returns undef
  Returntype : VRTrack::Individual object

=cut

sub new_by_name {
    my ($class,$vrtrack, $name) = @_;
    die "Need to call with a vrtrack handle and name" unless ($vrtrack && $name);
    if ( $vrtrack->isa('DBI::db') ) { croak "The interface has changed, expected vrtrack reference.\n"; }
    my $dbh = $vrtrack->{_dbh};
    my $sql = qq[select individual_id from individual where name = ?];
    my $sth = $dbh->prepare($sql);

    my $id;
    if ($sth->execute($name)){
        my $data = $sth->fetchrow_hashref;
        unless ($data){
            return undef;
        }
        $id = $data->{'individual_id'};
    }
    else{
        die(sprintf('Cannot retrieve individual by name $name: %s', $DBI::errstr));
    }
    return $class->new($vrtrack, $id);
}


=head2 create

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : individual name
  Example    : my $ind = VRTrack::Individual->create($vrtrack, $name)
  Description: Class method. Creates new Individual object in the database.
  Returntype : VRTrack::Individual object

=cut

sub create {
    my ($class,$vrtrack, $name) = @_;
    die "Need to call with a vrtrack handle and name" unless ($vrtrack && $name);
    if ( $vrtrack->isa('DBI::db') ) { croak "The interface has changed, expected vrtrack reference.\n"; }
    my $dbh = $vrtrack->{_dbh};

    my $hierarchy_name = $name;
    $hierarchy_name =~ s/\W+/_/g;

    my $sql = qq[INSERT INTO individual (individual_id, name, hierarchy_name) 
                 VALUES (NULL,?,?)];

                
    my $sth = $dbh->prepare($sql);
    my $id;
    if ($sth->execute( $name, $hierarchy_name)) {
        $id = $dbh->{'mysql_insertid'};
    }
    else {
        die( sprintf('DB load insert failed: %s %s', $name, $DBI::errstr));
    }

    return $class->new($vrtrack, $id);

}


###############################################################################
# Object methods
###############################################################################


=head2 dirty

  Arg [1]    : boolean for dirty status
  Example    : $obj->dirty(1);
  Description: Get/Set for object properties having been altered.
  Returntype : boolean

=cut

sub dirty {
    my ($self,$dirty) = @_;
    if (defined $dirty){
	$self->{_dirty} = $dirty ? 1 : 0;
    }
    return $self->{_dirty};
}


=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $ind->id();
               $ind->id('104');
  Description: Get/Set for database ID of a individual
  Returntype : Internal ID integer

=cut

sub id {
    my ($self,$id) = @_;
    if (defined $id and $id != $self->{'id'}){
        $self->{'id'} = $id;
	$self->dirty(1);
    }
    return $self->{'id'};
}


=head2 hierarchy_name

  Arg [1]    : directory name (optional)
  Example    : my $hname = $individual->hierarchy_name();
  Description: Get/set individual hierarchy name.  This is the directory name (without path) that the individual will be named in a file hierarchy.
  Returntype : string

=cut

sub hierarchy_name {
    my ($self,$name) = @_;
    if (defined $name and $name ne $self->{'hierarchy_name'}){
        $self->{'hierarchy_name'} = $name;
	$self->dirty(1);
    }
    return $self->{'hierarchy_name'};
}


=head2 name

  Arg [1]    : name (optional)
  Example    : my $name = $ind->name();
               $ind->name('Homo sapiens');
  Description: Get/Set for individual name
  Returntype : string

=cut

sub name {
    my ($self,$name) = @_;
    if (defined $name and $name ne $self->{'name'}){
        $self->{'name'} = $name;
	$self->dirty(1);
    }
    return $self->{'name'};
}


=head2 alias

  Arg [1]    : alias (optional)
  Example    : my $alias = $ind->alias();
               $ind->alias(1054);
  Description: Get/Set for individual name alias.  This is used for genotype checking, where the name we use for the individual may not be exactly the same as that specified in the genotyping files
  Returntype : string

=cut

sub alias {
    my ($self,$alias) = @_;
    if (defined $alias and $alias ne $self->{'alias'}){
        $self->{'alias'} = $alias;
	$self->dirty(1);
    }
    return $self->{'alias'};
}


=head2 sex

  Arg [1]    : sex (optional)
  Example    : my $sex = $ind->sex();
               $ind->sex('M');
  Description: Get/Set for individual sex
  Returntype : One of 'M', 'F', 'unknown'

=cut

sub sex {
    my ($self,$sex) = @_;
    if (defined $sex and $sex ne $self->{'sex'}){
        unless ($sex eq 'M' or $sex eq 'F'){
            $sex = "unknown";
        }
        $self->{'sex'} = $sex;
	$self->dirty(1);
    }
    return $self->{'sex'};
}


=head2 acc

  Arg [1]    : acc (optional)
  Example    : my $acc = $study->acc();
               $samp->acc('ERS000090');
  Description: Get/Set for individual accession, i.e. SRA/ERA sample id
  Returntype : string

=cut

sub acc {
    my ($self,$acc) = @_;
    if (defined $acc and $acc ne $self->{'acc'}){
        $self->{'acc'} = $acc;
	$self->dirty(1);
    }
    return $self->{'acc'};
}


=head2 population

  Arg [1]    : population name (optional)
  Example    : my $population = $ind->population();
               $ind->population('CEU');
  Description: Get/Set for population individual belongs to.  Lazy-loads population object from $self->population_id.  If a population name is supplied, then population_id is set to the corresponding population in the database.  If no such population exists, returns undef.  Use add_population to add a population in this case.
  Returntype : VRTrack::Population object

=cut

sub population {
    my ($self,$population) = @_;
    if ($population){
        # get existing population by name
        my $obj = $self->get_population_by_name($population);
        if ($obj){
            # Have we actually changed?
            if ($self->population_id != $obj->id){
                # do in this order because setting id clears the object cache
                $self->population_id($obj->id);
                $self->dirty(1);
            }
            $self->{'population'} = $obj;
        }
        else {
            # warn "No such population in the database";
            return undef; # explicitly return nothing.
        }
    }
    elsif ($self->{'population'}){
        # already got a population object.  We'll return it at the end.
    }
    else {  # lazy-load population from database
        if ($self->population_id){
            my $obj = VRTrack::Population->new($self->{vrtrack},$self->population_id);
            $self->{'population'} = $obj;
        }
    }
    return $self->{'population'};
}


=head2 add_population

  Arg [1]    : population name
  Example    : my $pop = $ind->add_population('NOD_500_SLX_1');
  Description: create a new population, and if successful, return the object
  Returntype : VRTrack::Library object

=cut

sub add_population {
    my ($self, $name) = @_;

    my $obj = $self->get_population_by_name($name);
    if ($obj){
        warn "Population $name is already present in the database\n";
        return undef;
    }
    else {
        my $pop = VRTrack::Population->create($self->{vrtrack}, $name);
        # populate caches
        $self->{'population_id'} = $pop->id;
        $self->{'population'} = $pop;
    }
    return $self->{'population'};
}


=head2 get_population_by_name

  Arg [1]    : population_name
  Example    : my $pop = $ind->get_population_by_name('CEU');
  Description: Retrieve a VRTrack::Population object by name
  Returntype : VRTrack::Population object

=cut

sub get_population_by_name {
    my ($self,$name) = @_;
    return VRTrack::Population->new_by_name($self->{vrtrack}, $name);
}


=head2 population_id

  Arg [1]    : population_id (optional)
  Example    : my $population_id = $ind->population_id();
               $ind->population_id(123);
  Description: Get/Set for individual internal population_id.  Setting this clears any cached population object.
  Returntype : integer

=cut

sub population_id {
    my ($self,$population_id) = @_;
    if (defined $population_id and $population_id != $self->{'population_id'}){
        delete $self->{'population'};
        $self->{'population_id'} = $population_id;
	$self->dirty(1);
    }
    return $self->{'population_id'};
}


=head2 species_id

  Arg [1]    : species_id (optional)
  Example    : my $species_id = $samp->species_id();
               $samp->species_id(123);
  Description: Get/Set for sample internal species_id
  Returntype : integer

=cut

sub species_id {
    my ($self,$species_id) = @_;
    if (defined $species_id and $species_id != $self->{'species_id'}){
        delete $self->{'species'};
        $self->{'species_id'} = $species_id;
	$self->dirty(1);
    }
    return $self->{'species_id'};
}


=head2 species

  Arg [1]    : species name (optional)
  Example    : my $species = $ind->species();
               $ind->species('Homo sapiens');
  Description: Get/Set for individual species.  Lazy-loads species object from $self->species_id.  If a species name is supplied, then species_id is set to the corresponding species in the database.  If no such species exists, returns undef.  Use add_species to add a species in this case.
  Returntype : VRTrack::Species object

=cut

sub species {
    my ($self,$species) = @_;
    if ($species){
        # get existing species by name
        my $obj = $self->get_species_by_name($species);
        if ($obj){
            # Have we actually changed?
            if ($self->species_id != $obj->id){
                # do in this order because setting id clears the object cache
                $self->species_id($obj->id);
                $self->dirty(1);
            }
            $self->{'species'} = $obj;
        }
        else {
            # warn "No such species in the database";
            return undef; # explicitly return nothing.
        }
    }
    elsif ($self->{'species'}){
        # already got a species object.  We'll return it at the end.
    }
    else {  # lazy-load species from database
        if ($self->species_id){
            my $obj = VRTrack::Species->new($self->{vrtrack},$self->species_id);
            $self->{'species'} = $obj;
        }
    }
    return $self->{'species'};
}


=head2 add_species

  Arg [1]    : species name
  Example    : my $spp = $ind->add_species('Homo sapiens');
  Description: create a new species, and if successful, return the object
  Returntype : VRTrack::Library object

=cut

sub add_species {
    my ($self, $name) = @_;

    my $obj = $self->get_species_by_name($name);
    if ($obj){
        warn "Species $name is already present in the database\n";
        return undef;
    }
    else {
        my $pop = VRTrack::Species->create($self->{vrtrack}, $name);
        # populate caches
        $self->{'species_id'} = $pop->id;
        $self->{'species'} = $pop;
    }
    return $self->{'species'};
}


=head2 get_species_by_name

  Arg [1]    : species_name
  Example    : my $pop = $ind->get_species_by_name('Homo sapiens');
  Description: Retrieve a VRTrack::Species object by name
  Returntype : VRTrack::Species object

=cut

sub get_species_by_name {
    my ($self,$name) = @_;
    return VRTrack::Species->new_by_name($self->{vrtrack}, $name);
}


=head2 update

  Arg [1]    : None
  Example    : $individual->update();
  Description: Update a individual whose properties you have changed.  If properties haven't changed (i.e. dirty flag is unset) do nothing.  
               Unsets the dirty flag on success.
  Returntype : 1 if successful, otherwise undef.

=cut

sub update {
    my ($self) = @_;
    my $success = undef;
    if ($self->dirty){
	my $dbh = $self->{_dbh};
	my $save_re = $dbh->{RaiseError};
	my $save_pe = $dbh->{PrintError};
	$dbh->{RaiseError} = 1; # raise exception if an error occurs
	$dbh->{PrintError} = 0; # don't print an error message

	eval {
	    my $updsql = qq[UPDATE individual SET name=?, hierarchy_name=?,alias=?, sex=?, species_id=?, population_id=? WHERE individual_id = ? ];
	    
	    $dbh->do ($updsql, undef, $self->name, $self->hierarchy_name, $self->alias, $self->sex,, $self->species_id, $self->population_id, $self->id);
	};

	if (!$@) {
	    $success = 1;
	}

	# restore attributes to original state
	$dbh->{PrintError} = $save_pe;
	$dbh->{RaiseError} = $save_re;

    }
    if ($success){
        $self->dirty(0);
    }

    return $success;
}

1;
