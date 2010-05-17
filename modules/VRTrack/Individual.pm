package VRTrack::Individual;

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

=head1 AUTHOR

jws@sanger.ac.uk (author)

=head1 METHODS

=cut

use strict;
use warnings;
use Carp qw(cluck confess);
use VRTrack::Species;
use VRTrack::Population;

use base qw(VRTrack::Named_obj
            VRTrack::Hierarchy_obj);


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
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    return $self;
}


=head2 fields_dispatch

  Arg [1]    : none
  Example    : my $fieldsref = $file->fields_dispatch();
  Description: Returns hashref dispatch table keyed on database field
               Used internally for new and update methods
  Returntype : hashref

=cut

sub fields_dispatch {
    my $self = shift;
    return {individual_id  => sub { $self->id(@_) },
	    name           => sub { $self->name(@_) },
	    hierarchy_name => sub { $self->hierarchy_name(@_) },
	    alias          => sub { $self->alias(@_) },
	    sex            => sub { $self->sex(@_) },
	    acc		   => sub { $self->acc(@_) },
	    species_id     => sub { $self->species_id(@_) },
	    population_id  => sub { $self->population_id(@_) }};
}


=head2 new_by_name

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : individual name
  Example    : my $ind = VRTrack::Individual->new($vrtrack, $name)
  Description: Class method. Returns Individual object by name.  If no such name is in the database, returns undef
  Returntype : VRTrack::Individual object

=cut


=head2 create

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : individual name
  Example    : my $ind = VRTrack::Individual->create($vrtrack, $name)
  Description: Class method. Creates new Individual object in the database.
  Returntype : VRTrack::Individual object

=cut

sub create {
    my ($self, $vrtrack, $name) = @_;
    my $hierarchy_name = $name;
    $hierarchy_name =~ s/\W+/_/g;
    return $self->SUPER::create($vrtrack, name => $name, hierarchy_name => $hierarchy_name);
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


=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $ind->id();
               $ind->id('104');
  Description: Get/Set for database ID of a individual
  Returntype : Internal ID integer

=cut


=head2 hierarchy_name

  Arg [1]    : directory name (optional)
  Example    : my $hname = $individual->hierarchy_name();
  Description: Get/set individual hierarchy name.  This is the directory name (without path) that the individual will be named in a file hierarchy.
  Returntype : string

=cut


=head2 name

  Arg [1]    : name (optional)
  Example    : my $name = $ind->name();
               $ind->name('Homo sapiens');
  Description: Get/Set for individual name
  Returntype : string

=cut


=head2 alias

  Arg [1]    : alias (optional)
  Example    : my $alias = $ind->alias();
               $ind->alias(1054);
  Description: Get/Set for individual name alias.  This is used for genotype checking, where the name we use for the individual may not be exactly the same as that specified in the genotyping files
  Returntype : string

=cut

sub alias {
    my $self = shift;
    return $self->_get_set('alias', 'string', @_);
}


=head2 sex

  Arg [1]    : sex (optional)
  Example    : my $sex = $ind->sex();
               $ind->sex('M');
  Description: Get/Set for individual sex
  Returntype : One of 'M', 'F', 'unknown'

=cut

sub sex {
    my $self = shift;
    if ($_[0]) {
	unless ($_[0] eq 'M' or $_[0] eq 'F') {
	    shift;
            unshift(@_, "unknown");
        }
    }
    return $self->_get_set('sex', 'string', @_);
}


=head2 acc

  Arg [1]    : acc (optional)
  Example    : my $acc = $study->acc();
               $samp->acc('ERS000090');
  Description: Get/Set for individual accession, i.e. SRA/ERA sample id
  Returntype : string

=cut

sub acc {
    my $self = shift;
    return $self->_get_set('acc', 'string', @_);
}


=head2 population

  Arg [1]    : population name (optional)
  Example    : my $population = $ind->population();
               $ind->population('CEU');
  Description: Get/Set for population individual belongs to.  Lazy-loads population object from $self->population_id.  If a population name is supplied, then population_id is set to the corresponding population in the database.  If no such population exists, returns undef.  Use add_population to add a population in this case.
  Returntype : VRTrack::Population object

=cut

sub population {
    my $self = shift;
    return $self->_get_set_child_object('get_population_by_name', 'VRTrack::Population', @_);
}


=head2 add_population

  Arg [1]    : population name
  Example    : my $pop = $ind->add_population('NOD_500_SLX_1');
  Description: create a new population, and if successful, return the object
  Returntype : VRTrack::Library object

=cut

sub add_population {
    my $self = shift;
    return $self->_create_child_object('get_population_by_name', 'VRTrack::Population', @_);
}


=head2 get_population_by_name

  Arg [1]    : population_name
  Example    : my $pop = $ind->get_population_by_name('CEU');
  Description: Retrieve a VRTrack::Population object by name
  Returntype : VRTrack::Population object

=cut

sub get_population_by_name {
    my ($self, $name) = @_;
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
    my $self = shift;
    return $self->_get_set('population_id', 'number', @_);
}


=head2 species_id

  Arg [1]    : species_id (optional)
  Example    : my $species_id = $samp->species_id();
               $samp->species_id(123);
  Description: Get/Set for sample internal species_id
  Returntype : integer

=cut

sub species_id {
    my $self = shift;
    return $self->_get_set('species_id', 'number', @_);
}


=head2 species

  Arg [1]    : species name (optional)
  Example    : my $species = $ind->species();
               $ind->species('Homo sapiens');
  Description: Get/Set for individual species.  Lazy-loads species object from $self->species_id.  If a species name is supplied, then species_id is set to the corresponding species in the database.  If no such species exists, returns undef.  Use add_species to add a species in this case.
  Returntype : VRTrack::Species object

=cut

sub species {
    my $self = shift;
    return $self->_get_set_child_object('get_species_by_name', 'VRTrack::Species', @_);
}


=head2 add_species

  Arg [1]    : species name
  Example    : my $spp = $ind->add_species('Homo sapiens');
  Description: create a new species, and if successful, return the object
  Returntype : VRTrack::Library object

=cut

sub add_species {
    my $self = shift;
    return $self->_create_child_object('get_species_by_name', 'VRTrack::Species', @_);
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


=head2 vrtrack

  Arg [1]    : vrtrack (optional)
  Example    : my $vrtrack = $obj->vrtrack();
               $obj->vrtrack($vrtrack);
  Description: Get/Set for vrtrack object.  NB you probably really shouldn't be setting vrtrack from outside this object unless you know what you're doing.
  Returntype : integer

=cut

1;
