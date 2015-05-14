package VRTrack::Species;

=head1 NAME

VRTrack::Species - Sequence Tracking Species object

=head1 SYNOPSIS
    my $spp = VRTrack::Species->new($vrtrack, $species_id);

    my $id      = $species->id();
    my $name    = $species->name();
    my $taxon   = $species->taxon_id();

=head1 DESCRIPTION

An object describing a species.  Species are usually attached to a
VRTrack::Sample by the species_id on the individual.

=head1 AUTHOR

jws@sanger.ac.uk (author)

=head1 METHODS

=cut

use strict;
use warnings;
use Carp qw(cluck confess);

use base qw(VRTrack::Named_obj);


###############################################################################
# Class methods
###############################################################################

=head2 new

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : species id
  Example    : my $spp = VRTrack::Species->new($vrtrack, $id)
  Description: Returns Species object by species_id
  Returntype : VRTrack::Species object

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
    return {species_id => sub { $self->id(@_) },
	    name       => sub { $self->name(@_) },
	    taxon_id   => sub { $self->taxon_id(@_) }};
}


=head2 new_by_name

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : species name
  Example    : my $spp = VRTrack::Species->new($vrtrack, $name)
  Description: Class method. Returns Species object by name.  If no such name is in the database, returns undef
  Returntype : VRTrack::Species object

=cut

=head2 new_by_taxon_id

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : taxon id
  Example    : my $spp = VRTrack::Species->new_by_taxon_id($vrtrack, $taxon)
  Description: Class method. Returns Species object by taxon id.  If no such taxon id is in the database, returns undef
  Returntype : VRTrack::Species object

=cut

sub new_by_taxon_id {
    my ($class, $vrtrack, $value) = @_;
    return $class->new_by_field_value($vrtrack, 'taxon_id', $value);
}

=head2 create

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : species name
  Example    : my $spp = VRTrack::Species->create($vrtrack, $name)
  Description: Class method. Creates new Species object in the database.
  Returntype : VRTrack::Species object

=cut

sub create {
    my ($self, $vrtrack, $name, $taxon_id) = @_;
    if(defined $taxon_id){
    	return $self->SUPER::create($vrtrack, name => $name, taxon_id => $taxon_id);
    }else{
    	return $self->SUPER::create($vrtrack, name => $name);
    }	
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
  Example    : my $id = $spp->id();
               $spp->id('104');
  Description: Get/Set for database ID of a species
  Returntype : Internal ID integer

=cut


=head2 name

  Arg [1]    : name (optional)
  Example    : my $name = $spp->name();
               $spp->name('Homo sapiens');
  Description: Get/Set for species name
  Returntype : string

=cut


=head2 taxon_id

  Arg [1]    : taxon_id (optional)
  Example    : my $taxon_id = $spp->taxon_id();
               $spp->taxon_id(1054);
  Description: Get/Set for species taxon id
  Returntype : integer

=cut

sub taxon_id {
    my $self = shift;
    return $self->_get_set('taxon_id', 'number', @_);
}


=head 2 genus

  Example    : my $genus = $species->genus();
  Description: Get the genus name of the species (for now, this is the first word)
  Returntype : string

=cut

sub genus {
    my $self = shift;
    my @whole_name = split(/ /, $self->name());
    return $whole_name[0];
}

=head 2 species_subspecies

  Example    : my $genus = $species->species_subspecies();
  Description: Gets everything after the genus name (for now, these are all the words after the first space)
  Returntype : string

=cut

sub species_subspecies {
    my $self = shift;
    $self->name() =~ m/(\S+)\s+(.*)/;  
    return $2;
}


=head2 update

  Arg [1]    : None
  Example    : $species->update();
  Description: Update a species whose properties you have changed.  If properties haven't changed (i.e. dirty flag is unset) do nothing.  
               Unsets the dirty flag on success.
  Returntype : 1 if successful, otherwise undef.

=cut



1;
