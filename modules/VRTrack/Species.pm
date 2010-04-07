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

=head1 CONTACT

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


=head2 update

  Arg [1]    : None
  Example    : $species->update();
  Description: Update a species whose properties you have changed.  If properties haven't changed (i.e. dirty flag is unset) do nothing.  
               Unsets the dirty flag on success.
  Returntype : 1 if successful, otherwise undef.

=cut

1;
