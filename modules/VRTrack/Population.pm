package VRTrack::Population;

=head1 NAME

VRTrack::Population - Sequence Tracking Population object

=head1 SYNOPSIS
    my $pop = VRTrack::Population->new($vrtrack, $population_id);

    my $id      = $population->id();
    my $name    = $population->name();

=head1 DESCRIPTION

An object describing a population.
Populations are usually attached to a VRTrack::Sample by the
population_id on the individual.

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
  Arg [2]    : population id
  Example    : my $pop = VRTrack::Population->new($vrtrack, $id)
  Description: Class method. Returns Population object by population_id
  Returntype : VRTrack::Population object

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
    return {population_id  => sub { $self->id(@_) },
	    name           => sub { $self->name(@_) }};
}


=head2 new_by_name

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : population name
  Example    : my $pop = VRTrack::Population->new($vrtrack, $name)
  Description: Class method. Returns Population object by name.  If no such name is in the database, returns undef
  Returntype : VRTrack::Population object

=cut


=head2 create

  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : population name
  Example    : my $pop = VRTrack::Population->create($vrtrack, $name)
  Description: Class method. Creates new Population object in the database.
  Returntype : VRTrack::Population object

=cut

sub create {
    my ($self, $vrtrack, $name) = @_;
    return $self->SUPER::create($vrtrack, name => $name);
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
  Example    : my $id = $pop->id();
               $pop->id('104');
  Description: Get/Set for database ID of a population
  Returntype : Internal ID integer

=cut


=head2 name

  Arg [1]    : name (optional)
  Example    : my $name = $pop->name();
               $pop->name('CEU');
  Description: Get/Set for population name
  Returntype : string

=cut


=head2 update

  Arg [1]    : None
  Example    : $population->update();
  Description: Update a population whose properties you have changed.  If properties haven't changed (i.e. dirty flag is unset) do nothing.  
               Unsets the dirty flag on success.
  Returntype : 1 if successful, otherwise undef.

=cut

1;
