package VRTrack::Assembly;

=head1 NAME

VRTrack::Assembly - Sequence Tracking Assembly object

=head1 SYNOPSIS
    my $assembly = VRTrack::Assembly->new($vrtrack, $assembly_id);

    my $id      = $assembly->id();
    my $name    = $assembly->name();

=head1 DESCRIPTION

An object describing a reference assembly sequence.
Used by VRTrack::Mapstats.

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

  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : assembly id
  Example    : my $assembly = VRTrack::Assembly->new($vrtrack, $id)
  Description: Returns Assembly object by assembly_id
  Returntype : VRTrack::Assembly object

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
    return {assembly_id => sub { $self->id(@_)},
            name        => sub { $self->name(@_)}};
}


=head2 new_by_name

  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : assembly name
  Example    : my $assembly = VRTrack::Assembly->new_by_name($vrtrack, $name)
  Description: Class method. Returns Assembly object by name and project_id.  If no such name is in the database, returns undef
  Returntype : VRTrack::Assembly object

=cut


=head2 create

  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : assembly name
  Example    : my $assembly = VRTrack::Assembly->create($vrtrack, $name)
  Description: Class method.  Creates new Assembly object in the database.
  Returntype : VRTrack::Assembly object

=cut

sub create {
    my ($self, $vrtrack, $value) = @_;
    return $self->SUPER::create($vrtrack, name => $value);
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
  Example    : my $id = $assembly->id();
               $assembly->id('104');
  Description: Get/Set for database ID of a assembly
  Returntype : Internal ID integer

=cut


=head2 name

  Arg [1]    : name (optional)
  Example    : my $name = $assembly->name();
               $assembly->name('NCBIm36');
  Description: Get/Set for assembly name
  Returntype : string

=cut


=head2 update

  Arg [1]    : None
  Example    : $assembly->update();
  Description: Update a assembly whose properties you have changed.  If properties haven't changed (i.e. dirty flag is unset) do nothing.  
               Unsets the dirty flag on success.
  Returntype : 1 if successful, otherwise undef.

=cut

1;
