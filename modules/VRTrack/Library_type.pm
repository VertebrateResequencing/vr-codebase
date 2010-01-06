package VRTrack::Library_type;

=head1 NAME

VRTrack::Library_type - Sequence Tracking Library_type object

=head1 SYNOPSIS
    my $library_type = VRTrack::Library_type->new($vrtrack, $library_type_id);

    my $id      = $library_type->id();
    my $name    = $library_type->name();

=head1 DESCRIPTION

An object describing a library_type, such as DSS or NOPCR.
Library_types are usually attached to a VRTrack::Library by the
library_type_id on the library.

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
  Arg [2]    : library_type id
  Example    : my $library_type = VRTrack::Library_type->new($vrtrack, $id)
  Description: Returns Library_type object by library_type_id
  Returntype : VRTrack::Library_type object

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
    return {library_type_id  => sub { $self->id(@_) },
	    name             => sub { $self->name(@_) }};
}


=head2 new_by_name

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : library_type name
  Example    : my $library_type = VRTrack::Library_type->new_by_name($vrtrack, $name)
  Description: Class method. Returns Library_type object by name and project_id.  If no such name is in the database, returns undef
  Returntype : VRTrack::Library_type object

=cut


=head2 create

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : library_type name
  Example    : my $library_type = VRTrack::Library_type->create($vrtrack, $name)
  Description: Class method.  Creates new Library_type object in the database.
  Returntype : VRTrack::Library_type object

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
  Example    : my $id = $library_type->id();
               $library_type->id('104');
  Description: Get/Set for database ID of a library_type
  Returntype : Internal ID integer

=cut


=head2 name

  Arg [1]    : name (optional)
  Example    : my $name = $library_type->name();
               $library_type->name('CEU');
  Description: Get/Set for library_type name
  Returntype : string

=cut


=head2 update

  Arg [1]    : None
  Example    : $library_type->update();
  Description: Update a library_type whose properties you have changed.  If properties haven't changed (i.e. dirty flag is unset) do nothing.  
               Unsets the dirty flag on success.
  Returntype : 1 if successful, otherwise undef.

=cut

1;
