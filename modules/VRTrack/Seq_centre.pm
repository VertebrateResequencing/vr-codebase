package VRTrack::Seq_centre;

=head1 NAME

VRTrack::Seq_centre - Sequence Tracking Seq_centre object

=head1 SYNOPSIS
    my $seq_centre = VRTrack::Seq_centre->new($vrtrack, $seq_centre_id);

    my $id      = $seq_centre->id();
    my $name    = $seq_centre->name();

=head1 DESCRIPTION

An object describing a sequencing_centre, such as SC or BROAD.
Seq_centres are usually attached to a VRTrack::Library by the
seq_centre_id on the library.

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
  Arg [2]    : seq_centre id
  Example    : my $seq_centre = VRTrack::Seq_centre->new($vrtrack, $id)
  Description: Returns Seq_centre object by seq_centre_id
  Returntype : VRTrack::Seq_centre object

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
    return {seq_centre_id  => sub { $self->id(@_) },
	    name           => sub { $self->name(@_) }};
}


=head2 new_by_name

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : seq_centre name
  Example    : my $seq_centre = VRTrack::Seq_centre->new_by_name($vrtrack, $name)
  Description: Class method. Returns Seq_centre object by seq_centre name.  If no such name is in the database, returns undef
  Returntype : VRTrack::Seq_centre object

=cut


=head2 create

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : seq_centre name
  Example    : my $seq_centre = VRTrack::Seq_centre->create($vrtrack, $name)
  Description: Class method.  Creates new Seq_centre object in the database.
  Returntype : VRTrack::Seq_centre object

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
  Example    : my $id = $seq_centre->id();
               $seq_centre->id('104');
  Description: Get/Set for database ID of a seq_centre
  Returntype : Internal ID integer

=cut


=head2 name

  Arg [1]    : name (optional)
  Example    : my $name = $seq_centre->name();
               $seq_centre->name('SC');
  Description: Get/Set for seq_centre name
  Returntype : string

=cut


=head2 update

  Arg [1]    : None
  Example    : $seq_centre->update();
  Description: Update a seq_centre whose properties you have changed.  If properties haven't changed (i.e. dirty flag is unset) do nothing.  
               Unsets the dirty flag on success.
  Returntype : 1 if successful, otherwise undef.

=cut

1;
