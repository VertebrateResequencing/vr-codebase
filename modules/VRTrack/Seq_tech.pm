package VRTrack::Seq_tech;

=head1 NAME

VRTrack::Seq_tech - Sequence Tracking Seq_tech object

=head1 SYNOPSIS
    my $seq_tech = VRTrack::Seq_tech->new($vrtrack, $seq_tech_id);

    my $id      = $seq_tech->id();
    my $name    = $seq_tech->name();

=head1 DESCRIPTION

An object describing a sequencing technology, such as SLX or SOLID.
Seq_techs are usually attached to a VRTrack::Library by the
seq_tech_id on the library.

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
  Arg [2]    : seq_tech id
  Example    : my $seq_tech = VRTrack::Seq_tech->new($vrtrack, $id)
  Description: Returns Seq_tech object by seq_tech_id
  Returntype : VRTrack::Seq_tech object

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
    return {seq_tech_id  => sub { $self->id(@_) },
	    name         => sub { $self->name(@_) }};
}


=head2 new_by_name

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : seq_tech name
  Example    : my $seq_tech = VRTrack::Seq_tech->new_by_name($vrtrack, $name)
  Description: Class method. Returns Seq_tech object by seq_tech name.  If no such name is in the database, returns undef
  Returntype : VRTrack::Seq_tech object

=cut


=head2 create

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : seq_tech name
  Example    : my $seq_tech = VRTrack::Seq_tech->create($vrtrack, $name)
  Description: Class method.  Creates new Seq_tech object in the database.
  Returntype : VRTrack::Seq_tech object

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
  Example    : my $id = $seq_tech->id();
               $seq_tech->id('104');
  Description: Get/Set for database ID of a seq_tech
  Returntype : Internal ID integer

=cut


=head2 name

  Arg [1]    : name (optional)
  Example    : my $name = $seq_tech->name();
               $seq_tech->name('SLX');
  Description: Get/Set for seq_tech name
  Returntype : string

=cut


=head2 update

  Arg [1]    : None
  Example    : $seq_tech->update();
  Description: Update a seq_tech whose properties you have changed.  If properties haven't changed (i.e. dirty flag is unset) do nothing.  
               Unsets the dirty flag on success.
  Returntype : 1 if successful, otherwise undef.

=cut

1;
