package VRTrack::Study;

=head1 NAME

VRTrack::Study - Sequence Tracking Study object

=head1 SYNOPSIS
    my $study = VRTrack::Study->new($vrtrack, $study_id);

    my $id      = $study->id();
    my $acc     = $study->acc();

=head1 DESCRIPTION

An object describing an accessioned study (i.e. an ERA/SRA study).
Studys are usually attached to a VRTrack::Project by study_id.

=head1 AUTHOR

jws@sanger.ac.uk (author)

=head1 METHODS

=cut

use strict;
use warnings;
use Carp qw(cluck confess);

use base qw(VRTrack::Table_obj);


###############################################################################
# Class methods
###############################################################################

=head2 new

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : study id
  Example    : my $study = VRTrack::Study->new($vrtrack, $id)
  Description: Returns Study object by study_id
  Returntype : VRTrack::Study object

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
    return {study_id => sub { $self->id(@_)},
            acc      => sub { $self->acc(@_)}};
}


=head2 new_by_acc

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : study acc
  Example    : my $study = VRTrack::Study->new_by_acc($vrtrack, $acc)
  Description: Class method. Returns Study object by acc.  If no such acc is in the database, returns undef
  Returntype : VRTrack::Study object

=cut

sub new_by_acc {
    my ($class, $vrtrack, $value) = @_;
    return $class->new_by_field_value($vrtrack, 'acc', $value);
}


=head2 create

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : study acc
  Example    : my $study = VRTrack::Study->create($vrtrack, $acc)
  Description: Class method.  Creates new Study object in the database.
  Returntype : VRTrack::Study object

=cut

sub create {
    my ($self, $vrtrack, $value) = @_;
    return $self->SUPER::create($vrtrack, acc => $value);
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
  Example    : my $id = $study->id();
               $study->id('1');
  Description: Get/Set for database ID of a study
  Returntype : Internal ID integer

=cut


=head2 acc

  Arg [1]    : acc (optional)
  Example    : my $acc = $study->acc();
               $samp->acc('ERS000090');
  Description: Get/Set for study accession, i.e. SRA/ERA study id
  Returntype : string

=cut

sub acc {
    my $self = shift;
    return $self->_get_set('acc', 'string', @_);
}


=head2 update

  Arg [1]    : None
  Example    : $study->update();
  Description: Update a study whose properties you have changed.  If properties haven't changed (i.e. dirty flag is unset) do nothing.  
               Unsets the dirty flag on success.
  Returntype : 1 if successful, otherwise undef.

=cut

1;
