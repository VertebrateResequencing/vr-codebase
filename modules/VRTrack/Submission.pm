package VRTrack::Submission;

=head1 NAME

VRTrack::Submission - Sequence Tracking Submission object

=head1 SYNOPSIS
    my $sub = VRTrack::Submission->new($vrtrack, $submission_id);

    my $id      = $submission->id();
    my $date    = $submission->date();
    my $name    = $submission->name();
    my $acc     = $submission->acc();

=head1 DESCRIPTION

An object describing an [ES]RA submission, so we can track which lanes have
been submitted, and to what submission.

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
  Arg [2]    : submission id
  Example    : my $sub = VRTrack::Submission->new($vrtrack, $id)
  Description: Returns Submission object by submission_id
  Returntype : VRTrack::Submission object

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
    return {submission_id  => sub { $self->id(@_) },
	    date           => sub { $self->date(@_) },
	    name           => sub { $self->name(@_) },
	    acc            => sub { $self->acc(@_) }};
}


=head2 new_by_name

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : submission name
  Example    : my $ind = VRTrack::Submission->new($vrtrack, $name)
  Description: Class method. Returns Submission object by name.  If no such name is in the database, returns undef
  Returntype : VRTrack::Submission object

=cut


=head2 create

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : submission name
  Example    : my $ind = VRTrack::Submission->create($vrtrack, $name)
  Description: Class method. Creates new Submission object in the database.
  Returntype : VRTrack::Submission object

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
  Example    : my $id = $sub->id();
               $sub->id('104');
  Description: Get/Set for database ID of a submission
  Returntype : Internal ID integer

=cut


=head2 date

  Arg [1]    : date (optional)
  Example    : my $date = $sub->date();
               $sub->date('20091023');
  Description: Get/Set for submission date
  Returntype : date string

=cut

sub date {
    my $self = shift;
    return $self->_get_set('date', 'string', @_);
}


=head2 acc

  Arg [1]    : acc (optional)
  Example    : my $acc = $sub->acc();
               $sub->acc('ERA000099');
  Description: Get/Set for submission acc
  Returntype : string

=cut

sub acc {
    my $self = shift;
    return $self->_get_set('acc', 'string', @_);
}


=head2 name

  Arg [1]    : name (optional)
  Example    : my $name = $sub->name();
               $sub->name('g1k-sc-20090506);
  Description: Get/Set for submission name.  This is the 'submission_id' field in the submission.xml file.
  Returntype : string

=cut


=head2 update

  Arg [1]    : None
  Example    : $submission->update();
  Description: Update a submission whose properties you have changed.  If properties haven't changed (i.e. dirty flag is unset) do nothing.  
               Unsets the dirty flag on success.
  Returntype : 1 if successful, otherwise undef.

=cut

1;
