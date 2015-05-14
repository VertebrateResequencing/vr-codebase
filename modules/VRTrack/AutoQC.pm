package VRTrack::AutoQC;

=head1 NAME

VRTrack::AutoQC - Sequence Tracking Mapstats AutoQC Test result object

=head1 SYNOPSIS
    my $autoqc = VRTrack::AutoQC->new($vrtrack, $id);

    $autoqc->id('104');

=head1 DESCRIPTION

An object describing the properties of an autoqc test result for a mapstats entry.

=head1 AUTHOR

mp15@sanger.ac.uk (author)

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
  Arg [2]    : autoqc id
  Example    : my $autoqc = VRTrack::AutoQC->new($vrtrack, $id)
  Description: Returns AutoQC object by autoqc_id
  Returntype : VRTrack::AutoQC object

=cut

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    return $self;
}

=head2 fields_dispatch

  Arg [1]    : none
  Example    : my $fieldsref = $lane->fields_dispatch();
  Description: Returns hashref dispatch table keyed on database field
               Used internally for new and update methods
  Returntype : hashref

=cut

sub fields_dispatch {
    my $self = shift;

    my %fields = (
               autoqc_id         => sub { $self->id(@_)},
	           mapstats_id       => sub { $self->mapstats_id(@_) },
               test              => sub { $self->test(@_)},
               result            => sub { $self->result(@_)},
               reason            => sub { $self->reason(@_)},
               );

    return \%fields;
}

=head2 create

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : autoqc test
  Arg [3]    : autoqc result
  Arg [4]    : autoqc reason
  Example    : my $autoqc = VRTrack::AutoQC->new($vrtrack,$test,$result,$reason);
  Description: Class method. Creates new AutoQC object in the database.
  Returntype : VRTrack::AutoQC object

=cut

sub create {
    my ($self, $vrtrack,$test,$result,$reason) = @_;
    return $self->SUPER::create($vrtrack, test => $test, result=> $result, reason => $reason);
}


###############################################################################
# Object methods
###############################################################################

=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $autoqc->id();
               $autoqc->id('104');
  Description: Get/Set for ID of an autoqc object
  Returntype : Internal ID integer

=cut

sub test {
    my $self = shift;
    return $self->_get_set('test', 'string', @_);
}
sub result {
    my $self = shift;
    return $self->_get_set('result', 'number', @_);
}
sub reason {
    my $self = shift;
    return $self->_get_set('reason', 'string', @_);
}
sub mapstats_id {
    my $self = shift;
    return $self->_get_set('mapstats_id', 'number', @_);
}


1;

