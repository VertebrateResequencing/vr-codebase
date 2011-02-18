package VRTrack::Exome_design; 

=head1 NAME

VRTrack::Exome_design - Sequence Tracking Exome_design object

=head1 SYNOPSIS
    my $exome_design = VRTrack::Exome_design->new($vrtrack, $exome_design_id);

    my $id = $exome_design->id();
    my $name = $exome_design->name();

=head1 DESCRIPTION

An object describing an exome project's baits and targets

=head1 AUTHOR

mh12@sanger.ac.uk (author)

=head1 METHODS

=cut

use strict;
use warnings;
use Carp qw(cluck confess);

use base qw(VRTrack::Named_obj);


###############################################################################
## Class methods
################################################################################

=head2 new

Arg [1] : vrtrack handle to seqtracking database
Arg [2] : exome_design_id
Example : my $exome_design = VRTrack::Exome_design->new($vrtrack, $id)
Description: Returns Exome_design object by exome_design_id
Returntype : VRTrack::Exome_design object

=cut

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    return $self;
}

=head2 fields_dispatch

  Arg [1]    : none
  Example    : my $fieldsref = $exome_design->fields_dispatch();
  Description: Returns hashref dispatch table keyed on database field
               Used internally for new and update methods
  Returntype : hashref

=cut

sub fields_dispatch {
    my $self = shift;
    return {exome_design_id => sub { $self->id(@_)},
            name            => sub { $self->name(@_)},
            bait_bases      => sub { $self->bait_bases(@_)},
            target_bases    => sub { $self->target_bases(@_)}};
}


=head2 new_by_name

Arg [1] : vrtrack handle to seqtracking database
Arg [2] : exome design name
Example : my $exome_design = VRTrack::Exome_design->new_by_name($vrtrack, $name)
Description: Class method. Returns Exome_design object by name and project_id. If no such name is in the database, returns undef
Returntype : VRTrack::Exome_design object

=cut


=head2 create

Arg [1] : vrtrack handle to seqtracking database
Arg [2] : exome design namee
Example : my $exome_design = VRTrack::Exome_design->create($vrtrack, $name)
Description: Class method. Creates new Exome_design object in the database.
Returntype : VRTrack::Exome_design object

=cut

sub create {
    my ($self, $vrtrack, $value) = @_;
    return $self->SUPER::create($vrtrack, name => $value);
}


###############################################################################
## Object methods
################################################################################


=head2 dirty

Arg [1] : boolean for dirty status
Example : $obj->dirty(1);
Description: Get/Set for object properties having been altered.
Returntype : boolean

=cut


=head2 id

Arg [1] : id (optional)
Example : my $id = $exome_design->id();
          $exome_design->id('104');
Description: Get/Set for database ID of an exome design
Returntype : Internal ID integer

=cut


=head2 name

Arg [1] : name (optional)
Example : my $name = $exome_design->name();
          $exome_design->name('uk10k.agilent.1');
Description: Get/Set forexome design name
Returntype : string

=cut


=head2 bait_bases

Arg [1] : bait_bases (optional)
Example : my $bait_bases = $obj->bait_bases();
          $obj->bait_bases(424242);
Description : Get/Set the number of bases in the bait set
Returntype : integer

=cut

sub bait_bases {
    my $self = shift;
    return $self->_get_set('bait_bases', 'number', @_);
}



=head2 target_bases

Arg [1] : target_bases (optional)
Example : my $target_bases = $obj->target_bases();
          $obj->target_bases(424242);
Description : Get/Set the number of bases in the target set
Returntype : integer

=cut

sub target_bases {
    my $self = shift;
    return $self->_get_set('target_bases', 'number', @_);
}

1;

