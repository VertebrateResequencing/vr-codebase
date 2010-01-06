package VRTrack::Hierarchy_obj; 

=head1 NAME

VRTrack::Hierarchy_obj - Sequence Tracking Hierarchy_obj object

=head1 SYNOPSIS

=head1 DESCRIPTION

This is the superclass of objects that we model on our disk hierarchy structure.
They all have a hierarchy_name() method.

It inherits from Table_obj.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;
use Carp qw(cluck confess);

use base qw(VRTrack::Table_obj);


=head2 new

  Arg [1]    : vrtrack handle
  Arg [2]    : obj id
  Example    : my $obj= $class->new($vrtrack, $id)
  Description: Returns core objects by id
  Returntype : $class object

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    return $self;
}


=head2 new_by_hierarchy_name

  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : hierarchy_name
  Example    : my $obj = VRTrack::Hierarchy_obj->new_by_hierarchy_name($vrtrack, $hierarchy_name)
  Description: Class method. Returns latest object by hierarchy_name.  If no
               such hierarchy_name is in the database, returns undef.
	       Dies if multiple hierarchy_names match.
  Returntype : VRTrack::Hierarchy_obj inheriting object

=cut

sub new_by_hierarchy_name {
    my ($class, $vrtrack, $hierarchy_name) = @_;
    confess "Need to call with a vrtrack handle, hierarchy_name" unless ($vrtrack && $hierarchy_name);
    return $class->new_by_field_value($vrtrack, 'hierarchy_name', $hierarchy_name);
}


=head2 hierarchy_name

  Arg [1]    : directory name (optional)
  Example    : my $hname = $project->hierarchy_name();
  Description: Get/set hierarchy name.  This is the directory name
               (without path) that the object will be named in a file hierarchy.
  Returntype : string

=cut

sub hierarchy_name {
    my $self = shift;
    return $self->_get_set('hierarchy_name', 'string', @_);
}

1;
