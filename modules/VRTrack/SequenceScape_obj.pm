package VRTrack::SequenceScape_obj; 

=head1 NAME

VRTrack::SequenceScape_obj - Sequence Tracking SequenceScape_obj object

=head1 SYNOPSIS

=head1 DESCRIPTION

This is the superclass of objects that have an analogue in the sequencescape
database, thus having an ssid() method to provide the relation.

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


=head2 new_by_ssid

  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : sequencescape id
  Arg [3]    : 'latest'(default)|datestamp string(in the format returned by
               changed()) (optional)
  Example    : my $obj = VRTrack::SequenceScape_obj->new_by_ssid($vrtrack, $ssid);
  Description: Class method. Returns latest object by ssid.  If no such ssid is
               in the database, returns undef.
  Returntype : VRTrack::SequenceScape_obj inheriting object

=cut

sub new_by_ssid {
    my ($class, $vrtrack, $ssid, @extra_args) = @_;
    confess "Need to call with a vrtrack handle, ssid" unless ($vrtrack && $ssid);
    return $class->new_by_field_value($vrtrack, 'ssid', $ssid, @extra_args);
}


=head2 ssid

  Arg [1]    : ssid (optional)
  Example    : my $ssid = $proj->ssid();
               $obj->ssid(104);
  Description: Get/Set for SequenceScape ID
  Returntype : int

=cut

sub ssid {
    my $self = shift;
    return $self->_get_set('ssid', 'number', @_);
}

1;
