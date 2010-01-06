package VRTrack::Image;

=head1 NAME

VRTrack::Image - Sequence Tracking Image object

=head1 SYNOPSIS
    my $image = VRTrack::Image->new($vrtrack, $image_id);

    my $id = $image->id();
    my $name = $image->name();
    my $img_bin = $image->image();

=head1 DESCRIPTION

An object describing an image file for use on the QC web displays.  The images
are generated as part of QC or as a mapping metric, and as such are properties
of a mapping.  The images is stored in the database, and retrieved as the 'image' property of the object.

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
  Arg [2]    : image id
  Example    : my $image = VRTrack::Image->new($vrtrack, $id)
  Description: Returns Image object by image_id
  Returntype : VRTrack::Image object

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
    return {image_id    => sub { $self->id(@_) },
	    mapstats_id => sub { $self->mapstats_id(@_) },
	    name        => sub { $self->name(@_) },
	    caption     => sub { $self->caption(@_) },
	    image       => sub { $self->image(@_) }};
}


=head2 create

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : image name
  Arg [3]    : image data - the actual image itself
  Example    : my $ind = VRTrack::Image->create($vrtrack, $name, $img)
  Description: Class method. Creates new Image object in the database.
  Returntype : VRTrack::Image object

=cut

sub create {
    my ($self, $vrtrack, $name, $img) = @_;
    return $self->SUPER::create($vrtrack, name => $name, image => $img);
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
  Example    : my $id = $image->id();
               $image->id('104');
  Description: Get/Set for ID of a image
  Returntype : Internal ID integer

=cut


=head2 name

  Arg [1]    : name (optional)
  Example    : my $name = $image->name();
               $image->name('gc.png');
  Description: Get/Set for image name
  Returntype : string

=cut


=head2 caption

  Arg [1]    : caption (optional)
  Example    : my $caption = $image->caption();
               $image->caption('percent GC');
  Description: Get/Set for image caption
  Returntype : string

=cut

sub caption {
    my $self = shift;
    return $self->_get_set('caption', 'string', @_);
}


=head2 image

  Arg [1]    : image (optional)
  Example    : my $img = $image->image();
               $image->image($my_image);
  Description: Get/Set for actual image data (i.e. the image itself)
  Returntype : binary image

=cut

sub image {
    my $self = shift;
    return $self->_get_set('image', 'string', @_);
}


=head2 mapstats_id

  Arg [1]    : mapstats_id (optional)
  Example    : my $mapstats_id = $image->mapstats_id();
               $image->mapstats_id(11);
  Description: Get/Set for image internal DB mapstats_id
  Returntype : integer

=cut

sub mapstats_id {
    my $self = shift;
    return $self->_get_set('mapstats_id', 'number', @_);
}


=head2 update

  Arg [1]    : None
  Example    : $image->update();
  Description: Update a image whose properties you have changed.  If properties haven't changed (i.e. dirty flag is unset) do nothing.  
	       Changes the changed datestamp to now() on the mysql server (i.e. you don't have to set changed yourself, and indeed if you do, it will be overridden).
               Unsets the dirty flag on success.
  Returntype : 1 if successful, otherwise undef.

=cut

1;
