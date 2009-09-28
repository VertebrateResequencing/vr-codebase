package VRTrack::Image;
# author: jws
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

jws@sanger.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;
use Carp;
no warnings 'uninitialized';

use constant DBI_DUPLICATE => '1062';

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
    my ($class,$vrtrack, $id) = @_;
    die "Need to call with a vrtrack handle and id" unless ($vrtrack && $id);
    if ( $vrtrack->isa('DBI::db') ) { croak "The interface has changed, expected vrtrack reference.\n"; }
    my $self = {};
    bless ($self, $class);
    my $dbh = $vrtrack->{_dbh};
    $self->{vrtrack} = $vrtrack;
    $self->{_dbh} = $dbh;

    my $sql = qq[select image_id, mapstats_id, name, caption, image from image where image_id = ?];
    my $sth = $self->{_dbh}->prepare($sql);

    if ($sth->execute($id)){
        my $data = $sth->fetchrow_hashref;
        unless ($data){
            return undef;
        }
        $self->id($data->{'image_id'});
        $self->mapstats_id($data->{'mapstats_id'});
        $self->name($data->{'name'});
        $self->caption($data->{'caption'});
        $self->image($data->{'image'});
	$self->dirty(0); # unset the dirty flag
    }
    else{
        die(sprintf('Cannot retrieve image: %s', $DBI::errstr));
    }

    return $self;
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
    my ($class,$vrtrack, $name, $img) = @_;
    die "Need to call with a vrtrack handle and name and image" unless ($vrtrack && $name && $img);
    if ( $vrtrack->isa('DBI::db') ) { croak "The interface has changed, expected vrtrack reference.\n"; }
    my $dbh = $vrtrack->{_dbh};
    my $sql = qq[INSERT INTO image (image_id, name, image) 
                 VALUES (NULL,?,?)];

                
    my $sth = $dbh->prepare($sql);
    my $id;
    if ($sth->execute( $name,$img)) {
        $id = $dbh->{'mysql_insertid'};
    }
    else {
        die( sprintf('DB load insert failed: %s %s', $name, $DBI::errstr));
    }
    return $class->new($vrtrack, $id);

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

sub dirty {
    my ($self,$dirty) = @_;
    if (defined $dirty){
	$self->{_dirty} = $dirty ? 1 : 0;
    }
    return $self->{_dirty};
}


=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $image->id();
               $image->id('104');
  Description: Get/Set for ID of a image
  Returntype : Internal ID integer

=cut

sub id {
    my ($self,$id) = @_;
    if (defined $id and $id != $self->{'id'}){
        $self->{'id'} = $id;
	$self->dirty(1);
    }
    return $self->{'id'};
}


=head2 name

  Arg [1]    : name (optional)
  Example    : my $name = $image->name();
               $image->name('gc.png');
  Description: Get/Set for image name
  Returntype : string

=cut

sub name {
    my ($self,$name) = @_;
    if (defined $name and $name ne $self->{'name'}){
        $self->{'name'} = $name;
	$self->dirty(1);
    }
    return $self->{'name'};
}


=head2 caption

  Arg [1]    : caption (optional)
  Example    : my $caption = $image->caption();
               $image->caption('percent GC');
  Description: Get/Set for image caption
  Returntype : string

=cut

sub caption {
    my ($self,$caption) = @_;
    if (defined $caption and $caption ne $self->{'caption'}){
        $self->{'caption'} = $caption;
	$self->dirty(1);
    }
    return $self->{'caption'};
}


=head2 image

  Arg [1]    : image (optional)
  Example    : my $img = $image->image();
               $image->image($my_image);
  Description: Get/Set for actual image data (i.e. the image itself)
  Returntype : binary image

=cut

sub image {
    my ($self,$image) = @_;
    if (defined $image and $image ne $self->{'image'}){
        $self->{'image'} = $image;
	$self->dirty(1);
    }
    return $self->{'image'};
}


=head2 mapstats_id

  Arg [1]    : mapstats_id (optional)
  Example    : my $mapstats_id = $image->mapstats_id();
               $image->mapstats_id(11);
  Description: Get/Set for image internal DB mapstats_id
  Returntype : integer

=cut

sub mapstats_id {
    my ($self,$mapstats_id) = @_;
    if (defined $mapstats_id and $mapstats_id != $self->{'mapstats_id'}){
        $self->{'mapstats_id'} = $mapstats_id;
	$self->dirty(1);
    }
    return $self->{'mapstats_id'};
}


=head2 update

  Arg [1]    : None
  Example    : $image->update();
  Description: Update a image whose properties you have changed.  If properties haven't changed (i.e. dirty flag is unset) do nothing.  
	       Changes the changed datestamp to now() on the mysql server (i.e. you don't have to set changed yourself, and indeed if you do, it will be overridden).
               Unsets the dirty flag on success.
  Returntype : 1 if successful, otherwise undef.

=cut

sub update {
    my ($self) = @_;
    my $success = undef;
    if ($self->dirty){
	my $dbh = $self->{_dbh};
	my $save_re = $dbh->{RaiseError};
	my $save_pe = $dbh->{PrintError};
	$dbh->{RaiseError} = 1; # raise exception if an error occurs
	$dbh->{PrintError} = 0; # don't print an error message

	eval {
	    my $updsql = qq[UPDATE image SET mapstats_id=?, name=?, caption=?, image=? WHERE image_id = ? ];
	    
	    $dbh->do ($updsql, undef,$self->mapstats_id, $self->name, $self->caption, $self->image, $self->id);
	};

	if (!$@) {
	    $success = 1;
	}

	# restore attributes to original state
	$dbh->{PrintError} = $save_pe;
	$dbh->{RaiseError} = $save_re;

    }
    if ($success){
        $self->dirty(0);
    }

    return $success;
}

1;
