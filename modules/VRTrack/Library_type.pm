package VRTrack::Library_type;
# author: jws
=head1 NAME

VRTrack::Library_type - Sequence Tracking Library_type object

=head1 SYNOPSIS
    my $library_type = VRTrack::Library_type->new($vrtrack, $library_type_id);

    my $id      = $library_type->id();
    my $name    = $library_type->name();

=head1 DESCRIPTION

An object describing a library_type, such as DSS or NOPCR.
Library_types are usually attached to a VRTrack::Library by the
library_type_id on the library.

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
  Arg [2]    : library_type id
  Example    : my $library_type = VRTrack::Library_type->new($vrtrack, $id)
  Description: Returns Library_type object by library_type_id
  Returntype : VRTrack::Library_type object

=cut

sub new {
    my ($class,$vrtrack, $id) = @_;
    die "Need to call with a vrtrack handle and id" unless ($vrtrack && $id);
    if ( $vrtrack->isa('DBI::db') ) { croak "The interface has changed, expected vrtrack reference.\n"; }
    my $dbh = $vrtrack->{_dbh};
    my $self = {};
    bless ($self, $class);
    $self->{_dbh} = $dbh;
    $self->{vrtrack} = $vrtrack;

    my $sql = qq[select library_type_id, name from library_type where library_type_id = ?];
    my $sth = $self->{_dbh}->prepare($sql);

    if ($sth->execute($id)){
        my $data = $sth->fetchrow_hashref;
        unless ($data){
            return undef;
        }
        $self->id($data->{'library_type_id'});
        $self->name($data->{'name'});
	$self->dirty(0); # unset the dirty flag
    }
    else{
        die(sprintf('Cannot retrieve library_type: %s', $DBI::errstr));
    }

    return $self;
}


=head2 new_by_name

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : library_type name
  Example    : my $library_type = VRTrack::Library_type->new_by_name($vrtrack, $name)
  Description: Class method. Returns Library_type object by name and project_id.  If no such name is in the database, returns undef
  Returntype : VRTrack::Library_type object

=cut

sub new_by_name {
    my ($class,$vrtrack, $name) = @_;
    die "Need to call with a vrtrack handle, name" unless ($vrtrack && $name);
    if ( $vrtrack->isa('DBI::db') ) { croak "The interface has changed, expected vrtrack reference.\n"; }
    my $dbh = $vrtrack->{_dbh};
    my $sql = qq[select library_type_id from library_type where name = ?];
    my $sth = $dbh->prepare($sql);

    my $id;
    if ($sth->execute($name)){
        my $data = $sth->fetchrow_hashref;
        unless ($data){
            return undef;
        }
        $id = $data->{'library_type_id'};
    }
    else{
        die(sprintf('Cannot retrieve library_type by $name: %s', $DBI::errstr));
    }
    return $class->new($vrtrack, $id);
}


=head2 create

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : library_type name
  Example    : my $library_type = VRTrack::Library_type->create($vrtrack, $name)
  Description: Class method.  Creates new Library_type object in the database.
  Returntype : VRTrack::Library_type object

=cut

sub create {
    my ($class,$vrtrack, $name) = @_;
    die "Need to call with a vrtrack handle and name" unless ($vrtrack && $name);
    if ( $vrtrack->isa('DBI::db') ) { croak "The interface has changed, expected vrtrack reference.\n"; }
    my $dbh = $vrtrack->{_dbh};
    my $sql = qq[INSERT INTO library_type (library_type_id, name) 
                 VALUES (NULL,?)];

                
    my $sth = $dbh->prepare($sql);
    my $id;
    if ($sth->execute( $name)) {
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
  Example    : my $id = $library_type->id();
               $library_type->id('104');
  Description: Get/Set for database ID of a library_type
  Returntype : Internal ID integer

=cut

sub id {
    my ($self,$id) = @_;
    if (defined $id and $id ne $self->{'id'}){
        $self->{'id'} = $id;
	$self->dirty(1);
    }
    return $self->{'id'};
}


=head2 name

  Arg [1]    : name (optional)
  Example    : my $name = $library_type->name();
               $library_type->name('CEU');
  Description: Get/Set for library_type name
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


=head2 update

  Arg [1]    : None
  Example    : $library_type->update();
  Description: Update a library_type whose properties you have changed.  If properties haven't changed (i.e. dirty flag is unset) do nothing.  
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
	    my $updsql = qq[UPDATE library_type SET name=? WHERE library_type_id = ? ];
	    
	    $dbh->do ($updsql, undef, $self->name,$self->id);
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
