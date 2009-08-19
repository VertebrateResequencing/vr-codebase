package VRTrack::Mapper;
# author: jws
=head1 NAME

VRTrack::Mapper - Sequence Tracking Mapper object

=head1 SYNOPSIS
    my $mapper = VRTrack::Mapper->new($dbh, $mapper_id);

    my $id      = $mapper->id();
    my $name    = $mapper->name();

=head1 DESCRIPTION

An object describing a sequence mapper, i.e. an aligner such as Maq or BWA.
Mappers are usually attached to a VRTrack::Mapstats object by the
mapper_id on the mapping.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use constant DBI_DUPLICATE => '1062';

###############################################################################
# Class methods
###############################################################################

=head2 new

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : mapper id
  Example    : my $mapper = VRTrack::Mapper->new($dbh, $id)
  Description: Returns Mapper object by mapper_id
  Returntype : VRTrack::Mapper object

=cut

sub new {
    my ($class,$dbh, $id) = @_;
    die "Need to call with a db handle and id" unless ($dbh && $id);
    my $self = {};
    bless ($self, $class);
    $self->{_dbh} = $dbh;

    my $sql = qq[select mapper_id, name, version from mapper where mapper_id = ?];
    my $sth = $self->{_dbh}->prepare($sql);

    if ($sth->execute($id)){
        my $data = $sth->fetchrow_hashref;
        unless ($data){
            return undef;
        }
        $self->id($data->{'mapper_id'});
        $self->name($data->{'name'});
        $self->version($data->{'version'});
	$self->dirty(0); # unset the dirty flag
    }
    else{
        die(sprintf('Cannot retrieve mapper: %s', $DBI::errstr));
    }

    return $self;
}


=head2 new_by_name_version

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : mapper name
  Arg [3]    : mapper version
  Example    : my $ind = VRTrack::Mapper->new($dbh, $name, $v)
  Description: Class method. Returns Mapper object by name.  If no such name is in the database, returns undef
  Returntype : VRTrack::Mapper object

=cut

sub new_by_name_version {
    my ($class,$dbh, $name, $version) = @_;
    die "Need to call with a db handle and name" unless ($dbh && $name && $version);
    my $sql = qq[select mapper_id from mapper where name = ? and version = ?];
    my $sth = $dbh->prepare($sql);

    my $id;
    if ($sth->execute($name, $version)){
        my $data = $sth->fetchrow_hashref;
        unless ($data){
            return undef;
        }
        $id = $data->{'mapper_id'};
    }
    else{
        die(sprintf('Cannot retrieve mapper by name %s version %s: %s', $name, $version,$DBI::errstr));
    }
    return $class->new($dbh, $id);
}


=head2 create

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : mapper name
  Arg [3]    : mapper version
  Example    : my $ind = VRTrack::Mapper->create($dbh, $name, $v)
  Description: Class method. Creates new Mapper object in the database.
  Returntype : VRTrack::Mapper object

=cut

sub create {
    my ($class,$dbh, $name, $version) = @_;
    die "Need to call with a db handle, name, version" unless ($dbh && $name && $version);

    my $sql = qq[INSERT INTO mapper (mapper_id, name, version) 
                 VALUES (NULL,?,?)];

                
    my $sth = $dbh->prepare($sql);
    my $id;
    if ($sth->execute( $name,$version)) {
        $id = $dbh->{'mysql_insertid'};
    }
    else {
        die( sprintf('DB load insert failed for %s %s: %s', $name, $version,$DBI::errstr));
    }
    return $class->new($dbh, $id);

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
  Example    : my $id = $mapper->id();
               $mapper->id('104');
  Description: Get/Set for database ID of a mapper
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
  Example    : my $name = $mapper->name();
               $mapper->name('SLX');
  Description: Get/Set for mapper name
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


=head2 version

  Arg [1]    : version (optional)
  Example    : my $version = $mapper->version();
               $mapper->version('1.01a');
  Description: Get/Set for mapper version
  Returntype : string

=cut

sub version {
    my ($self,$version) = @_;
    if (defined $version and $version ne $self->{'version'}){
        $self->{'version'} = $version;
	$self->dirty(1);
    }
    return $self->{'version'};
}


=head2 update

  Arg [1]    : None
  Example    : $mapper->update();
  Description: Update a mapper whose properties you have changed.  If properties haven't changed (i.e. dirty flag is unset) do nothing.  
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
	my $save_ac = $dbh->{AutoCommit};
	$dbh->{RaiseError} = 1; # raise exception if an error occurs
	$dbh->{PrintError} = 0; # don't print an error message
	$dbh->{AutoCommit} = 0; # disable auto-commit

	eval {
	    my $updsql = qq[UPDATE mapper SET name=?,version=? WHERE mapper_id = ? ];
	    
	    $dbh->do ($updsql, undef, $self->name,$self->version,$self->id);
	    $dbh->commit ( );
	};

	if ($@) {
	    warn "Transaction failed, rolling back. Error was:\n$@\n";
	    # roll back within eval to prevent rollback
	    # failure from terminating the script
	    eval { $dbh->rollback ( ); };
	}
	else {
	    $success = 1;
	}

	# restore attributes to original state
	$dbh->{AutoCommit} = $save_ac;
	$dbh->{PrintError} = $save_pe;
	$dbh->{RaiseError} = $save_re;

    }
    if ($success){
        $self->dirty(0);
    }

    return $success;
}

1;
