package VRTrack::Mapper;

=head1 NAME

VRTrack::Mapper - Sequence Tracking Mapper object

=head1 SYNOPSIS
    my $mapper = VRTrack::Mapper->new($vrtrack, $mapper_id);

    my $id      = $mapper->id();
    my $name    = $mapper->name();

=head1 DESCRIPTION

An object describing a sequence mapper, i.e. an aligner such as Maq or BWA.
Mappers are usually attached to a VRTrack::Mapstats object by the
mapper_id on the mapping.

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
  Arg [2]    : mapper id
  Example    : my $mapper = VRTrack::Mapper->new($vrtrack, $id)
  Description: Returns Mapper object by mapper_id
  Returntype : VRTrack::Mapper object

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
    return {mapper_id  => sub { $self->id(@_) },
	    name       => sub { $self->name(@_) },
	    version    => sub { $self->version(@_) }};
}


=head2 new_by_name_version

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : mapper name
  Arg [3]    : mapper version
  Example    : my $ind = VRTrack::Mapper->new($vrtrack, $name, $v)
  Description: Class method. Returns Mapper object by name.  If no such name is in the database, returns undef
  Returntype : VRTrack::Mapper object

=cut

sub new_by_name_version {
    my ($class, $vrtrack, $name, $version) = @_;
    confess "Need to call with a vrtrack handle and name" unless ($vrtrack && $name && $version);
    confess "The interface has changed, expected vrtrack reference." if $vrtrack->isa('DBI::db');
    
    my $dbh = $vrtrack->{_dbh};
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
        confess(sprintf('Cannot retrieve mapper by name %s version %s: %s', $name, $version,$DBI::errstr));
    }
    
    return $class->new($vrtrack, $id);
}


=head2 create

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : mapper name
  Arg [3]    : mapper version
  Example    : my $ind = VRTrack::Mapper->create($vrtrack, $name, $v)
  Description: Class method. Creates new Mapper object in the database.
  Returntype : VRTrack::Mapper object

=cut

sub create {
    my ($self, $vrtrack, $name, $version) = @_;
    return $self->SUPER::create($vrtrack, name => $name, version => $version);
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
  Example    : my $id = $mapper->id();
               $mapper->id('104');
  Description: Get/Set for database ID of a mapper
  Returntype : Internal ID integer

=cut


=head2 name

  Arg [1]    : name (optional)
  Example    : my $name = $mapper->name();
               $mapper->name('SLX');
  Description: Get/Set for mapper name
  Returntype : string

=cut

sub name {
    # (we can't be a Named_obj because we don't allow new_by_name)
    my $self = shift;
    return $self->_get_set('name', 'string', @_);
}


=head2 version

  Arg [1]    : version (optional)
  Example    : my $version = $mapper->version();
               $mapper->version('1.01a');
  Description: Get/Set for mapper version
  Returntype : string

=cut

sub version {
    my $self = shift;
    return $self->_get_set('version', 'string', @_);
}


=head2 update

  Arg [1]    : None
  Example    : $mapper->update();
  Description: Update a mapper whose properties you have changed.  If properties haven't changed (i.e. dirty flag is unset) do nothing.  
               Unsets the dirty flag on success.
  Returntype : 1 if successful, otherwise undef.

=cut

1;
