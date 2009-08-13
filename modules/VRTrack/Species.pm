package VRTrack::Species;
=head1 NAME

VRTrack::Species - Sequence Tracking Species object

=head1 SYNOPSIS
    my $spp = VRTrack::Species->new($dbh, $species_id);

    my $id      = $species->id();
    my $name    = $species->name();
    my $taxon   = $species->taxon_id();

=head1 DESCRIPTION

An object describing a species.  Species are usually attached to a
VRTrack::Sample by the species_id on the individual.

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
  Arg [2]    : species id
  Example    : my $spp = VRTrack::Species->new($dbh, $id)
  Description: Returns Species object by species_id
  Returntype : VRTrack::Species object

=cut

sub new {
    my ($class,$dbh, $id) = @_;
    die "Need to call with a db handle and id" unless ($dbh && $id);
    my $self = {};
    bless ($self, $class);
    $self->{_dbh} = $dbh;

    my $sql = qq[select species_id, name, taxon_id from species where species_id = ?];
    my $sth = $self->{_dbh}->prepare($sql);

    if ($sth->execute($id)){
        my $data = $sth->fetchrow_hashref;
        unless ($data){
            return undef;
        }
        $self->id($data->{'species_id'});
        $self->name($data->{'name'});
        $self->taxon_id($data->{'taxon_id'});
	$self->dirty(0); # unset the dirty flag
    }
    else{
        die(sprintf('Cannot retrieve species: %s', $DBI::errstr));
    }

    return $self;
}


=head2 new_by_name

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : species name
  Example    : my $spp = VRTrack::Species->new($dbh, $name)
  Description: Class method. Returns Species object by name.  If no such name is in the database, returns undef
  Returntype : VRTrack::Species object

=cut

sub new_by_name {
    my ($class,$dbh, $name) = @_;
    die "Need to call with a db handle and name" unless ($dbh && $name);
    my $sql = qq[select species_id from species where name = ?];
    my $sth = $dbh->prepare($sql);

    my $id;
    if ($sth->execute($name)){
        my $data = $sth->fetchrow_hashref;
        unless ($data){
            return undef;
        }
        $id = $data->{'species_id'};
    }
    else{
        die(sprintf('Cannot retrieve species by name $name: %s', $DBI::errstr));
    }
    return $class->new($dbh, $id);
}


=head2 create

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : species name
  Example    : my $spp = VRTrack::Species->create($dbh, $name)
  Description: Class method. Creates new Species object in the database.
  Returntype : VRTrack::Species object

=cut

sub create {
    my ($class,$dbh, $name) = @_;
    die "Need to call with a db handle and name" unless ($dbh && $name);

    my $sql = qq[INSERT INTO species (species_id, name) 
                 VALUES (NULL,?)];

                
    my $sth = $dbh->prepare($sql);
    my $id;
    if ($sth->execute( $name)) {
        $id = $dbh->{'mysql_insertid'};
    }
    else {
        die( sprintf('DB load insert failed: %s %s', $name, $DBI::errstr));
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
  Example    : my $id = $spp->id();
               $spp->id('104');
  Description: Get/Set for database ID of a species
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
  Example    : my $name = $spp->name();
               $spp->name('Homo sapiens');
  Description: Get/Set for species name
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


=head2 taxon_id

  Arg [1]    : taxon_id (optional)
  Example    : my $taxon_id = $spp->taxon_id();
               $spp->taxon_id(1054);
  Description: Get/Set for species taxon id
  Returntype : integer

=cut

sub taxon_id {
    my ($self,$taxon_id) = @_;
    if (defined $taxon_id and $taxon_id ne $self->{'taxon_id'}){
        $self->{'taxon_id'} = $taxon_id;
	$self->dirty(1);
    }
    return $self->{'taxon_id'};
}


=head2 update

  Arg [1]    : None
  Example    : $species->update();
  Description: Update a species whose properties you have changed.  If properties haven't changed (i.e. dirty flag is unset) do nothing.  
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
	    my $updsql = qq[UPDATE species SET name=?, taxon_id=? WHERE species_id = ? ];
	    
	    $dbh->do ($updsql, undef, $self->name, $self->taxon_id, $self->id);
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
