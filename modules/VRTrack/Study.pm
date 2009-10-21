package VRTrack::Study;
# author: jws
=head1 NAME

VRTrack::Study - Sequence Tracking Study object

=head1 SYNOPSIS
    my $study = VRTrack::Study->new($vrtrack, $study_id);

    my $id      = $study->id();
    my $acc     = $study->acc();

=head1 DESCRIPTION

An object describing an accessioned study (i.e. an ERA/SRA study).
Studys are usually attached to a VRTrack::Project by study_id.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;
use Carp;
no warnings 'uninitialized';

###############################################################################
# Class methods
###############################################################################

=head2 new

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : study id
  Example    : my $study = VRTrack::Study->new($vrtrack, $id)
  Description: Returns Study object by study_id
  Returntype : VRTrack::Study object

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

    my $sql = qq[select study_id, acc from study where study_id = ?];
    my $sth = $self->{_dbh}->prepare($sql);

    if ($sth->execute($id)){
        my $data = $sth->fetchrow_hashref;
        unless ($data){
            return undef;
        }
        $self->id($data->{'study_id'});
        $self->acc($data->{'acc'});
	$self->dirty(0); # unset the dirty flag
    }
    else{
        die(sprintf('Cannot retrieve study: %s', $DBI::errstr));
    }

    return $self;
}


=head2 new_by_acc

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : study acc
  Example    : my $study = VRTrack::Study->new_by_acc($vrtrack, $acc)
  Description: Class method. Returns Study object by acc and project_id.  If no such acc is in the database, returns undef
  Returntype : VRTrack::Study object

=cut

sub new_by_acc {
    my ($class,$vrtrack, $acc) = @_;
    die "Need to call with a vrtrack handle, acc" unless ($vrtrack && $acc);
    if ( $vrtrack->isa('DBI::db') ) { croak "The interface has changed, expected vrtrack reference.\n"; }
    my $dbh = $vrtrack->{_dbh};
    my $sql = qq[select study_id from study where acc = ?];
    my $sth = $dbh->prepare($sql);

    my $id;
    if ($sth->execute($acc)){
        my $data = $sth->fetchrow_hashref;
        unless ($data){
            return undef;
        }
        $id = $data->{'study_id'};
    }
    else{
        die(sprintf('Cannot retrieve study by $acc: %s', $DBI::errstr));
    }
    return $class->new($vrtrack, $id);
}


=head2 create

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : study acc
  Example    : my $study = VRTrack::Study->create($vrtrack, $acc)
  Description: Class method.  Creates new Study object in the database.
  Returntype : VRTrack::Study object

=cut

sub create {
    my ($class,$vrtrack, $acc) = @_;
    die "Need to call with a vrtrack handle and acc" unless ($vrtrack && $acc);
    if ( $vrtrack->isa('DBI::db') ) { croak "The interface has changed, expected vrtrack reference.\n"; }
    my $dbh = $vrtrack->{_dbh};
    my $sql = qq[INSERT INTO study (study_id, acc) 
                 VALUES (NULL,?)];

    my $sth = $dbh->prepare($sql);
    my $id;
    if ($sth->execute( $acc)) {
        $id = $dbh->{'mysql_insertid'};
    }
    else {
        die( sprintf('DB load insert failed: %s %s', $acc, $DBI::errstr));
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
  Example    : my $id = $study->id();
               $study->id('1');
  Description: Get/Set for database ID of a study
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


=head2 acc

  Arg [1]    : acc (optional)
  Example    : my $acc = $study->acc();
               $samp->acc('ERS000090');
  Description: Get/Set for study accession, i.e. SRA/ERA study id
  Returntype : string

=cut

sub acc {
    my ($self,$acc) = @_;
    if (defined $acc and $acc ne $self->{'acc'}){
        $self->{'acc'} = $acc;
	$self->dirty(1);
    }
    return $self->{'acc'};
}


=head2 update

  Arg [1]    : None
  Example    : $study->update();
  Description: Update a study whose properties you have changed.  If properties haven't changed (i.e. dirty flag is unset) do nothing.  
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
	    my $updsql = qq[UPDATE study SET acc=? WHERE study_id = ? ];
	    
	    $dbh->do ($updsql, undef, $self->acc,$self->id);
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
