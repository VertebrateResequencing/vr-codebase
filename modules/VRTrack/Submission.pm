package VRTrack::Submission;
# author: jws
=head1 NAME

VRTrack::Submission - Sequence Tracking Submission object

=head1 SYNOPSIS
    my $sub = VRTrack::Submission->new($vrtrack, $submission_id);

    my $id      = $submission->id();
    my $date    = $submission->date();
    my $name    = $submission->name();
    my $acc     = $submission->acc();

=head1 DESCRIPTION

An object describing an [ES]RA submission, so we can track which lanes have
been submitted, and to what submission.

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
  Arg [2]    : submission id
  Example    : my $sub = VRTrack::Submission->new($vrtrack, $id)
  Description: Returns Submission object by submission_id
  Returntype : VRTrack::Submission object

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

    my $sql = qq[select submission_id, date, name, acc from submission where submission_id = ?];
    my $sth = $self->{_dbh}->prepare($sql);

    if ($sth->execute($id)){
        my $data = $sth->fetchrow_hashref;
        unless ($data){
            return undef;
        }
        $self->id($data->{'submission_id'});
        $self->date($data->{'date'});
        $self->name($data->{'name'});
        $self->acc($data->{'acc'});
	$self->dirty(0); # unset the dirty flag
    }
    else{
        die(sprintf('Cannot retrieve submission: %s', $DBI::errstr));
    }

    return $self;
}


=head2 new_by_name

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : submission name
  Example    : my $ind = VRTrack::Submission->new($vrtrack, $name)
  Description: Class method. Returns Submission object by name.  If no such name is in the database, returns undef
  Returntype : VRTrack::Submission object

=cut

sub new_by_name {
    my ($class,$vrtrack, $name) = @_;
    die "Need to call with a vrtrack handle and name" unless ($vrtrack && $name);
    if ( $vrtrack->isa('DBI::db') ) { croak "The interface has changed, expected vrtrack reference.\n"; }
    my $dbh = $vrtrack->{_dbh};
    my $sql = qq[select submission_id from submission where name = ?];
    my $sth = $dbh->prepare($sql);

    my $id;
    if ($sth->execute($name)){
        my $data = $sth->fetchrow_hashref;
        unless ($data){
            return undef;
        }
        $id = $data->{'submission_id'};
    }
    else{
        die(sprintf('Cannot retrieve submission by name $name: %s', $DBI::errstr));
    }
    return $class->new($vrtrack, $id);
}


=head2 create

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : submission name
  Example    : my $ind = VRTrack::Submission->create($vrtrack, $name)
  Description: Class method. Creates new Submission object in the database.
  Returntype : VRTrack::Submission object

=cut

sub create {
    my ($class,$vrtrack, $name) = @_;
    die "Need to call with a vrtrack handle and name" unless ($vrtrack && $name);
    if ( $vrtrack->isa('DBI::db') ) { croak "The interface has changed, expected vrtrack reference.\n"; }
    my $dbh = $vrtrack->{_dbh};
    my $sql = qq[INSERT INTO submission (submission_id, name) 
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
  Example    : my $id = $sub->id();
               $sub->id('104');
  Description: Get/Set for database ID of a submission
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


=head2 date

  Arg [1]    : date (optional)
  Example    : my $date = $sub->date();
               $sub->date('20091023');
  Description: Get/Set for submission date
  Returntype : date string

=cut

sub date {
    my ($self,$date) = @_;
    if (defined $date and $date ne $self->{'date'}){
        $self->{'date'} = $date;
	$self->dirty(1);
    }
    return $self->{'date'};
}


=head2 acc

  Arg [1]    : acc (optional)
  Example    : my $acc = $sub->acc();
               $sub->acc('ERA000099');
  Description: Get/Set for submission acc
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


=head2 name

  Arg [1]    : name (optional)
  Example    : my $name = $sub->name();
               $sub->name('g1k-sc-20090506);
  Description: Get/Set for submission name.  This is the 'submission_id' field in the submission.xml file.
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
  Example    : $submission->update();
  Description: Update a submission whose properties you have changed.  If properties haven't changed (i.e. dirty flag is unset) do nothing.  
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
	    my $updsql = qq[UPDATE submission SET date=?, name=?, acc=? WHERE submission_id = ? ];
	    
	    $dbh->do ($updsql, undef, $self->date, $self->name, $self->acc,$self->id);
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
