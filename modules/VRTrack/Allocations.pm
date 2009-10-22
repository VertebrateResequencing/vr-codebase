package VRTrack::Allocations;
# author: jws
=head1 NAME

VRTrack::Allocations - Sequence Tracking Allocations object

=head1 SYNOPSIS
    my $allocs = VRTrack::Allocations->new($vrtrack);
    $allocs->add_allocation($study,$individual,$seq_centre);
    $allocs->get_study_centre_for_ind($individual);

=head1 DESCRIPTION

An object for managing sample/study/centre allocations.

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
  Example    : my $allocations = VRTrack::Allocations->new($vrtrack)
  Description: Returns Allocations object
  Returntype : VRTrack::Allocations object

=cut

sub new {
    my ($class,$vrtrack) = @_;
    die "Need to call with a vrtrack handle" unless $vrtrack;
    if ( $vrtrack->isa('DBI::db') ) { croak "The interface has changed, expected vrtrack reference.\n"; }
    my $dbh = $vrtrack->{_dbh};
    my $self = {};
    bless ($self, $class);
    $self->{_dbh} = $dbh;
    $self->{vrtrack} = $vrtrack;

    return $self;
}


=head2 add_allocation

  Arg [1]    : study id
  Arg [2]    : individual id
  Arg [3]    : seq centre id
  Example    : $vrallocs->add_allocation(5,3,1);
  Description: Adds an allocation to the database.  Checks if the passed ids exist, and fails if they aren't in the database - you should check this first and create new entries if required.
  Returntype : VRTrack::Allocations object

=cut

sub add_allocation {
    my ($self,$study,$ind,$centre) = @_;
    my $dbh = $self->{_dbh};
    my $vrtrack = $self->{vrtrack};
    my $checkstudy = VRTrack::Study->new($vrtrack,$study);
    my $checkind = VRTrack::Individual->new($vrtrack,$ind);
    my $checkcentre = VRTrack::Seq_centre->new($vrtrack,$centre);

    unless ($checkstudy && $checkind && $checkcentre){
        return 0;
    }
    
    my $sql = qq[INSERT INTO allocation (study_id, 
                                        individual_id, 
                                        seq_centre_id) 
                 VALUES (?,?,?)];
    my $sth = $dbh->prepare($sql);
    my $success = 0;
    if ($sth->execute($study,$ind,$centre )) {
        $success = 1;
    }
    else {
        die( sprintf('DB allocation insert failed: %s %s %s %s', $study,$ind,$centre, $DBI::errstr));
    }
    return $success;
}


=head2 is_allocation_in_database

  Arg [1]    : study id
  Arg [2]    : individual id
  Arg [3]    : seq centre id
  Example    : $vrallocs->is_allocation_in_database(5,3,1);
  Description: Checks to see if the allocation is already present in the database.
  Returntype : boolean

=cut

sub is_allocation_in_database {
    my ($self,$study,$ind,$centre) = @_;
    my $dbh = $self->{_dbh};

    my $sql = qq[SELECT study_id, individual_id, seq_centre_id from allocation  where study_id = ? and individual_id = ? and seq_centre_id = ?];
    my $sth = $dbh->prepare($sql);
    my $already_used = 0;
    if ($sth->execute($study,$ind,$centre )) {
        my $data = $sth->fetchrow_hashref;
        if ($data){
            $already_used = 1;
        }
    }
    else {
        die( sprintf('DB allocation retrieval failed: %s %s %s %s', $study,$ind,$centre, $DBI::errstr));
    }
    return $already_used;
}


=head2 get_centres_for_study_ind

  Arg [1]    : study id
  Arg [2]    : individual id
  Example    : my @centres = @{$vrallocs->get_centres_for_study_ind(5,3)};
  Description: get a list of sequencing centres that an individual has been allocated to in a specific study.
  Returntype : arrayref of VRTrack::Seq_centres

=cut

sub get_centres_for_study_ind {
    my ($self,$study,$ind) = @_;
    my $dbh = $self->{_dbh};
    my $vrtrack = $self->{vrtrack};

    my @centres;
    my $sql = qq[SELECT seq_centre_id from allocation 
                 where study_id = ? and individual_id = ?];
    my $sth = $dbh->prepare($sql);
    if ($sth->execute($study,$ind)) {
        foreach(@{$sth->fetchall_arrayref()}){
            push @centres, VRTrack::Seq_centre->new($vrtrack,$_->[0]);
        } 
    }
    else {
        die( sprintf('DB allocation retrieval failed: %s %s %s', $study,$ind, $DBI::errstr));
    }
    return \@centres;
}



1;

