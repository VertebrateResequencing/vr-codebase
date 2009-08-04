package Sfind::Project;
=head1 NAME

Sfind::Project - Sequence Tracking Project object

=head1 SYNOPSIS
    my $proj = Sfind::Project->new($dbh, $project_id);

    #get arrayref of sample objects in a project
    my $samples = $project->samples();
    
    my $id = $project->id();
    my $name = $project->name();

=head1 DESCRIPTION

An object describing the tracked properties of a project.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;
no warnings 'uninitialized';
use Sfind::Sample;


=head2 new

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : project id
  Example    : my $proj = Sfind::Sfind->new($dbh, $id)
  Description: Returns Project object by project_id
  Returntype : Sfind::Project object

=cut

sub new {
    my ($class,$dbh, $id) = @_;
    die "Need to call with a db handle and id" unless ($dbh && $id);
    my $self = {};
    bless ($self, $class);
    $self->{_dbh} = $dbh;

    my $sql = qq[select project_name from project_information where project_id=?];
    my $id_ref = $self->{_dbh}->selectrow_hashref($sql, undef, ($id));
    if ($id_ref){
	my $name = $id_ref->{project_name};
	#warn "Project name : $name\n";
	$self->id($id);
	$self->name($name);
    }
    else {
	return undef;
    }

    return $self;
}



=head2 samples

  Arg [1]    : None
  Example    : my $samples = $project->samples();
  Description: Returns a ref to an array of the sample objects that are associated with this project
  Returntype : ref to array of Sfind::Sample objects

=cut

sub samples {
    my ($self) = @_;

    unless ($self->{'samples'}){
	my @samples;
    	foreach my $id (@{$self->sample_ids()}){
	    my $obj = Sfind::Sample->new($self->{_dbh},$id, $self->id);
	    push @samples, $obj;
	}
	$self->{'samples'} = \@samples;
    }

    return $self->{'samples'};
}


=head2 sample_ids

  Arg [1]    : None
  Example    : my $sample_ids = $project->sample_ids();
  Description: Returns a ref to an array of the sample IDs that are associated with this project
  Returntype : ref to sorted array of integer sample IDs

=cut

sub sample_ids {
    my ($self) = @_;

    unless ($self->{'sample_ids'}){
	my $sql = qq[select distinct(sample_id) from requests where project_id=?];
	my @samples;
	my $sth = $self->{_dbh}->prepare($sql);

	$sth->execute($self->id);
	foreach(@{$sth->fetchall_arrayref()}){
	    push @samples, $_->[0];
	}
	@samples = sort {$a <=> $b} @samples;
	$self->{'sample_ids'} = \@samples;
    }
 
    return $self->{'sample_ids'};
}


=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $proj->id();
	       $proj->id('104');
  Description: Get/Set for ID of a project
  Returntype : SequenceScape ID (usu. integer)

=cut

sub id {
    my ($self,$id) = @_;
    if ($id){
	$self->{'id'} = $id;
    }
    return $self->{'id'};
}


=head2 name

  Arg [1]    : name (optional)
  Example    : my $name = $proj->name();
	       $proj->name('1000Genomes-A1-CEU');
  Description: Get/Set for project name
  Returntype : string

=cut

sub name {
    my ($self,$name) = @_;
    if ($name){
	$self->{'name'} = $name;
    }
    return $self->{'name'};
}


=head2 get_sample_by_id

  Arg [1]    : sample id from sequencescape
  Example    : my $sample = $proj->get_sample_by_id(1154);
  Description: retrieve sample object by sequencescape id
  Returntype : Sfind::Sample object

=cut

sub get_sample_by_id {
    my ($self, $id) = @_;
    my $obj = Sfind::Sample->new($self->{_dbh},$id,$self->id);
    return $obj;
}


=head2 get_sample_by_name

  Arg [1]    : sample name from sequencescape
  Example    : my $sample = $proj->get_sample_by_name('NA12878')
  Description: retrieve sample object by sequencescape name
  Returntype : Sfind::Sample object

=cut

sub get_sample_by_name {
    my ($self, $name) = @_;
    my $sql = qq[select distinct(sample_id) from requests where sample_name=? and project_id=?];
    my $id_ref = $self->{_dbh}->selectrow_hashref($sql, undef, ($name, $self->id));
    unless ($id_ref){
	warn "No sample with name $name\n";
	return undef;
    }

    my $id = $id_ref->{sample_id};
    return $self->get_sample_by_id($id);
}
1;
