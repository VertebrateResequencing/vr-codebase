package VRTrack::Project;
# author: jws
=head1 NAME

VRTrack::Project - Sequence Tracking Project object

=head1 SYNOPSIS
    my $proj = VRTrack::Project->new($dbh, $project_id);

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
use VRTrack::Core_obj;
use VRTrack::Sample;
our @ISA = qw(VRTrack::Core_obj);

=head2 fields_dispatch

  Arg [1]    : none
  Example    : my $fieldsref = $proj->fields_dispatch();
  Description: Returns hashref dispatch table keyed on database field
               Used internally for new and update methods
  Returntype : hashref

=cut

sub fields_dispatch {
    my $self = shift;
    my %fields = ( 
                'project_id'        => sub { $self->id(@_)},
                'ssid'              => sub { $self->ssid(@_)},
                'name'              => sub { $self->name(@_)},
                'hierarchy_name'    => sub { $self->hierarchy_name(@_)},
                'acc'               => sub { $self->acc(@_)},
                'note_id'           => sub { $self->note_id(@_)},
                'changed'           => sub { $self->changed(@_)},
                'latest'            => sub { $self->is_latest(@_)},
                );

    return \%fields;
}

###############################################################################
# Class methods
###############################################################################


=head2 new_by_name

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : project name
  Example    : my $project = VRTrack::Project->new_by_name($dbh, $name);
  Description: Class method. Returns latest Project object by name and project_id.  If no such name is in the database, returns undef
  Returntype : VRTrack::Project object

=cut

sub new_by_name {
    my ($class,$dbh, $name) = @_;
    die "Need to call with a db handle, name" unless ($dbh && $name);
    return $class->new_by_field_value($dbh, 'name',$name);
}


=head2 new_by_hierarchy_name

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : project hierarchy_name
  Example    : my $project = VRTrack::Project->new_by_hierarchy_name($dbh, $hierarchy_name)
  Description: Class method. Returns latest Project object by hierarchy_name.  If no such hierarchy_name is in the database, returns undef.  Dies if multiple hierarchy_names match.
  Returntype : VRTrack::Project object

=cut

sub new_by_hierarchy_name {
    my ($class,$dbh, $hierarchy_name) = @_;
    die "Need to call with a db handle, hierarchy_name" unless ($dbh && $hierarchy_name);
    return $class->new_by_field_value($dbh, 'hierarchy_name',$hierarchy_name);
}


=head2 new_by_ssid

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : project sequencescape id
  Example    : my $project = VRTrack::Project->new_by_ssid($dbh, $ssid);
  Description: Class method. Returns latest Project object by ssid.  If no such ssid is in the database, returns undef
  Returntype : VRTrack::Project object

=cut

sub new_by_ssid {
    my ($class,$dbh, $ssid) = @_;
    die "Need to call with a db handle, ssid" unless ($dbh && $ssid);
    return $class->new_by_field_value($dbh, 'ssid',$ssid);
}


=head2 create

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : project name
  Example    : my $project = VRTrack::Project->create($dbh, $name)
  Description: Class method.  Creates new Project object in the database.
  Returntype : VRTrack::Project object

=cut

sub create {
    my ($class,$dbh, $name) = @_;
    die "Need to call with a db handle and name" unless ($dbh && $name);

    my $hierarchy_name = $name;
    $hierarchy_name =~ s/\W+/_/g;

    # prevent adding a project with an existing name
    if ($class->is_name_in_database($dbh, $name, $hierarchy_name)){
        die "Already a project by name $name/$hierarchy_name";
    }

    # lock table, get max id, increment, and use as id
    $dbh->do (qq[LOCK TABLE project WRITE]);
    my $sql = qq[select max(project_id) as id from project];
    my $sth = $dbh->prepare($sql);
    my $next_id;
    if ($sth->execute()){
	my $data = $sth->fetchrow_hashref;
	unless ($data){
            $dbh->do (qq[UNLOCK TABLES]);
            die( sprintf("Can't retrieve next project id: %s", $DBI::errstr));
	}
        $next_id = 1;
        $next_id += $data->{'id'};
    }
    else{
	die(sprintf("Can't retrieve next project id: %s", $DBI::errstr));
    }

    # OK, have next project id to use for new project
    $sql = qq[INSERT INTO project (project_id,name,hierarchy_name,changed,latest) VALUES (?,?,?,now(),true)];

    $sth = $dbh->prepare($sql);
    my $obj;
    unless ($sth->execute( $next_id, $name, $hierarchy_name)) {
        die( sprintf('DB load insert failed: %s %s', $next_id, $DBI::errstr));
    }
    $dbh->do (qq[UNLOCK TABLES]);

    return $class->new($dbh, $next_id);
}


=head2 is_name_in_database

  Arg [1]    : project name
  Arg [2]    : hierarchy name
  Example    : if(VRTrack::Project->is_name_in_database($dbh, $name,$hname)
  Description: Class method. Checks to see if a name or hierarchy name is already used in the project table.
  Returntype : boolean

=cut

sub is_name_in_database {
    my ($class, $dbh, $name, $hname) = @_;
    die "Need to call with a db handle, name, hierarchy name" unless ($dbh && $name && $hname);
    my $sql = qq[select project_id from project where latest=true and (name = ? or hierarchy_name = ?) ];
    my $sth = $dbh->prepare($sql);

    my $already_used = 0;
    if ($sth->execute($name,$hname)){
        my $data = $sth->fetchrow_hashref;
        if ($data){
            $already_used = 1;
        }
    }
    else{
        die(sprintf('Cannot retrieve project by $name: %s', $DBI::errstr));
    }
    return $already_used;
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


=head2 samples

  Arg [1]    : None
  Example    : my $samples = $project->samples();
  Description: Returns a ref to an array of the sample objects that are associated with this project
  Returntype : ref to array of VRTrack::Sample objects

=cut

sub samples {
    my ($self) = @_;

    unless ($self->{'samples'}){
        my @samples;
        foreach my $id (@{$self->sample_ids()}){
            my $obj = VRTrack::Sample->new($self->{_dbh},$id);
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
  Returntype : ref to array of integer sample IDs

=cut

sub sample_ids {
    my ($self) = @_;

    unless ($self->{'sample_ids'}){
        my $sql = qq[select sample_id from sample where project_id=? and latest=true];
        my @samples;
        my $sth = $self->{_dbh}->prepare($sql);

        if ($sth->execute($self->id)){
            foreach(@{$sth->fetchall_arrayref()}){
                push @samples, $_->[0];
            }
        }
        else{
            die(sprintf('Cannot retrieve samples: %s', $DBI::errstr));
        }

        $self->{'sample_ids'} = \@samples;
    }
 
    return $self->{'sample_ids'};
}


=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $proj->id();
               $proj->id('104');
  Description: Get/Set for ID of a project
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


=head2 hierarchy_name

  Arg [1]    : directory name (optional)
  Example    : my $hname = $project->hierarchy_name();
  Description: Get/set project hierarchy name.  This is the directory name (without path) that the project will be named in a file hierarchy.
  Returntype : string

=cut

sub hierarchy_name {
    my ($self,$name) = @_;
    if (defined $name and $name ne $self->{'hierarchy_name'}){
        $self->{'hierarchy_name'} = $name;
	$self->dirty(1);
    }
    return $self->{'hierarchy_name'};
}


=head2 name

  Arg [1]    : name (optional)
  Example    : my $name = $proj->name();
               $proj->name('1000Genomes-B1-TOS');
  Description: Get/Set for project name
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


=head2 ssid

  Arg [1]    : ssid (optional)
  Example    : my $ssid = $proj->ssid();
               $proj->ssid(104);
  Description: Get/Set for project SequenceScape ID
  Returntype : string

=cut

sub ssid {
    my ($self,$ssid) = @_;
    if (defined $ssid and $ssid ne $self->{'ssid'}){
        $self->{'ssid'} = $ssid;
	$self->dirty(1);
    }
    return $self->{'ssid'};
}


=head2 acc

  Arg [1]    : acc (optional)
  Example    : my $acc = $proj->acc();
               $proj->acc('ERA000090');
  Description: Get/Set for project accession, i.e. [SE]RR project id
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


=head2 changed

  Arg [1]    : changed (optional)
  Example    : my $changed = $project->changed();
               $project->changed('20080810123000');
  Description: Get/Set for project changed
  Returntype : string

=cut

sub changed {
    my ($self,$changed) = @_;
    if (defined $changed and $changed ne $self->{'changed'}){
	$self->{'changed'} = $changed;
	$self->dirty(1);
    }
    return $self->{'changed'};
}


=head2 add_sample

  Arg [1]    : sample name
  Example    : my $newproj = $track->add_sample('NOD mouse 1');
  Description: create a new sample, and if successful, return the object
  Returntype : VRTrack::Sample object

=cut

sub add_sample {
    my ($self, $sname) = @_;
    # Sample names should be unique for a project
    # TODO: if ssid is defined, then it should also not be added twice
    my $obj = VRTrack::Sample->new_by_name_project($self->{_dbh},$sname,$self->id);
    if ($obj){
        warn "Sample $sname is already present in the database\n";
        return undef;
    }

    $obj = VRTrack::Sample->create($self->{_dbh},$sname);
    if ($obj){
        $obj->project_id($self->id);
        $obj->update;
    }
    delete $self->{'sample_ids'};
    delete $self->{'samples'};
    return $obj;
}


=head2 get_sample_by_name

  Arg [1]    : sample name
  Example    : my $sample = $track->get_sample_by_name('My sample');
  Description: retrieve sample object by name
  Returntype : VRTrack::Sample object

=cut

sub get_sample_by_name {
    my ($self, $name) = @_;
    #my $obj = VRTrack::Sample->new_by_name_project($self->{_dbh},$name, $self->id);
    my @match = grep {$_->name eq $name} @{$self->samples};
    if (scalar @match > 1){ # shouldn't happen
        die "More than one sample with name $name";
    }
    my $obj;
    if (@match){
        $obj = $match[0];
    }

    return $obj;
}


=head2 get_sample_by_id

  Arg [1]    : sample id 
  Example    : my $sample = $proj->get_sample_by_id(1154);
  Description: retrieve sample object by internal id
  Returntype : VRTrack::Sample object

=cut

sub get_sample_by_id {
    my ($self, $id) = @_;
    my @match = grep {$_->id == $id} @{$self->samples};
    if (scalar @match > 1){ # shouldn't happen
        die "More than one sample with id $id";
    }
    my $obj;
    if (@match){
        $obj = $match[0];
    }

    return $obj;
}


=head2 get_sample_by_ssid

  Arg [1]    : sample sequencescape id
  Example    : my $sample = $proj->get_sample_by_ssid(1154);
  Description: retrieve sample object by sequencescape id
  Returntype : VRTrack::Sample object

=cut

sub get_sample_by_ssid {
    my ($self, $ssid) = @_;
    my @match = grep {$_->ssid == $ssid} @{$self->samples};
    if (scalar @match > 1){ # shouldn't happen
        die "More than one sample with ssid $ssid";
    }
    my $obj;
    if (@match){
        $obj = $match[0];
    }

    return $obj;
}


=head2 descendants

  Arg [1]    : none
  Example    : my $desc_objs = $obj->descendants();
  Description: Returns a ref to an array of all objects that are descendants of this object
  Returntype : arrayref of objects

=cut

sub descendants {
    my ($self) = @_;
    my @desc;
    foreach (@{$self->samples}){
        push @desc, $_;
        push @desc, @{$_->descendants};
    }
    return \@desc;
}

1;
