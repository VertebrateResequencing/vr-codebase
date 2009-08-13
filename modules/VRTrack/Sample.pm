package VRTrack::Sample;
=head1 NAME

VRTrack::Sample - Sequence Tracking Sample object

=head1 SYNOPSIS
    my $samp = VRTrack::Sample->new($dbh, $sample_id);

    #get arrayref of library objects in a sample
    my $libs = $sample->libraries();
    
    my $id = $sample->id();
    my $name = $sample->name();

=head1 DESCRIPTION

An object describing the tracked properties of a sample.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;
no warnings 'uninitialized';
use VRTrack::Library;
use VRTrack::Individual;
use constant DBI_DUPLICATE => '1062';

###############################################################################
# Class methods
###############################################################################

=head2 new

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : sample id
  Example    : my $samp = VRTrack::Sample->new($dbh, $id)
  Description: Class method. Returns Sample object by sample_id
  Returntype : VRTrack::Sample object

=cut

sub new {
    my ($class,$dbh, $id) = @_;
    die "Need to call with a db handle and id" unless ($dbh && $id);
    my $self = {};
    bless ($self, $class);
    $self->{_dbh} = $dbh;

    my $sql = qq[select sample_id, project_id, ssid, name, acc, individual_id, changed, latest from sample where sample_id = ? and latest=true];
    my $sth = $self->{_dbh}->prepare($sql);

    if ($sth->execute($id)){
        my $data = $sth->fetchrow_hashref;
        unless ($data){
            return undef;
        }
        $self->id($data->{'sample_id'});
        $self->project_id($data->{'project_id'});
        $self->ssid($data->{'ssid'});
        $self->name($data->{'name'});
        $self->acc($data->{'acc'});
        $self->individual_id($data->{'individual_id'});
        $self->changed($data->{'changed'});
	$self->dirty(0); # unset the dirty flag
    }
    else{
        die(sprintf('Cannot retrieve sample: %s', $DBI::errstr));
    }

    return $self;
}


=head2 new_by_name_project

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : sample name
  Arg [3]    : project id
  Example    : my $sample = VRTrack::Sample->new_by_name_project($dbh, $name,$project_id)
  Description: Class method. Returns latest Sample object by name and project_id.  If no such name is in the database, returns undef
  Returntype : VRTrack::Sample object

=cut

sub new_by_name_project {
    my ($class,$dbh, $name, $project_id) = @_;
    die "Need to call with a db handle, name, project_id" unless ($dbh && $name && $project_id);
    my $sql = qq[select sample_id from sample where name = ? and project_id = ? and latest = true];
    my $sth = $dbh->prepare($sql);

    my $id;
    if ($sth->execute($name, $project_id)){
        my $data = $sth->fetchrow_hashref;
        unless ($data){
            return undef;
        }
        $id = $data->{'sample_id'};
    }
    else{
        die(sprintf('Cannot retrieve sample by $name, $project: %s', $DBI::errstr));
    }
    return $class->new($dbh, $id);
}


=head2 create

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : sample name
  Example    : my $sample = VRTrack::Sample->create($dbh, $name)
  Description: Class method.  Creates new Sample object in the database.
  Returntype : VRTrack::Sample object

=cut

sub create {
    my ($class,$dbh, $name) = @_;
    die "Need to call with a db handle and name" unless ($dbh && $name);
    $dbh->do (qq[LOCK TABLE sample WRITE]);
    my $sql = qq[select max(sample_id) as id from sample];
    my $sth = $dbh->prepare($sql);
    my $next_id;
    if ($sth->execute()){
	my $data = $sth->fetchrow_hashref;
	unless ($data){
            $dbh->do (qq[UNLOCK TABLES]);
            die( sprintf("Can't retrieve next sample id: %s", $DBI::errstr));
	}
        $next_id = $data->{'id'};
        $next_id++;
    }
    else{
	die(sprintf("Can't retrieve next sample id: %s", $DBI::errstr));
    }

    $sql = qq[INSERT INTO sample (sample_id, name, changed, latest) 
                 VALUES (?,?,now(),true)];

    $sth = $dbh->prepare($sql);
    unless ($sth->execute( $next_id, $name )) {
        $dbh->do (qq[UNLOCK TABLES]);
        die( sprintf('DB load insert failed: %s %s', $next_id, $DBI::errstr));
    }

    $dbh->do (qq[UNLOCK TABLES]);

    return $class->new($dbh, $next_id);
}


###############################################################################
# Object methods
###############################################################################

=head2 libraries

  Arg [1]    : None
  Example    : my $libraries = $sample->libraries();
  Description: Returns a ref to an array of the sample objects that are associated with this sample.
  Returntype : ref to array of VRTrack::Sample objects

=cut

sub libraries {
    my ($self) = @_;

    unless ($self->{'libraries'}){
        my @libraries;
        foreach my $id (@{$self->library_ids()}){
            my $obj = VRTrack::Library->new($self->{_dbh},$id);
            push @libraries, $obj;
        }
        $self->{'libraries'} = \@libraries;
    }

    return $self->{'libraries'};
}


=head2 library_ids

  Arg [1]    : None
  Example    : my $library_ids = $sample->library_ids();
  Description: Returns a ref to an array of the library IDs that are associated with this sample
  Returntype : ref to array of integer library IDs

=cut

sub library_ids {
    my ($self) = @_;

    unless ($self->{'library_ids'}){
        my $sql = qq[select distinct(library_id) from library where sample_id=? and latest=true];
        my @libraries;
        my $sth = $self->{_dbh}->prepare($sql);

        if ($sth->execute($self->id)){
            foreach(@{$sth->fetchall_arrayref()}){
                push @libraries, $_->[0];
            }
        }
        else{
            die(sprintf('Cannot retrieve libraries: %s', $DBI::errstr));
        }

        $self->{'library_ids'} = \@libraries;
    }
 
    return $self->{'library_ids'};
}


=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $samp->id();
               $samp->id('104');
  Description: Get/Set for ID of a sample
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


=head2 project_id

  Arg [1]    : project_id (optional)
  Example    : my $project_id = $samp->project_id();
               $samp->project_id('104');
  Description: Get/Set for ID of a sample
  Returntype : SequenceScape ID (usu. integer)

=cut

sub project_id {
    my ($self,$project_id) = @_;
    if (defined $project_id and $project_id ne $self->{'project_id'}){
        $self->{'project_id'} = $project_id;
	$self->dirty(1);
    }
    return $self->{'project_id'};
}


=head2 hierarchy_name

  Arg [1]    : directory name (optional)
  Example    : my $hname = $sample->hierarchy_name();
  Description: Get sample hierarchy name.  This is the directory name (without path) that the sample will be named in a file hierarchy.  Note that this is actually the individual hierarchy name, and is only gettable here for convenience.  Setting is done via the individual object.
  Returntype : string

=cut

sub hierarchy_name {
    my ($self) = @_;
    my $name = undef;
    if ($self->individual){
        $name = $self->individual->hierarchy_name();
    }
    return $name;
}


=head2 name

  Arg [1]    : name (optional)
  Example    : my $name = $samp->name();
               $samp->name('104');
  Description: Get/Set for sample name
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
  Example    : my $ssid = $samp->ssid();
               $samp->ssid(104);
  Description: Get/Set for sample SequenceScape ID
  Returntype : SequenceScape ID integer

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
  Example    : my $acc = $samp->acc();
               $samp->acc('ERS000090');
  Description: Get/Set for sample accession, i.e. [SE]RR sample id
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


=head2 individual_id

  Arg [1]    : individual_id (optional)
  Example    : my $individual_id = $samp->individual_id();
               $samp->individual_id(123);
  Description: Get/Set for sample internal individual_id
  Returntype : integer

=cut

sub individual_id {
    my ($self,$individual_id) = @_;
    if (defined $individual_id and $individual_id ne $self->{'individual_id'}){
        $self->{'individual_id'} = $individual_id;
	$self->dirty(1);
    }
    return $self->{'individual_id'};
}


=head2 individual

  Arg [1]    : individual name (optional)
  Example    : my $individual = $samp->individual();
               $samp->individual('NA19820');
  Description: Get/Set for sample individual.  Lazy-loads individual object from $self->individual_id.  If a individual name is supplied, then individual_id is set to the corresponding individual in the database.  If no such individual exists, returns undef.  Use add_individual to add a individual in this case.
  Returntype : VRTrack::Individual object

=cut

sub individual {
    my ($self,$individual) = @_;
    if ($individual){
        # get existing individual by name
        my $obj = $self->get_individual_by_name($individual);
        if ($obj){
            $self->{'individual'} = $obj;
            $self->{'individual_id'} = $obj->id;
            $self->dirty(1);
        }
        else {
            # warn "No such individual in the database";
            return undef; # explicitly return nothing.
        }
    }
    elsif ($self->{'individual'}){
        # already got a individual object.  We'll return it at the end.
    }
    else {  # lazy-load individual from database
        if ($self->individual_id){
            my $obj = VRTrack::Individual->new($self->{_dbh},$self->individual_id);
            $self->{'individual'} = $obj;
        }
    }
    return $self->{'individual'};
}


=head2 add_individual

  Arg [1]    : individual name
  Example    : my $ind = $samp->add_individual('NA19820');
  Description: create a new individual, and if successful, return the object
  Returntype : VRTrack::Library object

=cut

sub add_individual {
    my ($self, $name) = @_;

    my $obj = $self->get_individual_by_name($name);
    if ($obj){
        warn "Individual $name is already present in the database\n";
        return undef;
    }
    else {
        my $pop = VRTrack::Individual->create($self->{_dbh}, $name);
        # populate caches
        $self->{'individual_id'} = $pop->id;
        $self->{'individual'} = $pop;
        $self->dirty(1);
    }
    return $self->{'individual'};
}


=head2 get_individual_by_name

  Arg [1]    : individual_name
  Example    : my $ind = $samp->get_individual_by_name('NA19820');
  Description: Retrieve a VRTrack::Individual object by name
  Returntype : VRTrack::Individual object

=cut

sub get_individual_by_name {
    my ($self,$name) = @_;
    return VRTrack::Individual->new_by_name($self->{_dbh}, $name);
}


=head2 add_library

  Arg [1]    : library name
  Example    : my $newlib = $samp->add_library('NOD_500_SLX_1');
  Description: create a new library, and if successful, return the object
  Returntype : VRTrack::Library object

=cut

sub add_library {
    my ($self, $name) = @_;

    # Check for unwanted duplicates.  Should not allow two libraries to
    # have the same name even though this can happen in sequencescape.

    my $obj = $self->get_library_by_name($name);
    if ($obj){
        warn "Library $name is already present in the database\n";
        return undef;
    }
    $obj = VRTrack::Library->create($self->{_dbh}, $name);
    if ($obj){
        $obj->sample_id($self->id);
        $obj->update;
    }
    # clear caches
    delete $self->{'library_ids'};
    delete $self->{'libraries'};

    return $obj;
}


=head2 get_library_by_id

  Arg [1]    : library id from sequencescape
  Example    : my $library = $sam->get_library_by_id(1930);
  Description: retrieve library object by sequencescape id
  Returntype : VRTrack::Library object

=cut

sub get_library_by_id {
    my ($self, $id) = @_;
    my $obj = VRTrack::Library->new($self->{_dbh},$id);
    return $obj;
}


=head2 get_library_by_name

  Arg [1]    : library name
  Example    : my $library = $track->get_library_by_name('My library');
  Description: retrieve library object by name
  Returntype : VRTrack::Library object

=cut

sub get_library_by_name {
    my ($self, $name) = @_;
    # my $obj = VRTrack::Library->new_by_name($self->{_dbh},$name);
    # unless ($obj->sample_id == $self->sample_id){
    #    die "Library $name does not belong to sample ",$self->name,"\n";
    #}
    my @match = grep {$_->name eq $name} @{$self->libraries};
    if (scalar @match > 1){ # shouldn't happen
        die "More than one matching library with name $name";
    }
    my $obj;
    if (@match){
        $obj = $match[0];
    }
    return $obj;
}


=head2 changed

  Arg [1]    : changed (optional)
  Example    : my $changed = $sample->changed();
               $sample->changed('20080810123000');
  Description: Get/Set for sample changed
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


=head2 update

  Arg [1]    : None
  Example    : $sample->update();
  Description: Update a sample whose properties you have changed.  If properties haven't changed (i.e. dirty flag is unset) do nothing.  
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
	my $save_ac = $dbh->{AutoCommit};
	$dbh->{RaiseError} = 1; # raise exception if an error occurs
	$dbh->{PrintError} = 0; # don't print an error message
	$dbh->{AutoCommit} = 0; # disable auto-commit

	eval {
	    # Need to unset 'latest' flag on current latest file and add
	    # the new file details with the latest flag set
	    my $updsql = qq[UPDATE sample SET latest=false WHERE sample_id = ? and latest=true];
	    
	    my $addsql = qq[INSERT INTO sample (sample_id, project_id, ssid, name, acc, individual_id, changed, latest) 
			    VALUES (?,?,?,?,?,?,now(),true)];
	    $dbh->do ($updsql, undef,$self->id);
	    $dbh->do ($addsql, undef,$self->id, $self->project_id, $self->ssid, $self->name, $self->acc, $self->individual_id);
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
