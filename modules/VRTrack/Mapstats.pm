package VRTrack::Mapstats; 
=head1 NAME

VRTrack::Mapstats - Sequence Tracking Mapstats object

=head1 SYNOPSIS
    my $mapstats= VRTrack::Mapstats->new($dbh, $mapstats_id);

    my $id = $mapstats->id();
    my $qc_status = $mapstats->qc_status();

=head1 DESCRIPTION

An object describing a set of statistics for a specific lane mapping.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;
no warnings 'uninitialized';
use constant DBI_DUPLICATE => '1062';
use VRTrack::Image;
use VRTrack::Mapper;
use VRTrack::Assembly;

###############################################################################
# Class methods
###############################################################################

=head2 new

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : mapstats name
  Example    : my $mapstats= VRTrack::Mapstats->new($dbh, $name)
  Description: Returns Mapstats object by mapstats name
  Returntype : VRTrack::Mapstats object

=cut

sub new {
    my ($class,$dbh, $id) = @_;
    die "Need to call with a db handle and id" unless ($dbh && $id);
    my $self = {};
    bless ($self, $class);
    $self->{_dbh} = $dbh;

    my $sql = qq[select mapstats_id, lane_id, mapper_id, assembly_id, reads_mapped, reads_paired, bases_mapped, rmdup_reads_mapped, rmdup_bases_mapped, error_rate, mean_insert, sd_insert, changed, latest from mapstats where mapstats_id = ? and latest = true];
    my $sth = $self->{_dbh}->prepare($sql);

    if ($sth->execute($id)){
        my $data = $sth->fetchrow_hashref;
	unless ($data){
	    return undef;
	}
	$self->id($data->{'mapstats_id'});
	$self->lane_id($data->{'lane_id'});
	$self->mapper_id($data->{'mapper_id'});
	$self->assembly_id($data->{'assembly_id'});
	$self->reads_mapped($data->{'reads_mapped'});
	$self->reads_paired($data->{'reads_paired'});
	$self->bases_mapped($data->{'bases_mapped'});
	$self->rmdup_reads_mapped($data->{'rmdup_reads_mapped'});
	$self->rmdup_bases_mapped($data->{'rmdup_bases_mapped'});
	$self->error_rate($data->{'error_rate'});
	$self->mean_insert($data->{'mean_insert'});
	$self->sd_insert($data->{'sd_insert'});
        $self->changed($data->{'changed'});
	$self->dirty(0);    # unset the dirty flag
    }
    else{
	die(sprintf('Cannot retrieve mapstats: %s', $DBI::errstr));
    }

    return $self;
}


=head2 create

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : lane id that this mapping relates to
  Example    : my $mapstats = VRTrack::Mapstats->create($dbh, $lane_id)
  Description: Class method.  Creates new Mapstats object in the database.
  Returntype : VRTrack::Mapstats object

=cut

sub create {
    my ($class,$dbh, $lane_id) = @_;
    die "Need to call with a db handle and lane_id" unless ($dbh && $lane_id);
    $dbh->do (qq[LOCK TABLE mapstats WRITE]);
    my $sql = qq[select max(mapstats_id) as id from mapstats];
    my $sth = $dbh->prepare($sql);
    my $next_id;
    if ($sth->execute()){
	my $data = $sth->fetchrow_hashref;
	unless ($data){
            $dbh->do (qq[UNLOCK TABLES]);
            die( sprintf("Can't retrieve next mapstats id: %s", $DBI::errstr));
	}
        $next_id = $data->{'id'};
        $next_id++;
    }
    else{
	die(sprintf("Can't retrieve next mapstats id: %s", $DBI::errstr));
    }

    $sql = qq[INSERT INTO mapstats (mapstats_id, lane_id, changed, latest) 
                 VALUES (?,?,now(),true)];

    $sth = $dbh->prepare($sql);
    unless ($sth->execute( $next_id, $lane_id)) {
        $dbh->do (qq[UNLOCK TABLES]);
        die( sprintf('DB load insert failed: %s %s', $next_id, $DBI::errstr));
    }

    $dbh->do (qq[UNLOCK TABLES]);

    return $class->new($dbh, $next_id);
}


###############################################################################
# Object methods
###############################################################################

=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $mapstats->id();
	       $mapstats->id(1);
  Description: Get/Set for internal id of a mapping
  Returntype : integer

=cut

sub id {
    my ($self,$id) = @_;
    if (defined $id and $id != $self->{'id'}){
	$self->{'id'} = $id;
	$self->dirty(1);
    }
    return $self->{'id'};
}


=head2 lane_id

  Arg [1]    : lane_id (optional)
  Example    : my $lane_id = $mapstats->lane_id();
	       $mapstats->lane_id(104);
  Description: Get/Set for internal ID of the lane which was mapped
  Returntype : integer

=cut

sub lane_id {
    my ($self,$lane_id) = @_;
    if (defined $lane_id and $lane_id != $self->{'lane_id'}){
	$self->{'lane_id'} = $lane_id;
	$self->dirty(1);
    }
    return $self->{'lane_id'};
}


=head2 mapper_id

  Arg [1]    : mapper_id (optional)
  Example    : my $mapper_id = $mapstats->mapper_id();
               $mapstats->mapper_id(104);
  Description: Get/Set for internal ID of the mapper used for the mapping
  Returntype : integer

=cut

sub mapper_id {
    my ($self,$mapper_id) = @_;
    if (defined $mapper_id and $mapper_id != $self->{'mapper_id'}){
        $self->{'mapper_id'} = $mapper_id;
        $self->dirty(1);
    }
    return $self->{'mapper_id'};
}


=head2 mapper

  Arg [1]    : mapper name (optional)
  Arg [2]    : mapper version (optional)
  Example    : my $mapper = $mapstats->mapper();
               $mapstats->mapper('maq','0.7.1-6');
  Description: Get/Set for mapping mapper.  Lazy-loads mapper object from $self->mapper_id.  If a mapper name and version is supplied, then mapper_id is set to the corresponding mapper in the database.  If no such mapper exists, returns undef.  Use add_mapper to add a mapper in this case.
  Returntype : VRTrack::Mapper object

=cut

sub mapper {
    my ($self,$mapper, $version) = @_;
    if ($mapper && $version){
        # get existing mapper by name,version
        my $obj = $self->get_mapper_by_name_version($mapper,$version);
        if ($obj){
            $self->{'mapper'} = $obj;
            $self->{'mapper_id'} = $obj->id;
            $self->dirty(1);
        }
        else {
            # warn "No such mapper in the database";
            return undef; # explicitly return nothing.
        }
    }
    elsif ($self->{'mapper'}){
        # already got a mapper object.  We'll return it at the end.
    }
    else {  # lazy-load mapper from database
        if ($self->mapper_id){
            my $obj = VRTrack::Mapper->new($self->{_dbh},$self->mapper_id);
            $self->{'mapper'} = $obj;
        }
    }
    return $self->{'mapper'};
}


=head2 add_mapper

  Arg [1]    : mapper name
  Arg [2]    : mapper version
  Example    : my $map = $mapstats->add_mapper('maq','0.7.1-6');
  Description: create a new mapper, and if successful, return the object
  Returntype : VRTrack::Library object

=cut

sub add_mapper {
    my ($self, $name,$version) = @_;

    my $obj = $self->get_mapper_by_name_version($name,$version);
    if ($obj){
        warn "Mapper $name $version is already present in the database\n";
        return undef;
    }
    else {
        $obj = VRTrack::Mapper->create($self->{_dbh}, $name, $version);
        # populate caches
        $self->{'mapper_id'} = $obj->id;
        $self->{'mapper'} = $obj;
        $self->dirty(1);
    }
    return $self->{'mapper'};
}


=head2 get_mapper_by_name_version

  Arg [1]    : mapper_name
  Arg [2]    : mapper version
  Example    : my $map = $mapstats->get_mapper_by_name_version('maq','0.7.1-6');
  Description: Retrieve a VRTrack::Mapper object by name
  Returntype : VRTrack::Mapper object

=cut

sub get_mapper_by_name_version {
    my ($self,$name,$version) = @_;
    return VRTrack::Mapper->new_by_name_version($self->{_dbh}, $name,$version);
}


=head2 assembly_id

  Arg [1]    : assembly_id (optional)
  Example    : my $assembly_id = $mapstats->assembly_id();
               $mapstats->assembly_id(104);
  Description: Get/Set for internal ID of the assembly used for the mapping
  Returntype : integer

=cut

sub assembly_id {
    my ($self,$assembly_id) = @_;
    if (defined $assembly_id and $assembly_id != $self->{'assembly_id'}){
        $self->{'assembly_id'} = $assembly_id;
        $self->dirty(1);
    }
    return $self->{'assembly_id'};
}


=head2 assembly

  Arg [1]    : assembly name (optional)
  Example    : my $assembly = $mapstats->assembly();
               $mapstats->assembly('NCBIm36');
  Description: Get/Set for mapping assembly.  Lazy-loads assembly object from $self->assembly_id.  If a assembly name is supplied, then assembly_id is set to the corresponding assembly in the database.  If no such assembly exists, returns undef.  Use add_assembly to add a assembly in this case.
  Returntype : VRTrack::Assembly object

=cut

sub assembly {
    my ($self,$assembly) = @_;
    if ($assembly){
        # get existing assembly by name
        my $obj = $self->get_assembly_by_name($assembly);
        if ($obj){
            $self->{'assembly'} = $obj;
            $self->{'assembly_id'} = $obj->id;
            $self->dirty(1);
        }
        else {
            # warn "No such assembly in the database";
            return undef; # explicitly return nothing.
        }
    }
    elsif ($self->{'assembly'}){
        # already got a assembly object.  We'll return it at the end.
    }
    else {  # lazy-load assembly from database
        if ($self->assembly_id){
            my $obj = VRTrack::Assembly->new($self->{_dbh},$self->assembly_id);
            $self->{'assembly'} = $obj;
        }
    }
    return $self->{'assembly'};
}


=head2 add_assembly

  Arg [1]    : assembly name
  Example    : my $assy = $mapping->add_assembly('NCBIm36');
  Description: create a new assembly, and if successful, return the object
  Returntype : VRTrack::Library object

=cut

sub add_assembly {
    my ($self, $name) = @_;

    my $obj = $self->get_assembly_by_name($name);
    if ($obj){
        warn "Assembly $name is already present in the database\n";
        return undef;
    }
    else {
        $obj = VRTrack::Assembly->create($self->{_dbh}, $name);
        # populate caches
        $self->{'assembly_id'} = $obj->id;
        $self->{'assembly'} = $obj;
        $self->dirty(1);
    }
    return $self->{'assembly'};
}


=head2 get_assembly_by_name

  Arg [1]    : assembly_name
  Example    : my $map = $mapstats->get_assembly_by_name_version('maq');
  Description: Retrieve a VRTrack::Assembly object by name
  Returntype : VRTrack::Assembly object

=cut

sub get_assembly_by_name {
    my ($self,$name,$version) = @_;
    return VRTrack::Assembly->new_by_name($self->{_dbh}, $name);
}


=head2 reads_mapped

  Arg [1]    : number of mapped reads in mapping (optional)
  Example    : my $num_reads = $mapstats->reads_mapped();
	       $mapstats->reads_mapped(1_000_000);
  Description: Get/Set for total number of mapped reads in mapping
  Returntype : integer

=cut

sub reads_mapped {
    my ($self,$num_reads) = @_;
    if (defined $num_reads and $num_reads != $self->{'reads_mapped'}){
	$self->{'reads_mapped'} = $num_reads;
	$self->dirty(1);
    }
    return $self->{'reads_mapped'};
}


=head2 reads_paired

  Arg [1]    : number of paired mapped reads in mapping (optional)
  Example    : my $num_reads = $mapstats->reads_paired();
	       $mapstats->reads_paired(1_000_000);
  Description: Get/Set for number of paired mapped reads in mapping
  Returntype : integer

=cut

sub reads_paired {
    my ($self,$num_reads) = @_;
    if (defined $num_reads and $num_reads != $self->{'reads_paired'}){
	$self->{'reads_paired'} = $num_reads;
	$self->dirty(1);
    }
    return $self->{'reads_paired'};
}


=head2 bases_mapped

  Arg [1]    : number of mapped bases in mapping (optional)
  Example    : my $num_bases = $mapstats->bases_mapped();
	       $mapstats->bases_mapped(1_000_000);
  Description: Get/Set for total number of mapped bases in mapping
  Returntype : integer

=cut

sub bases_mapped {
    my ($self,$num_bases) = @_;
    if (defined $num_bases and $num_bases != $self->{'bases_mapped'}){
	$self->{'bases_mapped'} = $num_bases;
	$self->dirty(1);
    }
    return $self->{'bases_mapped'};
}


=head2 rmdup_reads_mapped

  Arg [1]    : number of mapped reads in mapping (optional)
  Example    : my $num_reads = $mapstats->rmdup_reads_mapped();
	       $mapstats->rmdup_reads_mapped(1_000_000);
  Description: Get/Set for number of mapped reads in mapping after removing duplicates
  Returntype : integer

=cut

sub rmdup_reads_mapped {
    my ($self,$num_reads) = @_;
    if (defined $num_reads and $num_reads != $self->{'rmdup_reads_mapped'}){
	$self->{'rmdup_reads_mapped'} = $num_reads;
	$self->dirty(1);
    }
    return $self->{'rmdup_reads_mapped'};
}


=head2 rmdup_bases_mapped

  Arg [1]    : number of mapped bases in mapping (optional)
  Example    : my $num_bases = $mapstats->rmdup_bases_mapped();
	       $mapstats->rmdup_bases_mapped(1_000_000);
  Description: Get/Set for number of mapped bases in mapping after removing duplicate reads
  Returntype : integer

=cut

sub rmdup_bases_mapped {
    my ($self,$num_bases) = @_;
    if (defined $num_bases and $num_bases != $self->{'rmdup_bases_mapped'}){
	$self->{'rmdup_bases_mapped'} = $num_bases;
	$self->dirty(1);
    }
    return $self->{'rmdup_bases_mapped'};
}


=head2 error_rate

  Arg [1]    : maq/mapping error rate (optional)
  Example    : my $err = $mapstats->error_rate();
	       $mapstats->error_rate(0.01);
  Description: Get/Set for mapper error rate in mapping
  Returntype : float

=cut

sub error_rate {
    my ($self,$err) = @_;
    if (defined $err and $err != $self->{'error_rate'}){
	$self->{'error_rate'} = $err;
	$self->dirty(1);
    }
    return $self->{'error_rate'};
}


=head2 mean_insert

  Arg [1]    : maq/mapping insert size (optional)
  Example    : my $ins = $mapstats->mean_insert();
	       $mapstats->mean_insert(243.4);
  Description: Get/Set for mean mapped insert size for pairs
  Returntype : float

=cut

sub mean_insert {
    my ($self,$ins) = @_;
    if (defined $ins and $ins != $self->{'mean_insert'}){
	$self->{'mean_insert'} = $ins;
	$self->dirty(1);
    }
    return $self->{'mean_insert'};
}


=head2 sd_insert

  Arg [1]    : maq/mapping insert sd (optional)
  Example    : my $ins = $mapstats->sd_insert();
	       $mapstats->sd_insert(15.64);
  Description: Get/Set for Standard Deviation of mapped insert size for pairs
  Returntype : float

=cut

sub sd_insert {
    my ($self,$ins) = @_;
    if (defined $ins and $ins != $self->{'sd_insert'}){
	$self->{'sd_insert'} = $ins;
	$self->dirty(1);
    }
    return $self->{'sd_insert'};
}


=head2 dirty

  Arg [1]    : boolean for dirty status
  Example    : $mapstats->dirty(1);
  Description: Get/Set for mapstats properties having been altered.
  Returntype : boolean

=cut

sub dirty {
    my ($self,$dirty) = @_;
    if (defined $dirty){
	$self->{_dirty} = $dirty ? 1 : 0;
    }
    return $self->{_dirty};
}


=head2 add_image

  Arg [1]    : image name
  Arg [2]    : image data (i.e. the binary image)
  Example    : my $newimage = $lib->add_image('gc.png',$gc_png_img);
  Description: create a new image, and if successful, return the object
  Returntype : VRTrack::Image object

=cut

sub add_image {
    my ($self, $name, $img) = @_;
    my $obj = VRTrack::Image->create($self->{_dbh}, $name, $img);
    if ($obj){
        $obj->mapstats_id($self->id);
        $obj->update;
    }
    # clear caches
    delete $self->{'image_ids'};
    delete $self->{'images'};

    return $obj;
}


=head2 images

  Arg [1]    : None
  Example    : my $images = $mapping->images();
  Description: Returns a ref to an array of the image objects that are associated with this mapping.
  Returntype : ref to array of VRTrack::Image objects

=cut

sub images {
    my ($self) = @_;
    unless ($self->{'images'}){
	my @images;
    	foreach my $id (@{$self->image_ids()}){
	    my $obj = VRTrack::Image->new($self->{_dbh},$id);
	    push @images, $obj;
	}
	$self->{'images'} = \@images;
    }

    return $self->{'images'};
}


=head2 image_ids

  Arg [1]    : None
  Example    : my $image_ids = $mapping->image_ids();
  Description: Returns a ref to an array of the image ids that are associated with this mapping.
  Returntype : ref to array of image ids

=cut

sub image_ids {
    my ($self) = @_;
    unless ($self->{'image_ids'}){
	my $sql = qq[select image_id from image where mapstats_id=?];
	my @images;
	my $sth = $self->{_dbh}->prepare($sql);

	if ($sth->execute($self->id)){
	    foreach(@{$sth->fetchall_arrayref()}){
		push @images, $_->[0];
	    }
	}
	else{
	    die(sprintf('Cannot retrieve images: %s', $DBI::errstr));
	}

	$self->{'image_ids'} = \@images;
    }
 
    return $self->{'image_ids'};
}


=head2 changed

  Arg [1]    : changed (optional)
  Example    : my $changed = $mapstats->changed();
	       $mapstats->changed('20090101102300');
  Description: Get/Set for mapstats changed
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


=head2 update

  Arg [1]    : None
  Example    : $mapstats->update();
  Description: Update a mapstats whose properties you have changed.  If properties haven't changed (i.e. dirty flag is unset) do nothing.  
	       Changes the changed datestamp to now() on the mysql server (i.e. you don't have to set changed yourself, and indeed if you do, it will be overridden).
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
	    # Need to unset 'latest' flag on current latest mapstats and add
	    # the new mapstats details with the latest flag set
	    my $updsql = qq[UPDATE mapstats SET latest=false WHERE mapstats_id = ? and latest=true];
	    
	    my $addsql = qq[INSERT INTO mapstats (mapstats_id, lane_id, mapper_id, assembly_id, reads_mapped, reads_paired, bases_mapped, rmdup_reads_mapped, rmdup_bases_mapped, error_rate, mean_insert, sd_insert, changed, latest) 
			    VALUES (?,?,?,?,?,?,?,?,?,?,?,?,now(),true)];
	    $dbh->do ($updsql, undef,$self->id);
	    $dbh->do ($addsql, undef,$self->id, $self->lane_id, $self->mapper_id, $self->assembly_id, $self->reads_mapped, $self->reads_paired, $self->bases_mapped, $self->rmdup_reads_mapped, $self->rmdup_bases_mapped, $self->error_rate, $self->mean_insert, $self->sd_insert);
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

    return $success;
}

1;
