package VRTrack::Mapstats; 
# author: jws
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

use VRTrack::Core_obj;
use VRTrack::Image;
use VRTrack::Mapper;
use VRTrack::Assembly;
use File::Basename;
our @ISA = qw(VRTrack::Core_obj);

=head2 fields_dispatch

  Arg [1]    : none
  Example    : my $fieldsref = $mapstats->fields_dispatch();
  Description: Returns hashref dispatch table keyed on database field
               Used internally for new and update methods
  Returntype : hashref

=cut

sub fields_dispatch {
    my $self = shift;
    my %fields = ( 
                'mapstats_id'        => sub { $self->id(@_)},
                'lane_id'            => sub { $self->lane_id(@_)},
                'mapper_id'          => sub { $self->mapper_id(@_)},
                'assembly_id'        => sub { $self->assembly_id(@_)},
                'raw_reads'          => sub { $self->raw_reads(@_)},
                'raw_bases'          => sub { $self->raw_bases(@_)},
                'clip_bases'         => sub { $self->clip_bases(@_)},
                'reads_mapped'       => sub { $self->reads_mapped(@_)},
                'reads_paired'       => sub { $self->reads_paired(@_)},
                'bases_mapped'       => sub { $self->bases_mapped(@_)},
                'rmdup_reads_mapped' => sub { $self->rmdup_reads_mapped(@_)},
                'rmdup_bases_mapped' => sub { $self->rmdup_bases_mapped(@_)},
                'error_rate'         => sub { $self->error_rate(@_)},
                'mean_insert'        => sub { $self->mean_insert(@_)},
                'sd_insert'          => sub { $self->sd_insert(@_)},
                'adapter_reads'      => sub { $self->adapter_reads(@_)},
                'gt_expected'        => sub { $self->genotype_expected(@_)},
                'gt_found'           => sub { $self->genotype_found(@_)},
                'gt_ratio'           => sub { $self->genotype_ratio(@_)},
                'note_id'            => sub { $self->note_id(@_)},
                'changed'            => sub { $self->changed(@_)},
                'latest'             => sub { $self->is_latest(@_)},
                );

    return \%fields;
}

###############################################################################
# Class methods
###############################################################################

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
        delete $self->{'mapper'};
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
            if ($self->mapper_id != $obj->id){
                $self->mapper_id($obj->id);
                $self->dirty(1);
            }
            $self->{'mapper'} = $obj;
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
        delete $self->{'assembly'};
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
            # Have we actually changed?
            if ($self->assembly_id != $obj->id){
                $self->assembly_id($obj->id);
                $self->dirty(1);
            }
            $self->{'assembly'} = $obj;
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


=head2 raw_reads

  Arg [1]    : number of raw reads used in mapping (optional)
  Example    : my $num_reads = $mapstats->raw_reads();
	       $mapstats->raw_reads(1_000_000);
  Description: Get/Set for total number of raw reads used in mapping
  Returntype : integer

=cut

sub raw_reads {
    my ($self,$num_reads) = @_;
    if (defined $num_reads and $num_reads != $self->{'raw_reads'}){
	$self->{'raw_reads'} = $num_reads;
	$self->dirty(1);
    }
    return $self->{'raw_reads'};
}


=head2 raw_bases

  Arg [1]    : number of raw bases used in mapping (optional)
  Example    : my $num_bases = $mapstats->raw_bases();
	       $mapstats->raw_bases(1_000_000);
  Description: Get/Set for total number of raw bases used in mapping
  Returntype : integer

=cut

sub raw_bases {
    my ($self,$num_bases) = @_;
    if (defined $num_bases and $num_bases != $self->{'raw_bases'}){
	$self->{'raw_bases'} = $num_bases;
	$self->dirty(1);
    }
    return $self->{'raw_bases'};
}


=head2 clip_bases

  Arg [1]    : number of raw bases after clipping
  Example    : my $bases_after_clipping = $mapstats->clip_bases();
	       $mapstats->clip_bases(1_000_000);
  Description: Get/Set for total number of raw bases after clipping
  Returntype : integer

=cut

sub clip_bases {
    my ($self,$num_bases) = @_;
    if (defined $num_bases and $num_bases != $self->{'clip_bases'}){
	$self->{'clip_bases'} = $num_bases;
	$self->dirty(1);
    }
    return $self->{'clip_bases'};
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


=head2 adapter_reads

  Arg [1]    : number of reads containing adapter sequence
  Example    : my $num_adap_reads = $mapstats->adapter_reads();
	       $mapstats->adapter_reads(1_000_000);
  Description: Get/Set for total number of reads containing adapter sequence
  Returntype : integer

=cut

sub adapter_reads {
    my ($self,$num_reads) = @_;
    if (defined $num_reads and $num_reads != $self->{'adapter_reads'}){
	$self->{'adapter_reads'} = $num_reads;
	$self->dirty(1);
    }
    return $self->{'adapter_reads'};
}


=head2 genotype_expected

  Arg [1]    : expected genotype in genotype checking (optional)
  Example    : my $gt_exp = $mapstats->genotype_expected();
  Description: Get/Set for expected genotype from genotype checking
  Returntype : string

=cut

sub genotype_expected {
    my ($self,$gt) = @_;
    if (defined $gt and $gt ne $self->{'gt_expected'}){
	$self->{'gt_expected'} = $gt;
	$self->dirty(1);
    }
    return $self->{'gt_expected'};
}


=head2 genotype_found

  Arg [1]    : found genotype in genotype checking (optional)
  Example    : my $gt_fnd = $mapstats->genotype_found();
  Description: Get/Set for found genotype from genotype checking
  Returntype : string

=cut

sub genotype_found {
    my ($self,$gt) = @_;
    if (defined $gt and $gt ne $self->{'gt_found'}){
	$self->{'gt_found'} = $gt;
	$self->dirty(1);
    }
    return $self->{'gt_found'};
}


=head2 genotype_ratio

  Arg [1]    : ratio between top two genotypes in genotype checking (optional)
  Example    : my $gt_ratio = $mapstats->genotype_ratio();
  Description: Get/Set for genotype ratio from genotype checking.
  Returntype : float

=cut

sub genotype_ratio {
    my ($self,$gt) = @_;
    if (defined $gt and $gt != $self->{'gt_ratio'}){
	$self->{'gt_ratio'} = $gt;
	$self->dirty(1);
    }
    return $self->{'gt_ratio'};
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


=head2 add_image_by_filename

  Arg [1]    : image file path
  Example    : my $newimage = $lib->add_image_by_filename('/tmp/gc.png');
  Description: create a new image, and if successful, return the object.  The image will be named for the name of the image file, e.g. 'gc.png' for '/tmp/gc.png'.
  Returntype : VRTrack::Image object

=cut

sub add_image_by_filename {
    my ($self, $imgfile) = @_;
    open( my $IMG, $imgfile ) or die "Can't read image $imgfile: $!\n";
    my $img = do { local( $/ ) ; <$IMG> } ;  # one-line slurp
    close $IMG;
    my $imgname = basename($imgfile);
    return $self->add_image($imgname, $img);
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


=head2 descendants

  Arg [1]    : none
  Example    : my $desc_objs = $obj->descendants();
  Description: Returns a ref to an array of all objects that are descendants of this object
  Returntype : arrayref of objects

=cut

sub descendants {
    my ($self) = @_;
    my @desc;
    foreach (@{$self->images}){
        push @desc, $_;
    }
    return \@desc;
}


1;
