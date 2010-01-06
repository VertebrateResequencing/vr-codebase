package VRTrack::Mapstats; 

=head1 NAME

VRTrack::Mapstats - Sequence Tracking Mapstats object

=head1 SYNOPSIS
    my $mapstats= VRTrack::Mapstats->new($vrtrack, $mapstats_id);

    my $id = $mapstats->id();
    my $qc_status = $mapstats->qc_status();

=head1 DESCRIPTION

An object describing a set of statistics for a specific lane mapping.

=head1 CONTACT

jws@sanger.ac.uk (author)

=head1 METHODS

=cut

use strict;
use warnings;
use Carp qw(cluck confess);
use VRTrack::Image;
use VRTrack::Mapper;
use VRTrack::Assembly;
use File::Basename;

use base qw(VRTrack::Core_obj);


=head2 fields_dispatch

  Arg [1]    : none
  Example    : my $fieldsref = $mapstats->fields_dispatch();
  Description: Returns hashref dispatch table keyed on database field
               Used internally for new and update methods
  Returntype : hashref

=cut

sub fields_dispatch {
    my $self = shift;
    
    my %fields = %{$self->SUPER::fields_dispatch()};
    %fields = (%fields,
	       mapstats_id        => sub { $self->id(@_)},
               lane_id            => sub { $self->lane_id(@_)},
               mapper_id          => sub { $self->mapper_id(@_)},
               assembly_id        => sub { $self->assembly_id(@_)},
               raw_reads          => sub { $self->raw_reads(@_)},
               raw_bases          => sub { $self->raw_bases(@_)},
               clip_bases         => sub { $self->clip_bases(@_)},
               reads_mapped       => sub { $self->reads_mapped(@_)},
               reads_paired       => sub { $self->reads_paired(@_)},
               bases_mapped       => sub { $self->bases_mapped(@_)},
               rmdup_reads_mapped => sub { $self->rmdup_reads_mapped(@_)},
               rmdup_bases_mapped => sub { $self->rmdup_bases_mapped(@_)},
               error_rate         => sub { $self->error_rate(@_)},
               mean_insert        => sub { $self->mean_insert(@_)},
               sd_insert          => sub { $self->sd_insert(@_)},
               adapter_reads      => sub { $self->adapter_reads(@_)},
               gt_expected        => sub { $self->genotype_expected(@_)},
               gt_found           => sub { $self->genotype_found(@_)},
               gt_ratio           => sub { $self->genotype_ratio(@_)});

    return \%fields;
}


###############################################################################
# Class methods
###############################################################################

=head2 create
   
     Arg [1]    : database handle to seqtracking database
     Arg [2]    : lane id that this mapping relates to
     Example    : my $mapstats = VRTrack::Mapstats->create($vrtrack, $lane_id)
     Description: Class method.  Creates new Mapstats object in the database.
     Returntype : VRTrack::Mapstats object
   
=cut


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


=head2 lane_id

  Arg [1]    : lane_id (optional)
  Example    : my $lane_id = $mapstats->lane_id();
	       $mapstats->lane_id(104);
  Description: Get/Set for internal ID of the lane which was mapped
  Returntype : integer

=cut

sub lane_id {
    my $self = shift;
    return $self->_get_set('lane_id', 'number', @_);
}


=head2 mapper_id

  Arg [1]    : mapper_id (optional)
  Example    : my $mapper_id = $mapstats->mapper_id();
               $mapstats->mapper_id(104);
  Description: Get/Set for internal ID of the mapper used for the mapping
  Returntype : integer

=cut

sub mapper_id {
    my $self = shift;
    return $self->_get_set('mapper_id', 'number', @_);
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
    my $self = shift;
    return $self->_get_set_child_object('get_mapper_by_name_version', 'VRTrack::Mapper', @_);
}


=head2 add_mapper

  Arg [1]    : mapper name
  Arg [2]    : mapper version
  Example    : my $map = $mapstats->add_mapper('maq','0.7.1-6');
  Description: create a new mapper, and if successful, return the object
  Returntype : VRTrack::Library object

=cut

sub add_mapper {
    my $self = shift;
    return $self->_create_child_object('get_mapper_by_name_version', 'VRTrack::Mapper', @_);
}


=head2 get_mapper_by_name_version

  Arg [1]    : mapper_name
  Arg [2]    : mapper version
  Example    : my $map = $mapstats->get_mapper_by_name_version('maq','0.7.1-6');
  Description: Retrieve a VRTrack::Mapper object by name
  Returntype : VRTrack::Mapper object

=cut

sub get_mapper_by_name_version {
    my ($self, $name, $version) = @_;
    return VRTrack::Mapper->new_by_name_version($self->{vrtrack}, $name, $version);
}


=head2 assembly_id

  Arg [1]    : assembly_id (optional)
  Example    : my $assembly_id = $mapstats->assembly_id();
               $mapstats->assembly_id(104);
  Description: Get/Set for internal ID of the assembly used for the mapping
  Returntype : integer

=cut

sub assembly_id {
    my $self = shift;
    return $self->_get_set('assembly_id', 'number', @_);
}


=head2 assembly

  Arg [1]    : assembly name (optional)
  Example    : my $assembly = $mapstats->assembly();
               $mapstats->assembly('NCBIm36');
  Description: Get/Set for mapping assembly.  Lazy-loads assembly object from $self->assembly_id.  If a assembly name is supplied, then assembly_id is set to the corresponding assembly in the database.  If no such assembly exists, returns undef.  Use add_assembly to add a assembly in this case.
  Returntype : VRTrack::Assembly object

=cut

sub assembly {
    my $self = shift;
    return $self->_get_set_child_object('get_assembly_by_name', 'VRTrack::Assembly', @_);
}


=head2 add_assembly

  Arg [1]    : assembly name
  Example    : my $assy = $mapping->add_assembly('NCBIm36');
  Description: create a new assembly, and if successful, return the object
  Returntype : VRTrack::Library object

=cut

sub add_assembly {
    my $self = shift;
    return $self->_create_child_object('get_assembly_by_name', 'VRTrack::Assembly', @_);
}


=head2 get_assembly_by_name

  Arg [1]    : assembly_name
  Example    : my $map = $mapstats->get_assembly_by_name_version('maq');
  Description: Retrieve a VRTrack::Assembly object by name
  Returntype : VRTrack::Assembly object

=cut

sub get_assembly_by_name {
    my ($self,$name,$version) = @_;
    return VRTrack::Assembly->new_by_name($self->{vrtrack}, $name);
}


=head2 raw_reads

  Arg [1]    : number of raw reads used in mapping (optional)
  Example    : my $num_reads = $mapstats->raw_reads();
	       $mapstats->raw_reads(1_000_000);
  Description: Get/Set for total number of raw reads used in mapping
  Returntype : integer

=cut

sub raw_reads {
    my $self = shift;
    return $self->_get_set('raw_reads', 'number', @_);
}


=head2 raw_bases

  Arg [1]    : number of raw bases used in mapping (optional)
  Example    : my $num_bases = $mapstats->raw_bases();
	       $mapstats->raw_bases(1_000_000);
  Description: Get/Set for total number of raw bases used in mapping
  Returntype : integer

=cut

sub raw_bases {
    my $self = shift;
    return $self->_get_set('raw_bases', 'number', @_);
}


=head2 clip_bases

  Arg [1]    : number of raw bases after clipping
  Example    : my $bases_after_clipping = $mapstats->clip_bases();
	       $mapstats->clip_bases(1_000_000);
  Description: Get/Set for total number of raw bases after clipping
  Returntype : integer

=cut

sub clip_bases {
    my $self = shift;
    return $self->_get_set('clip_bases', 'number', @_);
}


=head2 reads_mapped

  Arg [1]    : number of mapped reads in mapping (optional)
  Example    : my $num_reads = $mapstats->reads_mapped();
	       $mapstats->reads_mapped(1_000_000);
  Description: Get/Set for total number of mapped reads in mapping
  Returntype : integer

=cut

sub reads_mapped {
    my $self = shift;
    return $self->_get_set('reads_mapped', 'number', @_);
}


=head2 reads_paired

  Arg [1]    : number of paired mapped reads in mapping (optional)
  Example    : my $num_reads = $mapstats->reads_paired();
	       $mapstats->reads_paired(1_000_000);
  Description: Get/Set for number of paired mapped reads in mapping
  Returntype : integer

=cut

sub reads_paired {
    my $self = shift;
    return $self->_get_set('reads_paired', 'number', @_);
}


=head2 bases_mapped

  Arg [1]    : number of mapped bases in mapping (optional)
  Example    : my $num_bases = $mapstats->bases_mapped();
	       $mapstats->bases_mapped(1_000_000);
  Description: Get/Set for total number of mapped bases in mapping
  Returntype : integer

=cut

sub bases_mapped {
    my $self = shift;
    return $self->_get_set('bases_mapped', 'number', @_);
}


=head2 rmdup_reads_mapped

  Arg [1]    : number of mapped reads in mapping (optional)
  Example    : my $num_reads = $mapstats->rmdup_reads_mapped();
	       $mapstats->rmdup_reads_mapped(1_000_000);
  Description: Get/Set for number of mapped reads in mapping after removing duplicates
  Returntype : integer

=cut

sub rmdup_reads_mapped {
    my $self = shift;
    return $self->_get_set('rmdup_reads_mapped', 'number', @_);
}


=head2 rmdup_bases_mapped

  Arg [1]    : number of mapped bases in mapping (optional)
  Example    : my $num_bases = $mapstats->rmdup_bases_mapped();
	       $mapstats->rmdup_bases_mapped(1_000_000);
  Description: Get/Set for number of mapped bases in mapping after removing duplicate reads
  Returntype : integer

=cut

sub rmdup_bases_mapped {
    my $self = shift;
    return $self->_get_set('rmdup_bases_mapped', 'number', @_);
}


=head2 error_rate

  Arg [1]    : maq/mapping error rate (optional)
  Example    : my $err = $mapstats->error_rate();
	       $mapstats->error_rate(0.01);
  Description: Get/Set for mapper error rate in mapping
  Returntype : float

=cut

sub error_rate {
    my $self = shift;
    return $self->_get_set('error_rate', 'number', @_);
}


=head2 mean_insert

  Arg [1]    : maq/mapping insert size (optional)
  Example    : my $ins = $mapstats->mean_insert();
	       $mapstats->mean_insert(243.4);
  Description: Get/Set for mean mapped insert size for pairs
  Returntype : float

=cut

sub mean_insert {
    my $self = shift;
    return $self->_get_set('mean_insert', 'number', @_);
}


=head2 sd_insert

  Arg [1]    : maq/mapping insert sd (optional)
  Example    : my $ins = $mapstats->sd_insert();
	       $mapstats->sd_insert(15.64);
  Description: Get/Set for Standard Deviation of mapped insert size for pairs
  Returntype : float

=cut

sub sd_insert {
    my $self = shift;
    return $self->_get_set('sd_insert', 'number', @_);
}


=head2 adapter_reads

  Arg [1]    : number of reads containing adapter sequence
  Example    : my $num_adap_reads = $mapstats->adapter_reads();
	       $mapstats->adapter_reads(1_000_000);
  Description: Get/Set for total number of reads containing adapter sequence
  Returntype : integer

=cut

sub adapter_reads {
    my $self = shift;
    return $self->_get_set('adapter_reads', 'number', @_);
}


=head2 genotype_expected

  Arg [1]    : expected genotype in genotype checking (optional)
  Example    : my $gt_exp = $mapstats->genotype_expected();
  Description: Get/Set for expected genotype from genotype checking
  Returntype : string

=cut

sub genotype_expected {
    my $self = shift;
    return $self->_get_set('gt_expected', 'string', @_);
}


=head2 genotype_found

  Arg [1]    : found genotype in genotype checking (optional)
  Example    : my $gt_fnd = $mapstats->genotype_found();
  Description: Get/Set for found genotype from genotype checking
  Returntype : string

=cut

sub genotype_found {
    my $self = shift;
    return $self->_get_set('gt_found', 'string', @_);
}


=head2 genotype_ratio

  Arg [1]    : ratio between top two genotypes in genotype checking (optional)
  Example    : my $gt_ratio = $mapstats->genotype_ratio();
  Description: Get/Set for genotype ratio from genotype checking.
  Returntype : float

=cut

sub genotype_ratio {
    my $self = shift;
    return $self->_get_set('gt_ratio', 'number', @_);
}


=head2 add_image_by_filename

  Arg [1]    : image file path
  Example    : my $newimage = $lib->add_image_by_filename('/tmp/gc.png');
  Description: create a new image, and if successful, return the object.  The image will be named for the name of the image file, e.g. 'gc.png' for '/tmp/gc.png'.
  Returntype : VRTrack::Image object

=cut

sub add_image_by_filename {
    my ($self, $imgfile) = @_;
    open( my $IMG, $imgfile ) or confess "Can't read image $imgfile: $!\n";
    my $img = do { local( $/ ) ; <$IMG> } ;  # one-line slurp
    close $IMG;
    my $imgname = basename($imgfile);
    return $self->add_image($imgname, $img);
}


=head2 add_image

  Arg [1]    : image name
  Arg [2]    : image data (i.e. the binary image)
  Example    : my $newimage = $lib->add_image('gc.png',$gc_png_img);
  Description: create a new image, and if successful, return the object
  Returntype : VRTrack::Image object

=cut

sub add_image {
    my $self = shift;
    return $self->_add_child_object(undef, 'VRTrack::Image', @_);
}


=head2 images

  Arg [1]    : None
  Example    : my $images = $mapping->images();
  Description: Returns a ref to an array of the image objects that are associated with this mapping.
  Returntype : ref to array of VRTrack::Image objects

=cut

sub images {
    my $self = shift;
    return $self->_get_child_objects('VRTrack::Image');
}


=head2 image_ids

  Arg [1]    : None
  Example    : my $image_ids = $mapping->image_ids();
  Description: Returns a ref to an array of the image ids that are associated with this mapping.
  Returntype : ref to array of image ids

=cut

sub image_ids {
    my $self = shift;
    return $self->_get_child_ids('VRTrack::Image');
}


=head2 changed

  Arg [1]    : changed (optional)
  Example    : my $changed = $mapstats->changed();
	       $mapstats->changed('20090101102300');
  Description: Get/Set for mapstats changed
  Returntype : string

=cut


=head2 descendants

  Arg [1]    : none
  Example    : my $desc_objs = $obj->descendants();
  Description: Returns a ref to an array of all objects that are descendants of this object
  Returntype : arrayref of objects

=cut

sub _get_child_methods {
    return qw(images);
}

1;
