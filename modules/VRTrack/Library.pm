package VRTrack::Library; 

=head1 NAME

VRTrack::Library - Sequence Tracking Library object

=head1 SYNOPSIS
    my $lib = VRTrack::Library->new($vrtrack, $library_id);

    #get arrayref of lane objects in a library
    my $libs = $library->lanes();
    
    my $id = $library->id();
    my $name = $library->name();

=head1 DESCRIPTION

An object describing the tracked properties of a library.

=head1 AUTHOR

jws@sanger.ac.uk (author)

=head1 METHODS

=cut

use strict;
use warnings;
use Carp qw(cluck confess);
use VRTrack::Lane;
use VRTrack::Seq_request;
use VRTrack::Library_type;
use VRTrack::Seq_centre;
use VRTrack::Seq_tech;

use base qw(VRTrack::Core_obj
            VRTrack::Hierarchy_obj
	    VRTrack::Named_obj
	    VRTrack::SequenceScape_obj);


=head2 fields_dispatch

  Arg [1]    : none
  Example    : my $fieldsref = $lib->fields_dispatch();
  Description: Returns hashref dispatch table keyed on database field
               Used internally for new and update methods
  Returntype : hashref

=cut

sub fields_dispatch {
    my $self = shift;
    
    my %fields = %{$self->SUPER::fields_dispatch()};
    %fields = (%fields, 
               library_id        => sub { $self->id(@_)},
               sample_id         => sub { $self->sample_id(@_)},
               ssid              => sub { $self->ssid(@_)},
	       library_request_id        => sub { $self->library_request_id(@_)},
               hierarchy_name    => sub { $self->hierarchy_name(@_)},
               prep_status       => sub { $self->prep_status(@_)},
               qc_status         => sub { $self->qc_status(@_)},
               auto_qc_status    => sub { $self->auto_qc_status(@_)},
               fragment_size_from       => sub { $self->fragment_size_from(@_)},
	       fragment_size_to       => sub { $self->fragment_size_to(@_)},
               library_type_id   => sub { $self->library_type_id(@_)},
	       library_tag   	 => sub { $self->library_tag(@_)},
	       library_tag_group => sub { $self->library_tag_group(@_)},
	       library_tag_sequence  	 => sub { $self->library_tag_sequence(@_)},
               seq_centre_id     => sub { $self->seq_centre_id(@_)},
               seq_tech_id       => sub { $self->seq_tech_id(@_)},
               open              => sub { $self->open(@_)},
               name              => sub { $self->name(@_)});
    
    return \%fields;
}


###############################################################################
# Class methods
###############################################################################

=head2 new_by_name

  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : library name
  Example    : my $library = VRTrack::Library->new_by_name($vrtrack, $name)
  Description: Class method. Returns latest Library object by name.  If no such name is in the database, returns undef.  Dies if multiple names match.
  Returntype : VRTrack::Library object

=cut


=head2 new_by_hierarchy_name

  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : library hierarchy_name
  Example    : my $library = VRTrack::Library->new_by_hierarchy_name($vrtrack, $hierarchy_name)
  Description: Class method. Returns latest Library object by hierarchy_name.  If no such hierarchy_name is in the database, returns undef.  Dies if multiple hierarchy_names match.
  Returntype : VRTrack::Library object

=cut


=head2 new_by_ssid

  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : library sequencescape id
  Example    : my $library = VRTrack::Library->new_by_ssid($vrtrack, $ssid);
  Description: Class method. Returns latest Library object by ssid.  If no such ssid is in the database, returns undef
  Returntype : VRTrack::Library object

=cut


=head2 is_name_in_database

  Arg [1]    : library name
  Arg [2]    : hierarchy name
  Example    : if(VRTrack::Library->is_name_in_database($vrtrack, $name, $hname)
  Description: Class method. Checks to see if a name or hierarchy name is already used in the library table.
  Returntype : boolean

=cut


=head2 create
  
  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : name
  Example    : my $file = VRTrack::Library->create($vrtrack, $name)
  Description: Class method.  Creates new Library object in the database.
  Returntype : VRTrack::Library object
   
=cut


###############################################################################
# Object methods
###############################################################################

=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $lib->id();
               $lib->id(104);
  Description: Get/Set for internal db ID of a library
  Returntype : integer

=cut


=head2 sample_id

  Arg [1]    : sample_id (optional)
  Example    : my $sample_id = $lib->sample_id();
               $lib->sample_id(104);
  Description: Get/Set for ID of a library
  Returntype : Internal ID integer

=cut

sub sample_id {
    my $self = shift;
    return $self->_get_set('sample_id', 'number', @_);
}


=head2 ssid

  Arg [1]    : ssid (optional)
  Example    : my $ssid = $lib->ssid();
               $lib->ssid('104');
  Description: Get/Set for SequenceScape ID of a library
  Returntype : integer

=cut


=head2 library_request_id

  Arg [1]    : library_request_id (optional)
  Example    : my $sample_id = $library->library_request_id();
               $library->library_request_id(104);
  Description: Get/Set for ID of a library_request_id
  Returntype : Internal ID integer

=cut

sub library_request_id {
    my $self = shift;
    return $self->_get_set('library_request_id', 'number', @_);
}

=head2 hierarchy_name

  Arg [1]    : directory name (optional)
  Example    : my $hname = $library->hierarchy_name();
  Description: Get/set library hierarchy name.  This is the directory name (without path) that the library will be named in a file hierarchy.
  Returntype : string

=cut


=head2 name

  Arg [1]    : name (optional)
  Example    : my $name = $lib->name();
               $lib->name('104');
  Description: Get/Set for library name
  Returntype : string

=cut


=head2 prep_status

  Arg [1]    : prep_status (optional)
  Example    : my $prep_status = $lib->prep_status();
               $lib->prep_status('104');
  Description: Get/Set for library prep_status
  Returntype : string

=cut

sub prep_status {
    my $self = shift;
    $self->_check_status_value('prep_status', @_);
    return $self->_get_set('prep_status', 'string', @_);
}


=head2 auto_qc_status

  Arg [1]    : auto_qc_status (optional)
  Example    : my $qc_status = $lib->auto_qc_status();
	       $lib->auto_qc_status('passed');
  Description: Get/Set for library auto_qc_status
  Returntype : string

=cut

sub auto_qc_status {
    my $self = shift;
    $self->_check_status_value('auto_qc_status', @_);
    return $self->_get_set('auto_qc_status', 'string', @_);
}


=head2 qc_status

  Arg [1]    : qc_status (optional)
  Example    : my $qc_status = $lib->qc_status();
               $lib->qc_status('104');
  Description: Get/Set for library qc_status.  Checks database enum for allowed values before allowing qc_status to be set.
  Returntype : string

=cut

sub qc_status {
    my $self = shift;
    $self->_check_status_value('qc_status', @_);
    return $self->_get_set('qc_status', 'string', @_);
}

=head2 fragment_size_from

  Arg [1]    : fragment_size_from (optional)
  Example    : my $fragment_size_from = $lib->fragment_size_from();
               $lib->fragment_size_from(76);
  Description: Get/Set for library fragment_size_from
  Returntype : integer

=cut

sub fragment_size_from {
    my $self = shift;
    return $self->_get_set('fragment_size_from', 'number', @_);
}


=head2 fragment_size_to

  Arg [1]    : fragment_size_to (optional)
  Example    : my $fragment_size_to = $lib->fragment_size_to();
               $lib->fragment_size_to(76);
  Description: Get/Set for library fragment_size_to
  Returntype : integer

=cut

sub fragment_size_to {
    my $self = shift;
    return $self->_get_set('fragment_size_to', 'number', @_);
}


=head2 insert_size

  Arg [1]    : None
  Example    : my $insert_size = $lib->insert_size();
  Description: Get library insert_size.  Note that this is for backward
                compatibility only, and is the mean of fragment_size_from
                and fragment_size_to.  Use these methods if possible.
  Returntype : integer

=cut

sub insert_size {
    my $self = shift;
    my $f_from = $self->fragment_size_from();
    my $f_to = $self->fragment_size_to();
    my $insert_size = int(($f_from+$f_to)/2);
    return $insert_size;
}


=head2 library_type_id

  Arg [1]    : library_type_id (optional)
  Example    : my $library_type_id = $lib->library_type_id();
               $lib->library_type_id(1);
  Description: Get/Set for library library_type_id
  Return_type_id : string

=cut

sub library_type_id {
    my $self = shift;
    return $self->_get_set('library_type_id', 'number', @_);
}


=head2 library_tag

  Arg [1]    : library_tag (optional)
  Example    : my $library_tag_id = $lib->library_tag();
               $lib->library_tag(1);
  Description: Get/Set for library library_tag
  Return_tag : number

=cut

sub library_tag{
    my $self = shift;
    return $self->_get_set('library_tag', 'number', @_);
}

=head2 library_tag_group

  Arg [1]    : library_tag_group (optional)
  Example    : my $library_tag_group_id = $lib->library_tag_group();
               $lib->library_tag_group(1);
  Description: Get/Set for library library_tag_group
  Return_tag_group : number

=cut

sub library_tag_group{
    my $self = shift;
    return $self->_get_set('library_tag_group', 'number', @_);
}

=head2 library_tag_sequence

  Arg [1]    : library_tag_sequence (optional)
  Example    : my $library_tag_sequence_id = $lib->library_tag_sequence();
               $lib->library_tag_sequence(1);
  Description: Get/Set for library library_tag_sequence
  Return_tag_group : string

=cut

sub library_tag_sequence{
    my $self = shift;
    return $self->_get_set('library_tag_sequence', 'string', @_);
}





=head2 library_type

  Arg [1]    : library_type name (optional)
  Example    : my $library_type = $library->library_type();
               $library->library_type('DSS');
  Description: Get/Set for sample library_type.  Lazy-loads library_type object from $self->library_type_id.  If a library_type name is supplied, then library_type_id is set to the corresponding library_type in the database.  If no such library_type exists, returns undef.  Use add_library_type to add a library_type in this case.
  Returntype : VRTrack::Library_type object

=cut

sub library_type {
    my $self = shift;
    return $self->_get_set_child_object('get_library_type_by_name', 'VRTrack::Library_type', @_);
}


=head2 add_library_type

  Arg [1]    : library_type name
  Example    : my $library_type = $library->add_library_type('DSS');
  Description: create a new library_type, and if successful, return the object
  Returntype : VRTrack::Library_type object

=cut

sub add_library_type {
    my $self = shift;
    return $self->_create_child_object('get_library_type_by_name', 'VRTrack::Library_type', @_);
}


=head2 get_library_type_by_name

  Arg [1]    : library_type_name
  Example    : my $library_type = $samp->get_library_type_by_name('DSS');
  Description: Retrieve a VRTrack::Library_type object by name
               Note that the library_type object retrieved is not necessarily
               attached to this Library.  Use $lib->library_type for that.
  Returntype : VRTrack::Seq_centre object
  Returntype : VRTrack::Library_type object

=cut

sub get_library_type_by_name {
    my ($self,$name) = @_;
    return VRTrack::Library_type->new_by_name($self->{vrtrack}, $name);
}

=head2 sequencing requests

  Arg [1]    : None
  Example    : my $seq_requests = $library->seq_requests();
  Description: Returns a ref to an array of the Seq_Request objects that are associated with this sample.
  Returntype : ref to array of VRTrack::Seq_Request objects

=cut

sub seq_requests {
    my $self = shift;
    return $self->_get_child_objects('VRTrack::Seq_request');
}


=head2 seq_request_ids

  Arg [1]    : None
  Example    : my $seq_request_ids = $library->seq_request_ids();
  Description: Returns a ref to an array of the seq_request IDs that are associated with this sample
  Returntype : ref to array of integer seq_request IDs

=cut

sub seq_request_ids {
    my $self = shift;
    return $self->_get_child_ids('VRTrack::Seq_request');
}


=head2 add_seq_request

  Arg [1]    : sequencing request id
  Example    : my $newseqrequest = $samp->add_seq_request('12121212');
  Description: create a new seq_request , and if successful, return the object
  Returntype : VRTrack::Seq_request object

=cut

sub add_seq_request {
    my ($self,$ssid) = @_;
    return $self->_add_child_object('new_by_ssid','VRTrack::Seq_request',$ssid);
}


=head2 get_seq_request_by_id

  Arg [1]    : seq_request internal id
  Example    : my $lseq_request = $sam->get_lseq_request_by_id(1930);
  Description: retrieve seq_request  object by internal id
  Returntype : VRTrack::seq_request object

=cut

sub get_seq_request_by_id {
    my $self = shift;
    return $self->_get_child_by_field_value('seq_requests', 'id', @_);
}


=head2 get_seq_request_by_ssid

  Arg [1]    : seq_request sequencescape id
  Example    : my $seq_request = $sam->get_seq_request_by_ssid(1930);
  Description: retrieve seq_request object by sequencescape id
  Returntype : VRTrack::Seq_request object

=cut

sub get_seq_request_by_ssid {
    my $self = shift;
    return $self->_get_child_by_field_value('seq_requests', 'ssid', @_);
}

=head2 seq_centre_id

  Arg [1]    : seq_centre_id (optional)
  Example    : my $seq_centre_id = $lib->seq_centre_id();
               $lib->seq_centre_id(1);
  Description: Get/Set for library sequencing seq_centre_id
  Returntype : string

=cut

sub seq_centre_id {
    my $self = shift;
    return $self->_get_set('seq_centre_id', 'number', @_);
}


=head2 seq_centre

  Arg [1]    : seq_centre name (optional)
  Example    : my $seq_centre = $library->seq_centre();
               $library->seq_centre('SC');
  Description: Get/Set for sample seq_centre.  Lazy-loads seq_centre object from $self->seq_centre_id.  If a seq_centre name is supplied, then seq_centre_id is set to the corresponding seq_centre in the database.  If no such seq_centre exists, returns undef.  Use add_seq_centre to add a seq_centre in this case.
  Returntype : VRTrack::Seq_centre object

=cut

sub seq_centre {
    my $self = shift;
    return $self->_get_set_child_object('get_seq_centre_by_name', 'VRTrack::Seq_centre', @_);
}


=head2 add_seq_centre

  Arg [1]    : seq_centre name
  Example    : my $seq_centre = $library->add_seq_centre('SC');
  Description: create a new seq_centre, and if successful, return the object
  Returntype : VRTrack::Library object

=cut

sub add_seq_centre {
    my $self = shift;
    return $self->_create_child_object('get_seq_centre_by_name', 'VRTrack::Seq_centre', @_);
}


=head2 get_seq_centre_by_name

  Arg [1]    : seq_centre_name
  Example    : my $seq_centre = $library->get_seq_centre_by_name('SC');
  Description: Retrieve a VRTrack::Seq_centre object by name
               Note that the seq_centre object retrieved is not necessarily
               attached to this Library.  Use $lib->seq_centre for that.
  Returntype : VRTrack::Seq_centre object

=cut

sub get_seq_centre_by_name {
    my ($self,$name) = @_;
    return VRTrack::Seq_centre->new_by_name($self->{vrtrack}, $name);
}


=head2 seq_tech_id

  Arg [1]    : seq_tech_id (optional)
  Example    : my $seq_tech_id = $lib->seq_tech_id();
               $lib->seq_tech_id(4);
  Description: Get/Set for library seq_tech_id
  Returntype : string

=cut

sub seq_tech_id {
    my $self = shift;
    return $self->_get_set('seq_tech_id', 'number', @_);
}


=head2 seq_tech

  Arg [1]    : seq_tech name (optional)
  Example    : my $seq_tech = $lib->seq_tech();
               $lib->seq_tech('SLX');
  Description: Get/Set for sample seq_tech.  Lazy-loads seq_tech object from $self->seq_tech_id.  If a seq_tech name is supplied, then seq_tech_id is set to the corresponding seq_tech in the database.  If no such seq_tech exists, returns undef.  Use add_seq_tech to add a seq_tech in this case.
  Returntype : VRTrack::Seq_tech object

=cut

sub seq_tech {
    my $self = shift;
    return $self->_get_set_child_object('get_seq_tech_by_name', 'VRTrack::Seq_tech', @_);
}


=head2 add_seq_tech

  Arg [1]    : seq_tech name
  Example    : my $seq_tech = $lib->add_seq_tech('SLX');
  Description: create a new seq_tech, and if successful, return the object
  Returntype : VRTrack::Library object

=cut

sub add_seq_tech {
    my $self = shift;
    return $self->_create_child_object('get_seq_tech_by_name', 'VRTrack::Seq_tech', @_);
}


=head2 get_seq_tech_by_name

  Arg [1]    : seq_tech_name
  Example    : my $seq_tech = $lib->get_seq_tech_by_name('SLX');
  Description: Retrieve a VRTrack::Seq_tech object by name
               Note that the seq_tech object retrieved is not necessarily
               attached to this Library.  Use $lib->seq_tech for that.
  Returntype : VRTrack::Seq_tech object

=cut

sub get_seq_tech_by_name {
    my ($self,$name) = @_;
    return VRTrack::Seq_tech->new_by_name($self->{vrtrack}, $name);
}


=head2 open

  Arg [1]    : open (optional 0 or 1)
  Example    : my $open = $lib->open();
               $lib->open(1);
  Description: Get/Set for whether library is open for sequencing
  Returntype : boolean

=cut

sub open {
    my $self = shift;
    return $self->_get_set('open', 'boolean', @_);
}


=head2 changed

  Arg [1]    : timestamp (optional)
  Example    : my $changed = $lib->changed();
               $lib->changed('20080810123000');
  Description: Get/Set for library changed
  Returntype : string

=cut


=head2 library_multiplex_pools

  Arg [1]    : None
  Example    : my $library_multiplex_pools = $library->library_multiplex_pools();
  Description: Returns a ref to an array of the library_multiplex_pool objects that are associated with this library.
  Returntype : ref to array of VRTrack::Library_Multiplex_pool objects

=cut

sub library_multiplex_pools {
    my $self = shift;
    return $self->_get_child_objects('VRTrack::Library_Multiplex_pool');
}


=head2 library_multiplex_pool_ids

  Arg [1]    : None
  Example    : my $library_multiplex_pool_ids = $library->library_multiplex_pool_ids();
  Description: Returns a ref to an array of the library_multiplex_pool IDs that are associated with this library
  Returntype : ref to array of integer library_multiplex_pool IDs

=cut

sub library_multiplex_pool_ids {
    my $self = shift;
    return $self->_get_child_ids('VRTrack::Library_Multiplex_pool');
}



=head2 lanes

  Arg [1]    : None
  Example    : my $lanes = $library->lanes();
  Description: Returns a ref to an array of the lane objects that are associated with this library.
  Returntype : ref to array of VRTrack::Lane objects

=cut

sub lanes {
    my $self = shift;
    return $self->_get_child_objects('VRTrack::Lane');
}


=head2 lane_ids

  Arg [1]    : None
  Example    : my $lane_ids = $library->lane_ids();
  Description: Returns a ref to an array of the lane IDs that are associated with this library
  Returntype : ref to array of integer lane IDs

=cut

sub lane_ids {
    my $self = shift;
    return $self->_get_child_ids('VRTrack::Lane');
}


=head2 add_lane

  Arg [1]    : lane name
  Example    : my $newlane = $lib->add_lane('2631_3');
  Description: create a new lane, and if successful, return the object
  Returntype : VRTrack::Lane object

=cut

sub add_lane {
    my $self = shift;
    return $self->_add_child_object('new_by_name', 'VRTrack::Lane', @_);
}


=head2 get_lane_by_id

  Arg [1]    : internal lane id
  Example    : my $lane = $lib->get_lane_by_id(47);
  Description: retrieve lane object by internal db lane id
  Returntype : VRTrack::Lane object

=cut

sub get_lane_by_id {
    my $self = shift;
    return $self->_get_child_by_field_value('lanes', 'id', @_);
}


=head2 get_lane_by_name

  Arg [1]    : lane name
  Example    : my $lane = $track->get_lane_by_name('My lane');
  Description: retrieve lane object attached to this Library, by name
  Returntype : VRTrack::Lane object

=cut

sub get_lane_by_name {
    my $self = shift;
    return $self->_get_child_by_field_value('lanes', 'name', @_);
}


=head2 descendants

  Arg [1]    : none
  Example    : my $desc_objs = $obj->descendants();
  Description: Returns a ref to an array of all objects that are descendants of this object
  Returntype : arrayref of objects

=cut

sub _get_child_methods {
    return qw(library_multiplex_pools seq_requests lanes);
}


=head2 projected_passed_depth

  Arg [1]    : ref bp
  Example    : my $estDepth = $lib->projected_passed_depth(3000000000);
  Description: Returns the estimated mapped depth based on the qc_status of the lib lanes and latest mapstats
  Returntype : float

=cut

sub projected_passed_depth {
    my ($self, $ref_size) = @_;
    
    confess( "Size of reference must be a number!" ) unless $ref_size =~ /^\d+$/;
    
    my $lanes = $self->lanes();
    my $cum_depth = 0;
    foreach( @$lanes )
    {
        my $lane = $_;
        if( $lane->qc_status() eq 'passed' )
        {
            my $mapping = $lane->latest_mapping();
            if( $mapping )
            {
                $cum_depth += ( ( $mapping->rmdup_bases_mapped / $mapping->raw_bases() ) * $lane->raw_bases() ) /$ref_size;
            }
        }
    }
    
    $cum_depth = sprintf("%.2f", $cum_depth);
    return $cum_depth;
}

1;
