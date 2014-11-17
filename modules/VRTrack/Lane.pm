package VRTrack::Lane; 

=head1 NAME

VRTrack::Lane - Sequence Tracking Lane object

=head1 SYNOPSIS
    my $lane= VRTrack::Lane->new($vrtrack, $lane_id);

    my $id = $lane->id();
    my $qc_status = $lane->qc_status();

=head1 DESCRIPTION

An object describing the tracked properties of a lane.

=head1 AUTHOR

jws@sanger.ac.uk (author)

=head1 METHODS

=cut

use strict;
use warnings;
use Carp qw(cluck confess);
use VRTrack::Mapstats;
use VRTrack::File;
use VRTrack::Submission;
use File::Spec;

use base qw(VRTrack::Core_obj
            VRTrack::Hierarchy_obj
	    VRTrack::Named_obj);


=head2 fields_dispatch

  Arg [1]    : none
  Example    : my $fieldsref = $lane->fields_dispatch();
  Description: Returns hashref dispatch table keyed on database field
               Used internally for new and update methods
  Returntype : hashref

=cut

sub fields_dispatch {
    my $self = shift;
    
    my %fields = %{$self->SUPER::fields_dispatch()};
    %fields = (%fields,
               lane_id           => sub { $self->id(@_)},
	           library_id        => sub { $self->library_id(@_)},
               seq_request_id    => sub { $self->seq_request_id(@_)},
               name              => sub { $self->name(@_)},
               hierarchy_name    => sub { $self->hierarchy_name(@_)},
               acc               => sub { $self->acc(@_)},
               readlen           => sub { $self->read_len(@_)},
               paired            => sub { $self->is_paired(@_)},
               processed         => sub { $self->processed(@_)},
               raw_reads         => sub { $self->raw_reads(@_)},
               raw_bases         => sub { $self->raw_bases(@_)},
               npg_qc_status     => sub { $self->npg_qc_status(@_)},
               qc_status         => sub { $self->qc_status(@_)},
               auto_qc_status    => sub { $self->auto_qc_status(@_)},
               gt_status         => sub { $self->genotype_status(@_)},
               storage_path      => sub { $self->storage_path(@_)},
               submission_id     => sub { $self->submission_id(@_)},
               withdrawn         => sub { $self->is_withdrawn(@_)},
               manually_withdrawn => sub { $self->is_manually_withdrawn(@_)},
               run_date          => sub { $self->run_date(@_)});
    
    return \%fields;
}


###############################################################################
# Class methods
###############################################################################

=head2 new_by_name

  Arg [1]    : vrtrack handle
  Arg [2]    : lane name
  Example    : my $lane = VRTrack::Lane->new_by_name($vrtrack, $name)
  Description: Class method. Returns latest Lane object by name.  If no such name is in the database, returns undef.  Dies if multiple names match.
  Returntype : VRTrack::Lane object

=cut


=head2 new_by_hierarchy_name

  Arg [1]    : vrtrack handle
  Arg [2]    : lane hierarchy_name
  Example    : my $lane = VRTrack::Lane->new_by_hierarchy_name($vrtrack, $hierarchy_name)
  Description: Class method. Returns latest Lane object by hierarchy_name.  If no such hierarchy_name is in the database, returns undef.  Dies if multiple hierarchy_names match.
  Returntype : VRTrack::Lane object

=cut


=head2 is_name_in_database

  Arg [1]    : lane name
  Arg [2]    : hierarchy name
  Example    : if(VRTrack::Lane->is_name_in_database($vrtrack, $name, $hname)
  Description: Class method. Checks to see if a name or hierarchy name is already used in the lane table.
  Returntype : boolean

=cut


=head2 create
  
  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : name
  Example    : my $file = VRTrack::Lane->create($vrtrack, $name)
  Description: Class method.  Creates new Lane object in the database.
  Returntype : VRTrack::Lane object
   
=cut


###############################################################################
# Object methods
###############################################################################

=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $lane->id();
	       $lane->id(104);
  Description: Get/Set for internal db ID of a lane
  Returntype : integer

=cut


=head2 library_id

  Arg [1]    : library_id (optional)
  Example    : my $library_id = $lane->library_id();
	       $lane->library_id('104');
  Description: Get/Set for ID of a lane
  Returntype : Internal ID

=cut

sub library_id {
    my $self = shift;
    return $self->_get_set('library_id', 'number', @_);
}    

=head2 seq_request_id

  Arg [1]    : seq_request_id (optional)
  Example    : my $seq_request_id = $lane->seq_request_id();
	       $lane->seq_request_id('104');
  Description: Get/Set for ID of a lane
  Returntype : Internal ID

=cut

sub seq_request_id {
    my $self = shift;
    return $self->_get_set('seq_request_id', 'number', @_);
}


=head2 hierarchy_name

  Arg [1]    : directory name (optional)
  Example    : my $hname = $lane->hierarchy_name();
  Description: Get/set lane hierarchy name.  This is the directory name (without path) that the lane will be named in a file hierarchy.
  Returntype : string

=cut


=head2 name

  Arg [1]    : name (optional)
  Example    : my $name = $lane->name();
	       $lane->name('1044_1');
  Description: Get/Set for name of a lane
  Returntype : string

=cut


=head2 acc

  Arg [1]    : acc (optional)
  Example    : my $acc = $lane->acc();
	       $lane->acc('ERR0000538');
  Description: Get/Set for [ES]RA/DCC accession
  Returntype : string

=cut

sub acc {
    my $self = shift;
    return $self->_get_set('acc', 'string', @_);
}


=head2 read_len

  Arg [1]    : read_len (optional)
  Example    : my $read_len = $lane->read_len();
	       $lane->read_len(54);
  Description: Get/Set for lane read_len
  Returntype : integer

=cut

sub read_len {
    my $self = shift;
    return $self->_get_set('read_len', 'number', @_);
}


=head2 is_paired

  Arg [1]    : boolean for is_paired status
  Example    : $lane->is_paired(1);
  Description: Get/Set for lane being paired-end sequencing
  Returntype : boolean (undef if paired status had never been set)

=cut

sub is_paired {
    my $self = shift;
    return $self->_get_set('is_paired', 'boolean', @_);
}


=head2 is_withdrawn

  Arg [1]    : boolean for is_withdrawn status
  Example    : $lane->is_withdrawn(1);
  Description: Get/Set for whether lane has been withdrawn or not
  Returntype : boolean (undef if withdrawn status had never been set)

=cut

sub is_withdrawn {
    my $self = shift;
    return $self->_get_set('is_withdrawn', 'boolean', @_);
}


=head2 is_manually_withdrawn

  Arg [1]    : boolean for is_manually_withdrawn status
  Example    : $lane->is_manually_withdrawn(1);
  Description: Get/Set for whether lane has been manually withdrawn or not;
               The distinction between this and is_withdrawn is that a lane that
               is manually withdrawn won't be automatically unwithdrawn by some
               automated system that checks this value.
  Returntype : boolean (undef if withdrawn status had never been set)

=cut

sub is_manually_withdrawn {
    my $self = shift;
    my $withdrawn = $self->_get_set('is_manually_withdrawn', 'boolean', @_);
    if (defined $withdrawn) {
        $self->is_withdrawn($withdrawn);
    }
    return $withdrawn;
}


=head2 raw_reads

  Arg [1]    : raw_reads (optional)
  Example    : my $raw_reads = $lane->raw_reads();
	       $lane->raw_reads(100000);
  Description: Get/Set for number of raw reads in lane
  Returntype : integer

=cut

sub raw_reads {
    my $self = shift;
    return $self->_get_set('raw_reads', 'number', @_);
}


=head2 raw_bases

  Arg [1]    : raw_bases (optional)
  Example    : my $raw_bases = $lane->raw_bases();
	       $lane->raw_bases(100000);
  Description: Get/Set for number of raw reads in lane
  Returntype : integer

=cut

sub raw_bases {
    my $self = shift;
    return $self->_get_set('raw_bases', 'number', @_);
}


=head2 storage_path

  Arg [1]    : storage_path (optional)
  Example    : my $storage_path = $lane->storage_path();
	       $lane->storage_path('/abs/path/to/lane');
  Description: Get/Set the absolute path to where the files associated with
               this lane are stored on disc
  Returntype : string

=cut

sub storage_path {
    my $self = shift;
    if ($_[0]) {
	confess "The supplied storage_path '@_' was not absolute" unless File::Spec->file_name_is_absolute(@_);
    }
    return $self->_get_set('storage_path', 'string', @_);
}


=head2 submission_id

  Arg [1]    : submission_id (optional)
  Example    : my $submission_id = $lane->submission_id();
	       $lane->submission_id(3);
  Description: Get/Set for submission internal id
  Returntype : integer

=cut

sub submission_id {
    my $self = shift;
    return $self->_get_set('submission_id', 'number', @_);
}


=head2 submission

  Arg [1]    : submission name (optional)
  Example    : my $submission = $lane->submission();
               $lane->submission('g1k-sc-20080812-2');
  Description: Get/Set for sample submission.  Lazy-loads submission object from $self->submission_id.  If a submission name is supplied, then submission_id is set to the corresponding submission in the database.  If no such submission exists, returns undef.  Use add_submission to add a submission in this case.
  Returntype : VRTrack::Submission object

=cut

sub submission {
    my $self = shift;
    return $self->_get_set_child_object('get_submission_by_name', 'VRTrack::Submission', @_);
}


=head2 add_submission

  Arg [1]    : submission name
  Example    : my $sub = $lane->add_submission('NA19820');
  Description: create a new submission, and if successful, return the object
  Returntype : VRTrack::Library object

=cut

sub add_submission {
    my $self = shift;
    return $self->_create_child_object('get_submission_by_name', 'VRTrack::Submission', @_);
}


=head2 get_submission_by_name

  Arg [1]    : submission_name
  Example    : my $sub = $lane->get_submission_by_name('NA19820');
  Description: Retrieve a VRTrack::Submission object by name
               Note that the submission object retrieved is not necessarily
               attached to this Lane.
  Returntype : VRTrack::Submission object

=cut

sub get_submission_by_name {
    my ($self,$name) = @_;
    return VRTrack::Submission->new_by_name($self->{vrtrack}, $name);
}


=head2 auto_qc_status

  Arg [1]    : auto_qc_status (optional)
  Example    : my $qc_status = $lane->auto_qc_status();
	       $lane->auto_qc_status('passed');
  Description: Get/Set for lane auto_qc_status
  Returntype : string

=cut

sub auto_qc_status {
    my $self = shift;
    $self->_check_status_value('auto_qc_status', @_);
    return $self->_get_set('auto_qc_status', 'string', @_);
}


=head2 qc_status

  Arg [1]    : qc_status (optional)
  Example    : my $qc_status = $lane->qc_status();
	       $lane->qc_status('passed');
  Description: Get/Set for lane qc_status
  Returntype : string

=cut

sub qc_status {
    my $self = shift;
    $self->_check_status_value('qc_status', @_);
    return $self->_get_set('qc_status', 'string', @_);
}


=head2 npg_qc_status

  Arg [1]    : npg_qc_status (optional)
  Example    : my $npg_qc_status = $lane->npg_qc_status();
	       $lane->npg_qc_status('pass');
  Description: Get/Set for lane npg_qc_status.  This is the manual QC that sequencescape/npg perform on a lane.
  Returntype : string

=cut

sub npg_qc_status {
    my $self = shift;
    $self->_check_status_value('npg_qc_status', @_);
    return $self->_get_set('npg_qc_status', 'string', @_);
}


=head2 genotype_status

  Arg [1]    : genotype_status (optional)
  Example    : my $genotype_status = $lane->genotype_status();
	       $lane->genotype_status('confirmed');
  Description: Get/Set for lane genotype check status
  Returntype : string

=cut

sub genotype_status {
    my $self = shift;
    $self->_check_status_value('gt_status', @_);
    return $self->_get_set('genotype_status', 'string', @_);
}


=head2 run_date

  Arg [1]    : run_date (optional)
  Example    : my $run_date = $lane->run_date();
               $lane->run_date('20080810123000');
  Description: Get/Set for lane run_date
  Returntype : string

=cut

sub run_date {
    my $self = shift;
    return $self->_get_set('run_date', 'string', @_);
}


=head2 changed

  Arg [1]    : changed (optional)
  Example    : my $changed = $lane->changed();
               $lane->changed('20080810123000');
  Description: Get/Set for lane changed
  Returntype : string

=cut


=head2 processed

  Arg [1]    : processed (optional)
  Description: Don't use this method, use is_processed instead.
  Returntype : string

=cut

sub processed {
    my $self = shift;
    return $self->_get_set('processed', 'number', @_);
}


=head2 is_processed

  Arg [1]    : flag, one of the flags listed in Core_obj::allowed_processed_flags();
  Arg [2]    : processed: 0 or 1 (optional)
  Example    : my $processed = $lane->is_processed('qc');
               $lane->is_processed('qc', 1);
  Description: Get/Set for lane processed.
  Returntype : boolean

=cut

sub is_processed {
    my ($self, $flag, $processed) = @_;
    
    my %flags = $self->allowed_processed_flags();
    confess "The flag '$flag' not recognised" unless exists $flags{$flag};
    
    $flag = $flags{$flag};
    if (defined $processed) {
        $processed = $processed ? $self->{processed}|$flag : $self->{processed}&(~$flag);
        if (! defined $self->{processed} || $processed != $self->{processed}) {
            $self->{processed} = $processed;
            $self->dirty(1);
        }
    }
    
    return $self->{processed} & $flag ? 1 : 0;
}


=head2 latest_mapping

  Arg [1]    : None
  Example    : my $latest_mapping = $lane->latest_mapping();
  Description: Returns single most recent mapping on this lane.
  Returntype : VRTrack::Mapstat object

=cut

sub latest_mapping {
    my ($self) = @_;
    # sort by changed date, and then by row_id
    my @sorted_mappings = sort {$a->changed cmp $b->changed 
                                || $a->row_id <=> $b->row_id} 
                                @{$self->mappings};
    return $sorted_mappings[-1];
}


=head2 mappings

  Arg [1]    : None
  Example    : my $mappings = $lane->mappings();
  Description: Returns a ref to an array of the mappings that are associated with this lane.
  Returntype : ref to array of VRTrack::Mapstats objects

=cut

sub mappings {
    my $self = shift;
    return $self->_get_child_objects('VRTrack::Mapstats');
}


=head2 mappings_excluding_qc

  Arg [1]    : None
  Example    : my $mappings = $lane->mappings_excluding_qc();
  Description: Returns a ref to an array of the mappings that are associated with this lane excluding qc mappings
  Returntype : ref to array of VRTrack::Mapstats objects

=cut

sub mappings_excluding_qc {
    my $self = shift;
    my @filtered_mappings = map { $_->is_qc() == 1 ? () : $_ } @{$self->mappings()};
    return \@filtered_mappings;
}

=head2 qc_mappings

  Arg [1]    : None
  Example    : my $mappings = $lane->qc_mappings();
  Description: Returns a ref to an array of qc mappings that are associated with this lane
  Returntype : ref to array of VRTrack::Mapstats objects

=cut

sub qc_mappings {
    my $self = shift;
    my @filtered_mappings = map { $_->is_qc() == 0 ? () : $_ } @{$self->mappings()};
    return \@filtered_mappings;
}

=head2 mapping_ids

  Arg [1]    : None
  Example    : my $mapping_ids = $lane->mapping_ids();
  Description: Returns a ref to an array of the ids of the mappings that are associated with this lane
  Returntype : ref to array of mapping ids

=cut

sub mapping_ids {
    my $self = shift;
    return $self->_get_child_ids('VRTrack::Mapstats');
}
# _get_child_objects expects methods named after the child class, so we alias.
# alias is repeated twice to avoid warning
*mapstats_ids = \&mapping_ids;
*mapstats_ids = \&mapping_ids;


=head2 add_mapping

  Arg [1]    : none
  Example    : my $newmapping = $lib->add_mapping();
  Description: create a new mapping, and if successful, return the object
  Returntype : VRTrack::Mapstats object

=cut

sub add_mapping {
    my $self = shift;
    return $self->_add_child_object(undef, 'VRTrack::Mapstats');
}


=head2 files

  Arg [1]    : None
  Example    : my $files = $lane->files();
  Description: Returns a ref to an array of the file objects that are associated with this lane.
  Returntype : ref to array of VRTrack::File objects

=cut

sub files {
    my $self = shift;
    return $self->_get_child_objects('VRTrack::File');
}


=head2 file_ids

  Arg [1]    : None
  Example    : my $file_ids = $lane->file_ids();
  Description: Returns a ref to an array of the file ids that are associated with this lane
  Returntype : ref to array of file ids

=cut

sub file_ids {
    my $self = shift;
    return $self->_get_child_ids('VRTrack::File');
}


=head2 add_file

  Arg [1]    : file name
  Example    : my $newfile = $lib->add_file('1_s_1.fastq');
  Description: create a new file, and if successful, return the object
  Returntype : VRTrack::File object

=cut

sub add_file {
    my $self = shift;
    return $self->_add_child_object('new_by_name', 'VRTrack::File', @_);
}


=head2 get_file_by_name

  Arg [1]    : file name
  Example    : my $file = $lane->get_file_by_name('My file');
  Description: retrieve file object on this lane by name
  Returntype : VRTrack::File object

=cut

sub get_file_by_name {
    my $self = shift;
    return $self->_get_child_by_field_value('files', 'name', @_);
}


=head2 get_file_by_id

  Arg [1]    : internal file id
  Example    : my $file = $lib->get_file_by_id(47);
  Description: retrieve file object by internal db file id
  Returntype : VRTrack::File object

=cut

sub get_file_by_id {
    # old behaviour was to:
    # return VRTrack::File->new($self->{vrtrack}, $id);
    my $self = shift;
    return $self->_get_child_by_field_value('files', 'id', @_);
}


=head2 descendants

  Arg [1]    : none
  Example    : my $desc_objs = $obj->descendants();
  Description: Returns a ref to an array of all objects that are descendants of this object
  Returntype : arrayref of objects

=cut

sub _get_child_methods {
    return qw(files mappings);
}

1;
