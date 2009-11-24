package VRTrack::Lane; 
# author: jws
=head1 NAME

VRTrack::Lane - Sequence Tracking Lane object

=head1 SYNOPSIS
    my $lane= VRTrack::Lane->new($vrtrack, $lane_id);

    my $id = $lane->id();
    my $status = $lane->status();

=head1 DESCRIPTION

An object describing the tracked properties of a lane.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;
use Carp qw(cluck confess carp croak);
no warnings 'uninitialized';
use VRTrack::Mapstats;
use VRTrack::File;
use VRTrack::Submission;
use VRTrack::Core_obj;
use File::Spec;
our @ISA = qw(VRTrack::Core_obj);

=head2 fields_dispatch

  Arg [1]    : none
  Example    : my $fieldsref = $lane->fields_dispatch();
  Description: Returns hashref dispatch table keyed on database field
               Used internally for new and update methods
  Returntype : hashref

=cut

sub fields_dispatch {
    my $self = shift;
    my %fields = ( 
                'lane_id'           => sub { $self->id(@_)},
                'library_id'        => sub { $self->library_id(@_)},
                'name'              => sub { $self->name(@_)},
                'hierarchy_name'    => sub { $self->hierarchy_name(@_)},
                'acc'               => sub { $self->acc(@_)},
                'readlen'           => sub { $self->read_len(@_)},
                'paired'            => sub { $self->is_paired(@_)},
                'processed'         => sub { $self->processed(@_)},
                'raw_reads'         => sub { $self->raw_reads(@_)},
                'raw_bases'         => sub { $self->raw_bases(@_)},
                'npg_qc_status'     => sub { $self->npg_qc_status(@_)},
                'qc_status'         => sub { $self->qc_status(@_)},
                'auto_qc_status'    => sub { $self->auto_qc_status(@_)},
                'gt_status'         => sub { $self->genotype_status(@_)},
                'submission_id'     => sub { $self->submission_id(@_)},
                'withdrawn'         => sub { $self->is_withdrawn(@_)},
                'note_id'           => sub { $self->note_id(@_)},
                'changed'           => sub { $self->changed(@_)},
                'run_date'          => sub { $self->run_date(@_)},
                'latest'            => sub { $self->is_latest(@_)},
                );

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

sub new_by_name {
    my ($class,$vrtrack, $name) = @_;
    die "Need to call with a handle, name" unless ($vrtrack && $name);
    return $class->new_by_field_value($vrtrack, 'name',$name);
}


=head2 new_by_hierarchy_name

  Arg [1]    : vrtrack handle
  Arg [2]    : lane hierarchy_name
  Example    : my $lane = VRTrack::Lane->new_by_hierarchy_name($vrtrack, $hierarchy_name)
  Description: Class method. Returns latest Lane object by hierarchy_name.  If no such hierarchy_name is in the database, returns undef.  Dies if multiple hierarchy_names match.
  Returntype : VRTrack::Lane object

=cut

sub new_by_hierarchy_name {
    my ($class,$vrtrack, $hierarchy_name) = @_;
    die "Need to call with a vrtrack handle, hierarchy_name" unless ($vrtrack && $hierarchy_name);
    return $class->new_by_field_value($vrtrack, 'hierarchy_name',$hierarchy_name);
}


#   =head2 create
#   
#     Arg [1]    : vrtrack handle
#     Arg [2]    : lane name
#     Example    : my $lane = VRTrack::Lane->create($vrtrack, $name)
#     Description: Class method.  Creates new Lane object in the database.
#     Returntype : VRTrack::Lane object
#   
#   =cut
#   
#   sub create {
#       my ($class,$vrtrack, $name) = @_;
#       die "Need to call with a vrtrack handle and name" unless ($vrtrack && $name);
#       if ( $vrtrack->isa('DBI::db') ) { croak "The interface has changed, expected vrtrack reference.\n"; }
#       my $dbh = $vrtrack->{_dbh};
#   
#       my $hierarchy_name = $name;
#       $hierarchy_name =~ s/\W+/_/g;
#   
#       # prevent adding a lane with an existing name
#       if ($class->is_name_in_database($vrtrack, $name, $hierarchy_name)){
#           die "Already a lane by name $name/$hierarchy_name";
#       }
#   
#       $dbh->do (qq[LOCK TABLE lane WRITE]);
#       my $sql = qq[select max(lane_id) as id from lane];
#       my $sth = $dbh->prepare($sql);
#       my $next_id;
#       if ($sth->execute()){
#   	my $data = $sth->fetchrow_hashref;
#   	unless ($data){
#               $dbh->do (qq[UNLOCK TABLES]);
#               die( sprintf("Can't retrieve next lane id: %s", $DBI::errstr));
#   	}
#           $next_id = $data->{'id'};
#           $next_id++;
#       }
#       else{
#   	die(sprintf("Can't retrieve next lane id: %s", $DBI::errstr));
#       }
#   
#       $sql = qq[INSERT INTO lane (lane_id, name, hierarchy_name, changed, latest) 
#                    VALUES (?,?,?,now(),true)];
#   
#       $sth = $dbh->prepare($sql);
#       unless ($sth->execute( $next_id, $name,$hierarchy_name )) {
#           $dbh->do (qq[UNLOCK TABLES]);
#           die( sprintf('DB load insert failed: %s %s', $next_id, $DBI::errstr));
#       }
#   
#       $dbh->do (qq[UNLOCK TABLES]);
#   
#       return $class->new($vrtrack, $next_id);
#   }


=head2 is_name_in_database

  Arg [1]    : lane name
  Arg [2]    : hierarchy name
  Example    : if(VRTrack::Lane->is_name_in_database($vrtrack, $name, $hname)
  Description: Class method. Checks to see if a name or hierarchy name is already used in the lane table.
  Returntype : boolean

=cut

sub is_name_in_database {
    my ($class, $vrtrack, $name, $hname) = @_;
    die "Need to call with a vrtrack handle, name, hierarchy name" unless ($vrtrack && $name && $hname);
    if ( $vrtrack->isa('DBI::db') ) { croak "The interface has changed, expected vrtrack reference.\n"; }
    my $dbh = $vrtrack->{_dbh};
    my $sql = qq[select lane_id from lane where latest=true and (name = ? or hierarchy_name = ?) ];
    my $sth = $dbh->prepare($sql);

    my $already_used = 0;
    if ($sth->execute($name,$hname)){
        my $data = $sth->fetchrow_hashref;
        if ($data){
            $already_used = 1;
        }
    }
    else{
        die(sprintf('Cannot retrieve lane by $name: %s', $DBI::errstr));
    }
    return $already_used;
}


###############################################################################
# Object methods
###############################################################################

=head2 dirty

  Arg [1]    : boolean for dirty status
  Example    : $lane->dirty(1);
  Description: Get/Set for lane properties having been altered.
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
  Example    : my $id = $lane->id();
	       $lane->id(104);
  Description: Get/Set for internal db ID of a lane
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


=head2 library_id

  Arg [1]    : library_id (optional)
  Example    : my $library_id = $lane->library_id();
	       $lane->library_id('104');
  Description: Get/Set for ID of a lane
  Returntype : Internal ID

=cut

sub library_id {
    my ($self,$library_id) = @_;
    if (defined $library_id and $library_id != $self->{'library_id'}){
	$self->{'library_id'} = $library_id;
	$self->dirty(1);
    }
    return $self->{'library_id'};
}


=head2 hierarchy_name

  Arg [1]    : directory name (optional)
  Example    : my $hname = $lane->hierarchy_name();
  Description: Get/set lane hierarchy name.  This is the directory name (without path) that the lane will be named in a file hierarchy.
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
  Example    : my $name = $lane->name();
	       $lane->name('1044_1');
  Description: Get/Set for name of a lane
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


=head2 acc

  Arg [1]    : acc (optional)
  Example    : my $acc = $lane->acc();
	       $lane->acc('ERR0000538');
  Description: Get/Set for [ES]RA/DCC accession
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



=head2 read_len

  Arg [1]    : read_len (optional)
  Example    : my $read_len = $lane->read_len();
	       $lane->read_len(54);
  Description: Get/Set for lane read_len
  Returntype : integer

=cut

sub read_len {
    my ($self,$read_len) = @_;
    if (defined $read_len and $read_len != $self->{'read_len'}){
	$self->{'read_len'} = $read_len;
	$self->dirty(1);
    }
    return $self->{'read_len'};
}


=head2 is_paired

  Arg [1]    : boolean for is_paired status
  Example    : $lane->is_paired(1);
  Description: Get/Set for lane being paired-end sequencing
  Returntype : boolean (undef if paired status had never been set)

=cut

sub is_paired {
    my ($self,$is_paired) = @_;
    if (defined $is_paired){
	$self->{is_paired} = $is_paired ? 1 : 0;
	$self->dirty(1);
    }
    return $self->{is_paired};
}


=head2 is_withdrawn

  Arg [1]    : boolean for is_withdrawn status
  Example    : $lane->is_withdrawn(1);
  Description: Get/Set for whether lane has been withdrawn or not
  Returntype : boolean (undef if withdrawn status had never been set)

=cut

sub is_withdrawn {
    my ($self,$is_withdrawn) = @_;
    if (defined $is_withdrawn){
	$self->{is_withdrawn} = $is_withdrawn ? 1 : 0;
	$self->dirty(1);
    }
    return $self->{is_withdrawn};
}


=head2 raw_reads

  Arg [1]    : raw_reads (optional)
  Example    : my $raw_reads = $lane->raw_reads();
	       $lane->raw_reads(100000);
  Description: Get/Set for number of raw reads in lane
  Returntype : integer

=cut

sub raw_reads {
    my ($self,$raw_reads) = @_;
    if (defined $raw_reads and $raw_reads != $self->{'raw_reads'}){
	$self->{'raw_reads'} = $raw_reads;
	$self->dirty(1);
    }
    return $self->{'raw_reads'};
}


=head2 raw_bases

  Arg [1]    : raw_bases (optional)
  Example    : my $raw_bases = $lane->raw_bases();
	       $lane->raw_bases(100000);
  Description: Get/Set for number of raw reads in lane
  Returntype : integer

=cut

sub raw_bases {
    my ($self,$raw_bases) = @_;
    if (defined $raw_bases and $raw_bases != $self->{'raw_bases'}){
	$self->{'raw_bases'} = $raw_bases;
	$self->dirty(1);
    }
    return $self->{'raw_bases'};
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
    my ($self, $storage_path) = @_;
    if (defined $storage_path and $storage_path != $self->{'storage_path'}){
	die "The supplied storage_path '$storage_path' was not absolute" unless File::Spec->file_name_is_absolute($storage_path);
	$self->{'storage_path'} = $storage_path;
	$self->dirty(1);
    }
    return $self->{'storage_path'};
}


=head2 submission_id

  Arg [1]    : submission_id (optional)
  Example    : my $submission_id = $lane->submission_id();
	       $lane->submission_id(3);
  Description: Get/Set for submission internal id
  Returntype : integer

=cut

sub submission_id {
    my ($self,$submission_id) = @_;
    if (defined $submission_id and $submission_id != $self->{'submission_id'}){
	$self->{'submission_id'} = $submission_id;
	$self->dirty(1);
    }
    return $self->{'submission_id'};
}


=head2 submission

  Arg [1]    : submission name (optional)
  Example    : my $submission = $lane->submission();
               $lane->submission('g1k-sc-20080812-2');
  Description: Get/Set for sample submission.  Lazy-loads submission object from $self->submission_id.  If a submission name is supplied, then submission_id is set to the corresponding submission in the database.  If no such submission exists, returns undef.  Use add_submission to add a submission in this case.
  Returntype : VRTrack::Submission object

=cut

sub submission {
    my ($self,$submission) = @_;
    if ($submission){
        # get existing submission by name
        my $obj = $self->get_submission_by_name($submission);
        if ($obj){
            $self->{'submission'} = $obj;
            $self->{'submission_id'} = $obj->id;
            $self->dirty(1);
        }
        else {
            # warn "No such submission in the database";
            return undef; # explicitly return nothing.
        }
    }
    elsif ($self->{'submission'}){
        # already got a submission object.  We'll return it at the end.
    }
    else {  # lazy-load submission from database
        if ($self->submission_id){
            my $obj = VRTrack::Submission->new($self->{vrtrack},$self->submission_id);
            $self->{'submission'} = $obj;
        }
    }
    return $self->{'submission'};
}


=head2 add_submission

  Arg [1]    : submission name
  Example    : my $sub = $lane->add_submission('NA19820');
  Description: create a new submission, and if successful, return the object
  Returntype : VRTrack::Library object

=cut

sub add_submission {
    my ($self, $name) = @_;

    my $obj = $self->get_submission_by_name($name);
    if ($obj){
        warn "Submission $name is already present in the database\n";
        return undef;
    }
    else {
        $obj = VRTrack::Submission->create($self->{vrtrack}, $name);
        # populate caches
        $self->{'submission_id'} = $obj->id;
        $self->{'submission'} = $obj;
        $self->dirty(1);
    }
    return $self->{'submission'};
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
    my ($self,$auto_qc_status) = @_;
    if (defined $auto_qc_status and $auto_qc_status ne $self->{'auto_qc_status'}){
        my %allowed = map {$_ => 1} @{$self->list_enum_vals('lane','auto_qc_status')};
        unless ($allowed{lc($auto_qc_status)}){
            die "'$auto_qc_status' is not a defined auto_qc_status";
        }
	$self->{'auto_qc_status'} = $auto_qc_status;
	$self->dirty(1);
    }
    return $self->{'auto_qc_status'};
}


=head2 qc_status

  Arg [1]    : qc_status (optional)
  Example    : my $qc_status = $lane->qc_status();
	       $lane->qc_status('passed');
  Description: Get/Set for lane qc_status
  Returntype : string

=cut

sub qc_status {
    my ($self,$qc_status) = @_;
    if (defined $qc_status and $qc_status ne $self->{'qc_status'}){
        my %allowed = map {$_ => 1} @{$self->list_enum_vals('lane','qc_status')};
        unless ($allowed{lc($qc_status)}){
            die "'$qc_status' is not a defined qc_status";
        }
	$self->{'qc_status'} = $qc_status;
	$self->dirty(1);
    }
    return $self->{'qc_status'};
}


=head2 npg_qc_status

  Arg [1]    : npg_qc_status (optional)
  Example    : my $npg_qc_status = $lane->npg_qc_status();
	       $lane->npg_qc_status('pass');
  Description: Get/Set for lane npg_qc_status.  This is the manual QC that sequencescape/npg perform on a lane.
  Returntype : string

=cut

sub npg_qc_status {
    my ($self,$npg_qc_status) = @_;
    if (defined $npg_qc_status and $npg_qc_status ne $self->{'npg_qc_status'}){
        my %allowed = map {$_ => 1} @{$self->list_enum_vals('lane','npg_qc_status')};
        unless ($allowed{lc($npg_qc_status)}){
            die "'$npg_qc_status' is not a defined npg_qc_status";
        }
	$self->{'npg_qc_status'} = $npg_qc_status;
	$self->dirty(1);
    }
    return $self->{'npg_qc_status'};
}


=head2 genotype_status

  Arg [1]    : genotype_status (optional)
  Example    : my $genotype_status = $lane->genotype_status();
	       $lane->genotype_status('confirmed');
  Description: Get/Set for lane genotype check status
  Returntype : string

=cut

sub genotype_status {
    my ($self,$genotype_status) = @_;
    if (defined $genotype_status and $genotype_status ne $self->{'genotype_status'}){
        my %allowed = map {$_ => 1} @{$self->list_enum_vals('lane','gt_status')};
        unless ($allowed{lc($genotype_status)}){
            die "'$genotype_status' is not a defined genotype_status";
        }
	$self->{'genotype_status'} = $genotype_status;
	$self->dirty(1);
    }
    return $self->{'genotype_status'};
}


=head2 run_date

  Arg [1]    : run_date (optional)
  Example    : my $run_date = $lane->run_date();
               $lane->run_date('20080810123000');
  Description: Get/Set for lane run_date
  Returntype : string

=cut

sub run_date {
    my ($self,$run_date) = @_;
    if (defined $run_date and $run_date ne $self->{'run_date'}){
	$self->{'run_date'} = $run_date;
	$self->dirty(1);
    }
    return $self->{'run_date'};
}


=head2 changed

  Arg [1]    : changed (optional)
  Example    : my $changed = $lane->changed();
               $lane->changed('20080810123000');
  Description: Get/Set for lane changed
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

=head2 processed

  Arg [1]    : processed (optional)
  Description: Don't use this method, use is_processed instead.
  Returntype : string

=cut

sub processed {
    my ($self,$processed) = @_;
    if (defined $processed and $processed ne $self->{'processed'}){
	$self->{'processed'} = $processed;
	$self->dirty(1);
    }
    return $self->{'processed'};
}


=head2 is_processed

  Arg [1]    : flag, one of the flags listed in Core_obj::allowed_processed_flags();
  Arg [2]    : processed: 0 or 1 (optional)
  Example    : my $processed = $lane->is_processed('qc');
               $lane->processed('qc',1);
  Description: Get/Set for lane processed
  Returntype : 1 or 0

=cut

sub is_processed {
    my ($self,$flag,$processed) = @_;

    my %flags = VRTrack::Core_obj->allowed_processed_flags();
    if ( !exists($flags{$flag}) ) { croak qq[The flag "$flag" not recognised.\n]; }

    $flag = $flags{$flag};
    if ( defined $processed )
    {
        $processed = $processed ? $self->{processed}|$flag : $self->{processed}&(~$flag);
        if ( $processed != $self->{'processed'} )
        {
            $self->{'processed'} = $processed;
            $self->dirty(1);
        }
    }
    return $self->{'processed'} & $flag ? 1 : 0;
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
    my ($self) = @_;
    unless ($self->{'mappings'}){
	my @mappings;
    	foreach my $id (@{$self->mapping_ids()}){
	    my $obj = VRTrack::Mapstats->new($self->{vrtrack},$id);
	    push @mappings, $obj;
	}
	$self->{'mappings'} = \@mappings;
    }

    return $self->{'mappings'};
}


=head2 mapping_ids

  Arg [1]    : None
  Example    : my $mapping_ids = $lane->mapping_ids();
  Description: Returns a ref to an array of the ids of the mappings that are associated with this lane
  Returntype : ref to array of mapping ids

=cut

sub mapping_ids {
    my ($self) = @_;
    unless ($self->{'mapping_ids'}){
	my $sql = qq[select distinct(mapstats_id) from mapstats where lane_id=? and latest=true];
	my @files;
	my $sth = $self->{_dbh}->prepare($sql);

	if ($sth->execute($self->id)){
	    foreach(@{$sth->fetchall_arrayref()}){
		push @files, $_->[0];
	    }
	}
	else{
	    die(sprintf('Cannot retrieve files: %s', $DBI::errstr));
	}

	$self->{'mapping_ids'} = \@files;
    }
 
    return $self->{'mapping_ids'};
}


=head2 files

  Arg [1]    : None
  Example    : my $files = $lane->files();
  Description: Returns a ref to an array of the file objects that are associated with this lane.
  Returntype : ref to array of VRTrack::File objects

=cut

sub files {
    my ($self) = @_;
    unless ($self->{'files'}){
	my @files;
    	foreach my $id (@{$self->file_ids()}){
	    my $obj = VRTrack::File->new($self->{vrtrack},$id);
	    push @files, $obj;
	}
	$self->{'files'} = \@files;
    }

    return $self->{'files'};
}


=head2 file_ids

  Arg [1]    : None
  Example    : my $file_ids = $lane->file_ids();
  Description: Returns a ref to an array of the file ids that are associated with this lane
  Returntype : ref to array of file ids

=cut

sub file_ids {
    my ($self) = @_;
    unless ($self->{'file_ids'}){
	my $sql = qq[select file_id from file where lane_id=? and latest = true];
	my @files;
	my $sth = $self->{_dbh}->prepare($sql);

	if ($sth->execute($self->id)){
	    foreach(@{$sth->fetchall_arrayref()}){
		push @files, $_->[0];
	    }
	}
	else{
	    die(sprintf('Cannot retrieve files: %s', $DBI::errstr));
	}

	$self->{'file_ids'} = \@files;
    }
 
    return $self->{'file_ids'};
}


=head2 add_file

  Arg [1]    : file name
  Example    : my $newfile = $lib->add_file('1_s_1.fastq');
  Description: create a new file, and if successful, return the object
  Returntype : VRTrack::File object

=cut

sub add_file {
    my ($self, $name) = @_;
    $name or die "Must call add_file with file name";

    # File names should not be added twice - this can't be caught by the
    # database, as we expect there will be multiple rows for the same file.
    my $obj = VRTrack::File->new_by_name($self->{vrtrack},$name);
    if ($obj){
        warn "File $name is already present in the database\n";
        return undef;
    }
    $obj = VRTrack::File->create($self->{vrtrack}, $name);
    if ($obj){
        $obj->lane_id($self->id);
        $obj->update;
    }
    # clear caches
    delete $self->{'file_ids'};
    delete $self->{'files'};

    return $obj;

}


=head2 add_mapping

  Arg [1]    : none
  Example    : my $newmapping = $lib->add_mapping();
  Description: create a new mapping, and if successful, return the object
  Returntype : VRTrack::Mapstats object

=cut

sub add_mapping {
    my ($self) = @_;
    my $obj = VRTrack::Mapstats->create($self->{vrtrack});
    if ($obj){
        $obj->lane_id($self->id);
        $obj->update;
    }
    # clear caches
    delete $self->{'mapping_ids'};
    delete $self->{'mappings'};
    return $obj;
}


=head2 get_file_by_name

  Arg [1]    : file name
  Example    : my $file = $lane->get_file_by_name('My file');
  Description: retrieve file object on this lane by name
  Returntype : VRTrack::File object

=cut

sub get_file_by_name {
    my ($self, $name) = @_;
    #my $obj = VRTrack::File->new_by_name($self->{vrtrack},$name);
    my @match = grep {$_->name eq $name} @{$self->files};
    if (scalar @match > 1){ # shouldn't happen
        die "More than one matching file with name $name";
    }
    my $obj;
    if (@match){
        $obj = $match[0];
    }

    return $obj;
}


=head2 get_file_by_id

  Arg [1]    : internal file id
  Example    : my $file = $lib->get_file_by_id(47);
  Description: retrieve file object by internal db file id
  Returntype : VRTrack::File object

=cut

sub get_file_by_id {
    my ($self, $id) = @_;
    my $obj = VRTrack::File->new($self->{vrtrack},$id);
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
    foreach (@{$self->files}){
        push @desc, $_;
    }
    foreach (@{$self->mappings}){
        push @desc, $_;
        push @desc, @{$_->descendants};
    }
    return \@desc;
}

1;
