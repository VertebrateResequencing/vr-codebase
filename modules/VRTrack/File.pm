package VRTrack::File; 
# author: jws
=head1 NAME

VRTrack::File - Sequence Tracking File object

=head1 SYNOPSIS
    my $file= VRTrack::File->new($dbh, $file_id);

    my $id = $file->id();
    my $qc_status = $file->qc_status();

=head1 DESCRIPTION

An object describing the tracked properties of a file.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;
no warnings 'uninitialized';
use VRTrack::Core_obj;
our @ISA = qw(VRTrack::Core_obj);

=head2 fields_dispatch

  Arg [1]    : none
  Example    : my $fieldsref = $file->fields_dispatch();
  Description: Returns hashref dispatch table keyed on database field
               Used internally for new and update methods
  Returntype : hashref

=cut

sub fields_dispatch {
    my $self = shift;
    my %fields = ( 
                'file_id'           => sub { $self->id(@_)},
                'lane_id'           => sub { $self->lane_id(@_)},
                'name'              => sub { $self->name(@_)},
                'hierarchy_name'    => sub { $self->hierarchy_name(@_)},
                'qc_status'         => sub { $self->qc_status(@_)},
                'type'              => sub { $self->type(@_)},
                'readlen'           => sub { $self->read_len(@_)},
                'raw_reads'         => sub { $self->raw_reads(@_)},
                'raw_bases'         => sub { $self->raw_bases(@_)},
                'recalibrated'      => sub { $self->is_recalibrated(@_)},
                'md5'               => sub { $self->md5(@_)},
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
  Arg [2]    : file name
  Example    : my $file = VRTrack::File->new_by_name($dbh, $name)
  Description: Class method. Returns latest File object by name.  If no such name is in the database, returns undef.  Dies if multiple names match.
  Returntype : VRTrack::File object

=cut

sub new_by_name {
    my ($class,$dbh, $name) = @_;
    die "Need to call with a db handle, name" unless ($dbh && $name);
    return $class->new_by_field_value($dbh, 'name',$name);
}


=head2 new_by_hierarchy_name

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : file hierarchy_name
  Example    : my $file = VRTrack::File->new_by_hierarchy_name($dbh, $hierarchy_name)
  Description: Class method. Returns latest File object by hierarchy_name.  If no such hierarchy_name is in the database, returns undef.  Dies if multiple hierarchy_names match.
  Returntype : VRTrack::File object

=cut

sub new_by_hierarchy_name {
    my ($class,$dbh, $hierarchy_name) = @_;
    die "Need to call with a db handle, hierarchy_name" unless ($dbh && $hierarchy_name);
    return $class->new_by_field_value($dbh, 'hierarchy_name',$hierarchy_name);
}


=head2 create

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : file name
  Example    : my $file = VRTrack::File->create($dbh, $name)
  Description: Class method.  Creates new File object in the database.
  Returntype : VRTrack::File object

=cut

sub create {
    my ($class,$dbh, $name) = @_;
    die "Need to call with a db handle and name" unless ($dbh && $name);

    # prevent adding a file with an existing name
    if ($class->is_name_in_database($dbh, $name, $name)){
        die "Already a file by name $name";
    }

    $dbh->do (qq[LOCK TABLE file WRITE]);
    my $sql = qq[select max(file_id) as id from file];
    my $sth = $dbh->prepare($sql);
    my $next_id;
    if ($sth->execute()){
	my $data = $sth->fetchrow_hashref;
	unless ($data){
            $dbh->do (qq[UNLOCK TABLES]);
            die( sprintf("Can't retrieve next file id: %s", $DBI::errstr));
	}
        $next_id = $data->{'id'};
        $next_id++;
    }
    else{
	die(sprintf("Can't retrieve next file id: %s", $DBI::errstr));
    }

    $sql = qq[INSERT INTO file (file_id, name, changed, latest) 
                 VALUES (?,?,now(),true)];

    $sth = $dbh->prepare($sql);
    unless ($sth->execute( $next_id, $name )) {
        $dbh->do (qq[UNLOCK TABLES]);
        die( sprintf('DB load insert failed: %s %s', $next_id, $DBI::errstr));
    }

    $dbh->do (qq[UNLOCK TABLES]);

    return $class->new($dbh, $next_id);
}


=head2 is_name_in_database

  Arg [1]    : file name
  Arg [2]    : hierarchy name
  Example    : if(VRTrack::File->is_name_in_database($dbh, $name,$hname)
  Description: Class method. Checks to see if a name or hierarchy name is already used in the file table.
  Returntype : boolean

=cut

sub is_name_in_database {
    my ($class, $dbh, $name, $hname) = @_;
    die "Need to call with a db handle, name, hierarchy name" unless ($dbh && $name && $hname);
    my $sql = qq[select file_id from file where latest=true and (name = ? or hierarchy_name = ?) ];
    my $sth = $dbh->prepare($sql);

    my $already_used = 0;
    if ($sth->execute($name,$hname)){
        my $data = $sth->fetchrow_hashref;
        if ($data){
            $already_used = 1;
        }
    }
    else{
        die(sprintf('Cannot retrieve file by $name: %s', $DBI::errstr));
    }
    return $already_used;
}


###############################################################################
# Object methods
###############################################################################

=head2 dirty

  Arg [1]    : boolean for dirty status
  Example    : $file->dirty(1);
  Description: Get/Set for file properties having been altered.
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
  Example    : my $id = $lib->id();
               $lib->id(104);
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


=head2 is_recalibrated

  Arg [1]    : boolean for recalibrated status
  Example    : $file->is_recalibrated(1);
  Description: Get/Set for file having been recalibrated
  Returntype : boolean

=cut

sub is_recalibrated {
    my ($self,$recalibrated) = @_;
    if (defined $recalibrated){
	$self->{_recalibrated} = $recalibrated ? 1 : 0;
    }
    return $self->{_recalibrated};
}


=head2 hierarchy_name

  Arg [1]    : directory name (optional)
  Example    : my $hname = $file->hierarchy_name();
  Description: Get/set file hierarchy name.  This is the directory name (without path) that the file will be named in a file hierarchy.
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
  Example    : my $name = $file->name();
	       $file->name('1111_s_3_1.fastq');
  Description: Get/Set for name of a file
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


=head2 md5

  Arg [1]    : md5 (optional)
  Example    : my $md5 = $file->md5();
	       $file->md5('1027759a77eb562fcab253c9f01d7661');
  Description: Get/Set for md5 of a file
  Returntype : string

=cut

sub md5 {
    my ($self,$md5) = @_;
    if (defined $md5 and $md5 ne $self->{'md5'}){
	$self->{'md5'} = $md5;
	$self->dirty(1);
    }
    return $self->{'md5'};
}


=head2 lane_id

  Arg [1]    : lane_id (optional)
  Example    : my $lane_id = $file->lane_id();
	       $file->lane_id('104');
  Description: Get/Set for ID of a file
  Returntype : SequenceScape ID (usu. integer)

=cut

sub lane_id {
    my ($self,$lane_id) = @_;
    if (defined $lane_id and $lane_id != $self->{'lane_id'}){
	$self->{'lane_id'} = $lane_id;
	$self->dirty(1);
    }
    return $self->{'lane_id'};
}


=head2 qc_status

  Arg [1]    : qc_status (optional) ['pending','passed','failed']
  Example    : my $qc_status = $file->qc_status();
	       $file->qc_status('passed');

               NB before a file has been 'imported', i.e. copied into a disk
               hierarchy, the qc_status should be NULL.  On import, it should
               be set to 'no_qc'.
  Description: Get/Set for file qc_status
  Returntype : string

=cut

sub qc_status {
    my ($self,$qc_status) = @_;
    if (defined $qc_status and $qc_status ne $self->{'qc_status'}){
        my %allowed = map {$_ => 1} @{$self->list_enum_vals('file','qc_status')};
        unless ($allowed{lc($qc_status)}){
            die "'$qc_status' is not a defined qc_status";
        }
	$self->{'qc_status'} = $qc_status;
	$self->dirty(1);
    }
    return $self->{'qc_status'};
}


=head2 changed

  Arg [1]    : changed (optional)
  Example    : my $changed = $file->changed();
	       $file->changed('20090101102300');
  Description: Get/Set for file changed
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


=head2 type

  Arg [1]    : type (optional) [0,1,2]
  Example    : my $type = $file->type();
	       $file->type(1);
  Description: Get/Set for file type - 0 is single-end, 1 fwd, 2 rev
  Returntype : integer

=cut

sub type {
    my ($self,$type) = @_;
    if (defined $type and $type != $self->{'type'}){
	$self->{'type'} = $type;
	$self->dirty(1);
    }
    return $self->{'type'};
}


=head2 raw_reads

  Arg [1]    : number of sequences in file
  Example    : my $num_seqs = $file->raw_reads();
	       $file->raw_reads(1_000_000);
  Description: Get/Set for number of sequences in a fasta file
  Returntype : integer

=cut

sub raw_reads {
    my ($self,$num_seqs) = @_;
    if (defined $num_seqs and $num_seqs != $self->{'seq_count'}){
	$self->{'seq_count'} = $num_seqs;
	$self->dirty(1);
    }
    return $self->{'seq_count'};
}



=head2 raw_bases

  Arg [1]    : number of basepairs in file
  Example    : my $num_bp = $file->raw_bases();
	       $file->raw_bases(1_000_000);
  Description: Get/Set for total number of unfiltered basepairs in file
  Returntype : integer

=cut

sub raw_bases {
    my ($self,$num_bp) = @_;
    if (defined $num_bp and $num_bp != $self->{'raw_bases'}){
	$self->{'raw_bases'} = $num_bp;
	$self->dirty(1);
    }
    return $self->{'raw_bases'};
}


=head2 read_len

  Arg [1]    : read length
  Example    : my $readlen = $file->read_len();
	       $file->read_len(54);
  Description: Get/Set for unclipped read length in file
  Returntype : integer

=cut

sub read_len {
    my ($self,$readlen) = @_;
    if (defined $readlen and $readlen != $self->{'readlen'}){
	$self->{'readlen'} = $readlen;
	$self->dirty(1);
    }
    return $self->{'readlen'};
}


1;
