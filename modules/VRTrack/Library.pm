package VRTrack::Library; 
# author: jws
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

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;
use Carp;
no warnings 'uninitialized';
use VRTrack::Lane;
use VRTrack::Library_type;
use VRTrack::Seq_centre;
use VRTrack::Seq_tech;
use VRTrack::Core_obj;
our @ISA = qw(VRTrack::Core_obj);

=head2 fields_dispatch

  Arg [1]    : none
  Example    : my $fieldsref = $lib->fields_dispatch();
  Description: Returns hashref dispatch table keyed on database field
               Used internally for new and update methods
  Returntype : hashref

=cut

sub fields_dispatch {
    my $self = shift;
    my %fields = ( 
                'library_id'        => sub { $self->id(@_)},
                'sample_id'         => sub { $self->sample_id(@_)},
                'ssid'              => sub { $self->ssid(@_)},
                'name'              => sub { $self->name(@_)},
                'hierarchy_name'    => sub { $self->hierarchy_name(@_)},
                'prep_status'       => sub { $self->prep_status(@_)},
                'qc_status'         => sub { $self->qc_status(@_)},
                'insert_size'       => sub { $self->insert_size(@_)},
                'library_type_id'   => sub { $self->library_type_id(@_)},
                'seq_centre_id'     => sub { $self->seq_centre_id(@_)},
                'seq_tech_id'       => sub { $self->seq_tech_id(@_)},
                'open'              => sub { $self->open(@_)},
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

  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : library name
  Example    : my $library = VRTrack::Library->new_by_name($vrtrack, $name)
  Description: Class method. Returns latest Library object by name.  If no such name is in the database, returns undef.  Dies if multiple names match.
  Returntype : VRTrack::Library object

=cut

sub new_by_name {
    my ($class,$vrtrack, $name) = @_;
    die "Need to call with a vrtrack handle, name" unless ($vrtrack && $name);
    return $class->new_by_field_value($vrtrack, 'name',$name);
}


=head2 new_by_hierarchy_name

  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : library hierarchy_name
  Example    : my $library = VRTrack::Library->new_by_hierarchy_name($vrtrack, $hierarchy_name)
  Description: Class method. Returns latest Library object by hierarchy_name.  If no such hierarchy_name is in the database, returns undef.  Dies if multiple hierarchy_names match.
  Returntype : VRTrack::Library object

=cut

sub new_by_hierarchy_name {
    my ($class,$vrtrack, $hierarchy_name) = @_;
    die "Need to call with a vrtrack handle, hierarchy_name" unless ($vrtrack && $hierarchy_name);
    return $class->new_by_field_value($vrtrack, 'hierarchy_name',$hierarchy_name);
}


=head2 new_by_ssid

  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : library sequencescape id
  Example    : my $library = VRTrack::Library->new_by_ssid($vrtrack, $ssid);
  Description: Class method. Returns latest Library object by ssid.  If no such ssid is in the database, returns undef
  Returntype : VRTrack::Library object

=cut

sub new_by_ssid {
    my ($class,$vrtrack, $ssid) = @_;
    die "Need to call with a vrtrack handle, ssid" unless ($vrtrack && $ssid);
    return $class->new_by_field_value($vrtrack, 'ssid',$ssid);
}


=head2 create

  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : library name
  Example    : my $library = VRTrack::Library->create($vrtrack, $name)
  Description: Class method.  Creates new Library object in the database.
               Dies if the name or hierarchy_name already exists.
  Returntype : VRTrack::Library object

=cut

sub create {
    my ($class,$vrtrack, $name) = @_;
    die "Need to call with a vrtrack handle and name" unless ($vrtrack && $name);
    if ( $vrtrack->isa('DBI::db') ) { croak "The interface has changed, expected vrtrack reference.\n"; }
    my $dbh = $vrtrack->{_dbh};

    my $hierarchy_name = $name;
    $hierarchy_name =~ s/\W+/_/g;

    # prevent adding a library with an existing name
    if ($class->is_name_in_database($vrtrack, $name, $hierarchy_name)){
        die "Already a library by name $name/$hierarchy_name";
    }

    $dbh->do (qq[LOCK TABLE library WRITE]);
    my $sql = qq[select max(library_id) as id from library];
    my $sth = $dbh->prepare($sql);
    my $next_id;
    if ($sth->execute()){
	my $data = $sth->fetchrow_hashref;
	unless ($data){
            $dbh->do (qq[UNLOCK TABLES]);
            die( sprintf("Can't retrieve next library id: %s", $DBI::errstr));
	}
        $next_id = $data->{'id'};
        $next_id++;
    }
    else{
	die(sprintf("Can't retrieve next library id: %s", $DBI::errstr));
    }

    $sql = qq[INSERT INTO library (library_id, name, hierarchy_name, changed, latest) 
                 VALUES (?,?,?,now(),true)];

    $sth = $dbh->prepare($sql);
    unless ($sth->execute( $next_id, $name, $hierarchy_name )) {
        $dbh->do (qq[UNLOCK TABLES]);
        die( sprintf('DB load insert failed: %s %s', $next_id, $DBI::errstr));
    }

    $dbh->do (qq[UNLOCK TABLES]);

    return $class->new($vrtrack, $next_id);
}


=head2 is_name_in_database

  Arg [1]    : library name
  Arg [2]    : hierarchy name
  Example    : if(VRTrack::Library->is_name_in_database($vrtrack, $name, $hname)
  Description: Class method. Checks to see if a name or hierarchy name is already used in the library table.
  Returntype : boolean

=cut

sub is_name_in_database {
    my ($class, $vrtrack, $name, $hname) = @_;
    die "Need to call with a vrtrack handle, name, hierarchy name" unless ($vrtrack && $name && $hname);
    my $dbh = $vrtrack->{_dbh};
    my $sql = qq[select library_id from library where latest=true and (name = ? or hierarchy_name = ?) ];
    my $sth = $dbh->prepare($sql);

    my $already_used = 0;
    if ($sth->execute($name,$hname)){
        my $data = $sth->fetchrow_hashref;
        if ($data){
            $already_used = 1;
        }
    }
    else{
        die(sprintf('Cannot retrieve library by $name: %s', $DBI::errstr));
    }
    return $already_used;
}



###############################################################################
# Object methods
###############################################################################

=head2 lanes

  Arg [1]    : None
  Example    : my $lanes = $library->lanes();
  Description: Returns a ref to an array of the lane objects that are associated with this library.
  Returntype : ref to array of VRTrack::Lane objects

=cut

sub lanes {
    my ($self) = @_;
    unless ($self->{'lanes'}){
        my @lanes;
        foreach my $id (@{$self->lane_ids()}){
            my $obj = VRTrack::Lane->new($self->{vrtrack},$id);
            push @lanes, $obj if $obj;
        }
        $self->{'lanes'} = \@lanes;
    }

    return $self->{'lanes'};
}


=head2 lane_ids

  Arg [1]    : None
  Example    : my $lane_ids = $library->lane_ids();
  Description: Returns a ref to an array of the lane IDs that are associated with this library
  Returntype : ref to array of integer lane IDs

=cut

sub lane_ids {
    my ($self) = @_;
    unless ($self->{'lane_ids'}){
        my $sql = qq[select lane_id from lane where library_id=? and latest=true];
        my @lanes;
        my $sth = $self->{_dbh}->prepare($sql);

        if ($sth->execute($self->id)){
            foreach(@{$sth->fetchall_arrayref()}){
                push @lanes, $_->[0];
            }
        }
        else{
            die(sprintf('Cannot retrieve lanes: %s', $DBI::errstr));
        }

        $self->{'lane_ids'} = \@lanes;
    }
 
    return $self->{'lane_ids'};
}


=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $lib->id();
               $lib->id(104);
  Description: Get/Set for internal db ID of a library
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


=head2 sample_id

  Arg [1]    : sample_id (optional)
  Example    : my $sample_id = $lib->sample_id();
               $lib->sample_id(104);
  Description: Get/Set for ID of a library
  Returntype : Internal ID integer

=cut

sub sample_id {
    my ($self,$sample_id) = @_;
    if (defined $sample_id and $sample_id != $self->{'sample_id'}){
        $self->{'sample_id'} = $sample_id;
        $self->dirty(1);
    }
    return $self->{'sample_id'};
}


=head2 ssid

  Arg [1]    : ssid (optional)
  Example    : my $ssid = $lib->ssid();
               $lib->ssid('104');
  Description: Get/Set for SequenceScape ID of a library
  Returntype : integer

=cut

sub ssid {
    my ($self,$ssid) = @_;
    if (defined $ssid and $ssid != $self->{'ssid'}){
        $self->{'ssid'} = $ssid;
        $self->dirty(1);
    }
    return $self->{'ssid'};
}


=head2 hierarchy_name

  Arg [1]    : directory name (optional)
  Example    : my $hname = $library->hierarchy_name();
  Description: Get/set library hierarchy name.  This is the directory name (without path) that the library will be named in a file hierarchy.
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
  Example    : my $name = $lib->name();
               $lib->name('104');
  Description: Get/Set for library name
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


=head2 prep_status

  Arg [1]    : prep_status (optional)
  Example    : my $prep_status = $lib->prep_status();
               $lib->prep_status('104');
  Description: Get/Set for library prep_status
  Returntype : string

=cut

sub prep_status {
    my ($self,$prep_status) = @_;
    if (defined $prep_status and $prep_status ne $self->{'prep_status'}){
        my %allowed = map {$_ => 1} @{$self->list_enum_vals('library','prep_status')};
        unless ($allowed{lc($prep_status)}){
            die "'$prep_status' is not a defined qc_status";
        }
        $self->{'prep_status'} = $prep_status;
        $self->dirty(1);
    }
    return $self->{'prep_status'};
}


=head2 qc_status

  Arg [1]    : qc_status (optional)
  Example    : my $qc_status = $lib->qc_status();
               $lib->qc_status('104');
  Description: Get/Set for library qc_status.  Checks database enum for allowed values before allowing qc_status to be set.
  Returntype : string

=cut

sub qc_status {
    my ($self,$qc_status) = @_;
    if (defined $qc_status and $qc_status ne $self->{'qc_status'}){
        my %allowed = map {$_ => 1} @{$self->list_enum_vals('library','qc_status')};
        unless ($allowed{lc($qc_status)}){
            die "'$qc_status' is not a defined qc_status";
        }

        $self->{'qc_status'} = $qc_status;
        $self->dirty(1);
    }
    return $self->{'qc_status'};
}


=head2 insert_size

  Arg [1]    : insert_size (optional)
  Example    : my $insert_size = $lib->insert_size();
               $lib->insert_size(76);
  Description: Get/Set for library insert_size
  Returntype : integer

=cut

sub insert_size {
    my ($self,$insert_size) = @_;
    if (defined $insert_size and $insert_size != $self->{'insert_size'}){
        $self->{'insert_size'} = $insert_size;
        $self->dirty(1);
    }
    return $self->{'insert_size'};
}


=head2 library_type_id

  Arg [1]    : library_type_id (optional)
  Example    : my $library_type_id = $lib->library_type_id();
               $lib->library_type_id(1);
  Description: Get/Set for library library_type_id
  Return_type_id : string

=cut

sub library_type_id {
    my ($self,$library_type_id) = @_;
    if (defined $library_type_id and $library_type_id != $self->{'library_type_id'}){
        $self->{'library_type_id'} = $library_type_id;
        $self->dirty(1);
    }
    return $self->{'library_type_id'};
}


=head2 library_type

  Arg [1]    : library_type name (optional)
  Example    : my $library_type = $samp->library_type();
               $samp->library_type('DSS');
  Description: Get/Set for sample library_type.  Lazy-loads library_type object from $self->library_type_id.  If a library_type name is supplied, then library_type_id is set to the corresponding library_type in the database.  If no such library_type exists, returns undef.  Use add_library_type to add a library_type in this case.
  Returntype : VRTrack::Library_type object

=cut

sub library_type {
    my ($self,$name) = @_;
    if ($name){
        # get existing library_type by name
        my $obj = $self->get_library_type_by_name($name);
        if ($obj){
            # Have we actually changed?
            if ($self->library_type_id != $obj->id){
                $self->library_type_id($obj->id);
                $self->dirty(1);
            }
            $self->{'library_type'} = $obj;
        }
        else {
            # warn "No such library_type in the database";
            return undef; # explicitly return nothing.
        }
    }
    elsif ($self->{'library_type'}){
        # already got a library_type object.  We'll return it at the end.
    }
    else {  # lazy-load library_type from database
        if ($self->library_type_id){
            my $obj = VRTrack::Library_type->new($self->{vrtrack},$self->library_type_id);
            $self->{'library_type'} = $obj;
        }
    }
    return $self->{'library_type'};
}


=head2 add_library_type

  Arg [1]    : library_type name
  Example    : my $library_type = $samp->add_library_type('DSS');
  Description: create a new library_type, and if successful, return the object
  Returntype : VRTrack::Library object

=cut

sub add_library_type {
    my ($self, $name) = @_;

    my $obj = $self->get_library_type_by_name($name);
    if ($obj){
        warn "Library_type $name is already present in the database\n";
        return undef;
    }
    else {
        $obj = VRTrack::Library_type->create($self->{vrtrack}, $name);
        # populate caches
        $self->{'library_type_id'} = $obj->id;
        $self->{'library_type'} = $obj;
        $self->dirty(1);
    }
    return $self->{'library_type'};
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


=head2 seq_centre_id

  Arg [1]    : seq_centre_id (optional)
  Example    : my $seq_centre_id = $lib->seq_centre_id();
               $lib->seq_centre_id(1);
  Description: Get/Set for library sequencing seq_centre_id
  Returntype : string

=cut

sub seq_centre_id {
    my ($self,$seq_centre_id) = @_;
    if (defined $seq_centre_id and $seq_centre_id != $self->{'seq_centre_id'}){
        $self->{'seq_centre_id'} = $seq_centre_id;
        $self->dirty(1);
    }
    return $self->{'seq_centre_id'};
}


=head2 seq_centre

  Arg [1]    : seq_centre name (optional)
  Example    : my $seq_centre = $samp->seq_centre();
               $samp->seq_centre('SC');
  Description: Get/Set for sample seq_centre.  Lazy-loads seq_centre object from $self->seq_centre_id.  If a seq_centre name is supplied, then seq_centre_id is set to the corresponding seq_centre in the database.  If no such seq_centre exists, returns undef.  Use add_seq_centre to add a seq_centre in this case.
  Returntype : VRTrack::Seq_centre object

=cut

sub seq_centre {
    my ($self,$name) = @_;
    if ($name){
        # get existing seq_centre by name
        my $obj = $self->get_seq_centre_by_name($name);
        if ($obj){
            # Have we actually changed?
            if ($self->seq_centre_id != $obj->id){
                $self->seq_centre_id($obj->id);
                $self->dirty(1);
            }
            $self->{'seq_centre'} = $obj;
        }
        else {
            # warn "No such seq_centre in the database";
            return undef; # explicitly return nothing.
        }
    }
    elsif ($self->{'seq_centre'}){
        # already got a seq_centre object.  We'll return it at the end.
    }
    else {  # lazy-load seq_centre from database
        if ($self->seq_centre_id){
            my $obj = VRTrack::Seq_centre->new($self->{vrtrack},$self->seq_centre_id);
            $self->{'seq_centre'} = $obj;
        }
    }
    return $self->{'seq_centre'};
}


=head2 add_seq_centre

  Arg [1]    : seq_centre name
  Example    : my $seq_centre = $samp->add_seq_centre('SC');
  Description: create a new seq_centre, and if successful, return the object
  Returntype : VRTrack::Library object

=cut

sub add_seq_centre {
    my ($self, $name) = @_;

    my $obj = $self->get_seq_centre_by_name($name);
    if ($obj){
        warn "Seq_centre $name is already present in the database\n";
        return undef;
    }
    else {
        $obj = VRTrack::Seq_centre->create($self->{vrtrack}, $name);
        # populate caches
        $self->{'seq_centre_id'} = $obj->id;
        $self->{'seq_centre'} = $obj;
        $self->dirty(1);
    }
    return $self->{'seq_centre'};
}


=head2 get_seq_centre_by_name

  Arg [1]    : seq_centre_name
  Example    : my $seq_centre = $samp->get_seq_centre_by_name('SC');
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
    my ($self,$seq_tech_id) = @_;
    if (defined $seq_tech_id and $seq_tech_id != $self->{'seq_tech_id'}){
        $self->{'seq_tech_id'} = $seq_tech_id;
        $self->dirty(1);
    }
    return $self->{'seq_tech_id'};
}


=head2 seq_tech

  Arg [1]    : seq_tech name (optional)
  Example    : my $seq_tech = $samp->seq_tech();
               $samp->seq_tech('SLX');
  Description: Get/Set for sample seq_tech.  Lazy-loads seq_tech object from $self->seq_tech_id.  If a seq_tech name is supplied, then seq_tech_id is set to the corresponding seq_tech in the database.  If no such seq_tech exists, returns undef.  Use add_seq_tech to add a seq_tech in this case.
  Returntype : VRTrack::Seq_tech object

=cut

sub seq_tech {
    my ($self,$name) = @_;
    if ($name){
        # get existing seq_tech by name
        my $obj = $self->get_seq_tech_by_name($name);
        if ($obj){
            # Have we actually changed?
            if ($self->seq_tech_id != $obj->id){
                $self->seq_tech_id($obj->id);
                $self->dirty(1);
            }
            $self->{'seq_tech'} = $obj;
        }
        else {
            # warn "No such seq_tech in the database";
            return undef; # explicitly return nothing.
        }
    }
    elsif ($self->{'seq_tech'}){
        # already got a seq_tech object.  We'll return it at the end.
    }
    else {  # lazy-load seq_tech from database
        if ($self->seq_tech_id){
            my $obj = VRTrack::Seq_tech->new($self->{vrtrack},$self->seq_tech_id);
            $self->{'seq_tech'} = $obj;
        }
    }
    return $self->{'seq_tech'};
}


=head2 add_seq_tech

  Arg [1]    : seq_tech name
  Example    : my $seq_tech = $samp->add_seq_tech('SLX');
  Description: create a new seq_tech, and if successful, return the object
  Returntype : VRTrack::Library object

=cut

sub add_seq_tech {
    my ($self, $name) = @_;

    my $obj = $self->get_seq_tech_by_name($name);
    if ($obj){
        warn "Seq_tech $name is already present in the database\n";
        return undef;
    }
    else {
        $obj = VRTrack::Seq_tech->create($self->{vrtrack}, $name);
        # populate caches
        $self->{'seq_tech_id'} = $obj->id;
        $self->{'seq_tech'} = $obj;
        $self->dirty(1);
    }
    return $self->{'seq_tech'};
}


=head2 get_seq_tech_by_name

  Arg [1]    : seq_tech_name
  Example    : my $seq_tech = $samp->get_seq_tech_by_name('SLX');
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
    my ($self,$open) = @_;
    if (defined $open and $open != $self->{'open'}){
        $self->{'open'} = $open ? 1 : 0;
        $self->dirty(1);
    }
    return $self->{'open'};
}


=head2 changed

  Arg [1]    : timestamp (optional)
  Example    : my $changed = $lib->changed();
               $lib->changed('20080810123000');
  Description: Get/Set for library changed
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
  Example    : $lib->dirty(1);
  Description: Get/Set for library properties having been altered.
  Returntype : boolean

=cut

sub dirty {
    my ($self,$dirty) = @_;
    if (defined $dirty){
        $self->{_dirty} = $dirty ? 1 : 0;
    }
    return $self->{_dirty};
}


=head2 add_lane

  Arg [1]    : lane name
  Example    : my $newlane = $lib->add_lane('2631_3');
  Description: create a new lane, and if successful, return the object
  Returntype : VRTrack::Lane object

=cut

sub add_lane {
    my ($self, $name) = @_;
    $name or die "Must call add_lane with lane name";
    # Lane names should not be added twice - this can't be caught by the
    # database, as we expect there will be multiple rows for the same lane.
    my $obj = VRTrack::Lane->new_by_name($self->{vrtrack},$name);
    if ($obj){
        warn "Lane $name is already present in the database\n";
        return undef;
    }
    $obj = VRTrack::Lane->create($self->{vrtrack}, $name);
    if ($obj){
        $obj->library_id($self->id);
        $obj->update;
    }
    # clear caches
    delete $self->{'lane_ids'};
    delete $self->{'lanes'};

    return $obj;

}


=head2 get_lane_by_id

  Arg [1]    : internal lane id
  Example    : my $lane = $lib->get_lane_by_id(47);
  Description: retrieve lane object by internal db lane id
  Returntype : VRTrack::Lane object

=cut

sub get_lane_by_id {
    my ($self, $id) = @_;
    my @match = grep {$_->id == $id} @{$self->lanes};
    if (scalar @match > 1){ # shouldn't happen
        die "More than one matching lane with id $id";
    }
    my $obj;
    if (@match){
        $obj = $match[0];
    }
    return $obj;
}


=head2 get_lane_by_name

  Arg [1]    : lane name
  Example    : my $lane = $track->get_lane_by_name('My lane');
  Description: retrieve lane object attached to this Library, by name
  Returntype : VRTrack::Lane object

=cut

sub get_lane_by_name {
    my ($self, $name) = @_;
    #my $obj = VRTrack::Lane->new_by_name($self->{vrtrack},$name);
    my @match = grep {$_->name eq $name} @{$self->lanes};
    if (scalar @match > 1){ # shouldn't happen
        die "More than one matching lane with name $name";
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
    # Requests aren't linked into the API yet
    #foreach (@{$self->requests}){
    #    push @desc, $_;
    #}
    foreach (@{$self->lanes}){
        push @desc, $_;
        push @desc, @{$_->descendants};
    }
    return \@desc;
}
1;
