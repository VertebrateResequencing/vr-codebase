package VRTrack::Core_obj; 
# author: jws
=head1 NAME

VRTrack::Core_obj - Sequence Tracking Core_obj object

=head1 SYNOPSIS

=head1 DESCRIPTION

This is the superclass of the core objects (Project, Sample, Library, etc) in
VRTrack, and implements some common methods.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;
use Carp;
no warnings 'uninitialized';


=head2 new

  Arg [1]    : vrtrack handle
  Arg [2]    : obj id
  Example    : my $obj= $class->new($vrtrack, $id)
  Description: Returns core objects by id
  Returntype : $class object

=cut

sub new {
    my ($class, $vrtrack, $id) = @_;
    confess "Need to call with a vrtrack reference and id" unless ($vrtrack && $id);
    if ( $vrtrack->isa('DBI::db') ) { croak "The interface has changed, expected vrtrack reference.\n"; }

    my $self = {};
    bless ($self, $class);
    $self->{_dbh} = $vrtrack->{_dbh};
    $self->{vrtrack} = $vrtrack;
    
    # database tables for core objects are named the same as the class,
    $class =~/VRTrack::(\w+)$/;
    my $table = lc($1);
    $table or die "Unrecognised classname $class\n";

    my $fieldsref = $self->fields_dispatch;
    my $sql = qq[select row_id,].(join ", ",keys %$fieldsref).qq[ from $table ];
    $sql .= qq[where ${table}_id = ? and latest = true];

    my $sth = $self->{_dbh}->prepare($sql);

    if ($sth->execute($id)){
        my $data = $sth->fetchrow_hashref;
	unless ($data){
	    return undef;
	}
        foreach (keys %$fieldsref){
            $fieldsref->{$_}->($data->{$_});
        }

        # handle row_id as a special case outside fields_dispatch, as otherwise
        # we'd need to remove it from fields_dispatch before update.
        # Note also it is hard-coded in the select statement above
        $self->row_id($data->{'row_id'});

	$self->dirty(0);    # unset the dirty flag
    }
    else{
	die(sprintf('Cannot retrieve $table: %s', $DBI::errstr));
    }

    return $self;
}


=head2 new_by_field_value

  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : field name
  Arg [3]    : field value
  Example    : my $obj = $class->new_by_field_value($vrtrack, 'name',$name)
  Description: Class method. Returns latest $class object by field name and value.  If no such value is in the database, returns undef.
               Dies if there is more than one matching record.
  Returntype : $class object

=cut

sub new_by_field_value {
    my ($class,$vrtrack, $field, $value) = @_;
    die "Need to call with a vrtrack handle, field name, field value" unless ($vrtrack && $field && defined $value);
    if ( $vrtrack->isa('DBI::db') ) { croak "The interface has changed, expected vrtrack reference.\n"; }
    my $dbh = $vrtrack->{_dbh};
    
    # database tables for core objects are named the same as the class,
    $class =~/VRTrack::(\w+)$/;
    my $table = lc($1);
    $table or die "Unrecognised classname $class\n";

    # check field exists
    my $colnames = $dbh->selectcol_arrayref(qq[select column_name from information_schema.columns where table_name='$table']);
    my %cols = map { $_ => 1 } @$colnames;
    unless (exists($cols{lc($field)})){
        die "No such column $field in $table table\n";
    }

    # retrieve lane_id
    my $sql = qq[select ${table}_id from $table where $field = ? and latest = true];
    my $sth = $dbh->prepare($sql);
    my $id;
    if ($sth->execute($value)){
        my $data = $sth->fetchall_arrayref({}); #return array of hashes
        unless (@$data){
            return undef;
        }
        if (scalar @$data > 1){
            die "$field = $value is not a unique identifier for $table\n";
        }
        $id = $data->[0]{"${table}_id"};
    }
    else{
        die(sprintf('Cannot retrieve $class by %s = %s: %s', ($field,$value,$DBI::errstr)));
    }
    return $class->new($vrtrack, $id);
}


=head2 create

  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : file name
  Example    : my $file = VRTrack::File->create($vrtrack, $name)
  Description: Class method.  Creates new File object in the database.
  Returntype : VRTrack::File object

=cut

sub create 
{
    my ($class,$vrtrack, $name) = @_;
    die "Need to call with a vrtrack handle" unless $vrtrack;
    if ( $vrtrack->isa('DBI::db') ) { croak "The interface has changed, expected vrtrack reference.\n"; }
    my $dbh = $vrtrack->{_dbh};

    # prevent adding an object with an existing name, if name supplied. In case of mapstats, the name is void
    if ($name && $class->is_name_in_database($vrtrack, $name, $name)){
        die "Already a file by name $name";
    }

    if ( !($class=~/([^:]+)$/) ) { croak "Could not determine the class name [$class]."; } 
    my $table = lc($1);

    $vrtrack->transaction_start();

    # insert a fake record to obtain a unique id (row_id)
    my $query = qq[INSERT INTO $table SET ${table}_id=0];
    my $sth   = $dbh->prepare($query) or croak qq[The query "$query" failed: $!];
    my $rv    = $sth->execute or croak qq[The query "$query" failed: $!];

    # now update the inserted the record
    my $next_id = $dbh->last_insert_id(undef,undef,$table,'row_id') or croak "No last_insert_id? $!";

    if ( $name )
    {
        my $hierarchy_name;

        my $fieldsref = $class->fields_dispatch();
        if ( exists($fieldsref->{hierarchy_name}) )
        {
            $hierarchy_name = $name;
            $hierarchy_name =~ s/\W+/_/g;
        }

        $name = qq[name='$name', ];
        if ( $hierarchy_name )
        {
            $name .= qq[hierarchy_name='$hierarchy_name', ];
        }
    }
    $query = qq[UPDATE $table SET ${table}_id=$next_id, $name changed=now(), latest=true WHERE row_id=$next_id];
    $sth   = $dbh->prepare($query) or croak qq[The query "$query" failed: $!];
    $sth->execute or croak qq[The query "$query" failed: $!];

    $vrtrack->transaction_commit();

    return $class->new($vrtrack, $next_id);
}



=head2 update

  Arg [1]    : None
  Example    : $obj->update();
  Description: Update a object whose properties you have changed.  If properties haven't changed (i.e. dirty flag is unset) do nothing.  
	       Changes the changed datestamp to now() on the mysql server (i.e. you don't have to set changed yourself, and indeed if you do, it will be overridden).
  Returntype : 1 if successful, otherwise undef.

=cut

sub update {
    my ($self) = @_;

    # database tables for core objects are named the same as the class,
    my $class = ref($self);
    $class =~/VRTrack::(\w+)$/;
    my $table = lc($1);

    my $success = undef;
    if ($self->dirty){
	    my $dbh = $self->{_dbh};
        $self->{vrtrack}->transaction_start();

        my $fieldsref = $self->fields_dispatch;
        my @fields = grep {!/^changed$/ && !/^latest$/} keys %$fieldsref;
        my $row_id; # new row_id if update works
            eval {
                # Need to unset 'latest' flag on current latest obj and add
                # the new obj details with the latest flag set
                my $updsql = qq[UPDATE $table SET latest=false WHERE ${table}_id = ? and latest=true];

                # build insert statement from update fields
                my $addsql = qq[INSERT INTO $table ( ].(join ", ", @fields);
                        $addsql .= qq[, changed, latest ) ];
                $addsql .= qq[ VALUES ( ].('?,' x scalar @fields).qq[now(),true) ];
                $dbh->do ($updsql, undef,$self->id);
                $dbh->do ($addsql, undef, map {$_->()} @$fieldsref{@fields});
                $row_id = $dbh->{'mysql_insertid'};
                $self->{vrtrack}->transaction_commit();
            };

        if ($@) {
            warn "Transaction failed, rolling back. Error was:\n$@\n";
            $self->{vrtrack}->transaction_rollback();
        }
        else {
            # Remember to update the row_id to the _new_ row_id, as this is an autoinc field
            $self->row_id($row_id);
            $success = 1;
        }
    }

    return $success;
}


=head2 delete

  Arg [1]    : none
  Example    : $obj->delete();
  Description: Deletes the current object and all its descendant objects from the database.
  Returntype : 1 if successful, otherwise undef.

=cut

sub delete {
    my ($self) = @_;
    my $success;
    my $objs_to_delete = [$self];
    push @$objs_to_delete, @{$self->descendants};
    my %rows_from_table;
    foreach (@$objs_to_delete){
        my $class=ref($_);
        $class =~/VRTrack::(\w+)$/;
        my $table = lc($1);
        $table or die "Unrecognised classname $class\n";
        push @{$rows_from_table{$table}}, $_->id;
    }

    # now delete them all, one table at a time
    my $dbh = $self->{_dbh};
    $self->{vrtrack}->transaction_start();

    eval {
        foreach my $table (keys %rows_from_table){
            my $rows = join ",", @{$rows_from_table{$table}};
            my $delsql = qq[delete from $table WHERE ${table}_id in ($rows)];
            $dbh->do ($delsql);
        }
        $self->{vrtrack}->transaction_commit();
    };

    if ($@) {
        warn "Transaction failed, rolling back. Error was:\n$@\n";
        $self->{vrtrack}->transaction_rollback();
    }
    else {
        $success = 1;
    }

    return $success;
}



=head2 note_id

  Arg [1]    : note_id (optional)
  Example    : my $note_id = $obj->note_id();
               $obj->note_id(104);
  Description: Get/Set for internal ID of the text note on this object
  Returntype : integer

=cut

sub note_id {
    my ($self,$note_id) = @_;
    if (defined $note_id and $note_id != $self->{'note_id'}){
        $self->{'note_id'} = $note_id;
        $self->dirty(1);
    }
    return $self->{'note_id'};
}


=head2 is_latest

  Arg [1]    : boolean for is_latest status
  Example    : $obj->is_latest(1);
  Description: Get/Set for object being the latest
  Returntype : boolean

=cut

sub is_latest {
    my ($self,$is_latest) = @_;
    if (defined $is_latest){
	$self->{is_latest} = $is_latest ? 1 : 0;
    }
    return $self->{is_latest};
}


=head2 row_id

  Arg [1]    : row_id (optional)
  Example    : my $row_id = $obj->row_id();
               $obj->row_id(104);
  Description: Get/Set for internal ID of the row_id of this object.
               Note that this is a database auto_increment value, so if set, is not written to the database.
  Returntype : integer

=cut

sub row_id {
    my ($self,$row_id) = @_;
    if (defined $row_id and $row_id != $self->{'row_id'}){
        $self->{'row_id'} = $row_id;
    }
    return $self->{'row_id'};
}


=head2 list_enum_vals

  Arg [1]    : table name
  Arg [2]    : column name
  Example    : my $vals = $obj->list_enum_vals('library','qc_status');
  Description: retrieves the list of allowed enum values for a column in lowercase.  Dies if the column is not of type enum
  Returntype : array ref

=cut

sub list_enum_vals {
    my ($self, $table, $column) = @_;
    my $dbh = $self->{_dbh};
    my $row = $dbh->selectrow_hashref("SHOW COLUMNS FROM $table LIKE ?", undef,$column);
    my $type = lc($row->{Type});
    unless ($type =~ /^enum/){
        die "$table:$column is not of type enum";
    }
    $type =~ s/^enum\('//;
    $type =~ s/'\)$//;
    my @vals = split /','/, $type;
    return \@vals;
}


=head2 vrtrack

  Arg [1]    : vrtrack (optional)
  Example    : my $vrtrack = $obj->vrtrack();
               $obj->vrtrack($vrtrack);
  Description: Get/Set for vrtrack object.  NB you probably really shouldn't be setting vrtrack from outside this object unless you know what you're doing.
  Returntype : integer

=cut

sub vrtrack {
    my ($self,$vrtrack) = @_;
    if (defined $vrtrack and $vrtrack != $self->{'vrtrack'}){
        $self->{_dbh} = $vrtrack->{_dbh};
        $self->{vrtrack} = $vrtrack;
    }
    return $self->{'vrtrack'};
}


=head2 allowed_processed_flags

  Example    : my %flags = VRTrack::Core_obj->allowed_processed_flags();
  Description: List allowed flags: ( import=>1, qc=>2, mapped=>4 )
  Returntype : Hash

=cut

sub allowed_processed_flags
{
    my %flags = ( import=>1, qc=>2, mapped=>4 );
    return %flags;
}

1;
