package VRTrack::Core_obj;

=head1 NAME

VRTrack::Core_obj - Sequence Tracking Core_obj object

=head1 SYNOPSIS

=head1 DESCRIPTION

This is the superclass of the core objects (Project, Sample, Library, Request,
Lane, File, Mapstats) in VRTrack, and implements some common methods.

The primary feature of a core object is that it has a history and concept of
a latest version.

=head1 CONTACT

jws@sanger.ac.uk (author)

=head1 METHODS

=cut

use strict;
use warnings;
use Carp qw(cluck confess);

use base qw(VRTrack::Table_obj);


=head2 new

  Arg [1]    : vrtrack handle
  Arg [2]    : obj id
  Example    : my $obj= $class->new($vrtrack, $id)
  Description: Returns core objects by id
  Returntype : $class object

=cut

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    return $self;
}

sub _initialize {
    my ($self, $id) = @_;
    
    my $table = $self->_class_to_table;
    
    my $fields = $self->fields_dispatch;
    my $sql = qq[select row_id,].(join ", ", keys %$fields).qq[ from $table ];
    $sql .= qq[where ${table}_id = ? and latest = true];
    
    my $sth = $self->{_dbh}->prepare($sql);
    
    if ($sth->execute($id)){
        my $data = $sth->fetchrow_hashref;
	unless ($data){
	    return;
	}
        foreach (keys %$fields){
            $fields->{$_}->($data->{$_});
        }
	
        # handle row_id as a special case outside fields_dispatch, as otherwise
        # we'd need to remove it from fields_dispatch before update.
        # Note also it is hard-coded in the select statement above
        $self->row_id($data->{'row_id'});
	
	$self->dirty(0); # unset the dirty flag
    }
    else {
	confess(sprintf('Cannot retrieve $table: %s', $DBI::errstr));
    }
}


=head2 fields_dispatch

  Arg [1]    : none
  Example    : my $fieldsref = $obj->fields_dispatch();
  Description: Returns hash dispatch table keyed on database field.
               Used internally for new and update methods. Should be extended
	       by inheriting classes to add their unique table columns and
	       methods. NB: inheriting classes should initialise their hash by
	       calling SUPER::fields_dispatch.
  Returntype : hash ref

=cut

sub fields_dispatch {
    my $self = shift;
    
    my %fields = (note_id => sub { $self->note_id(@_)},
                  changed => sub { $self->changed(@_)},
                  latest  => sub { $self->is_latest(@_)});
    
    return \%fields;
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
    my $class = shift;
    return $class->SUPER::new_by_field_value(@_);
}

sub _get_id_by_field_value {
    my ($self, $dbh, $table, $field, $value) = @_;
    
    my $sql = qq[select ${table}_id from $table where $field = ? and latest = true];
    my $sth = $dbh->prepare($sql);
    my $id;
    if ($sth->execute($value)) {
        my $data = $sth->fetchall_arrayref({}); # return array of hashes
        unless (@$data) {
            return;
        }
        if (scalar @$data > 1) {
            confess "$field = $value is not a unique identifier for $table\n";
        }
        $id = $data->[0]{"${table}_id"};
    }
    else {
        confess(sprintf('Cannot retrieve $class by %s = %s: %s', ($field, $value, $DBI::errstr)));
    }
    
    return $id;
}


=head2 create

  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : name
  Example    : my $obj = VRTrack::Core_obj->create($vrtrack, $name)
  Description: Class method.  Creates new object in the database.
  Returntype : VRTrack::Core_obj inheriting object

=cut

sub create {
    my ($class, $vrtrack, $name) = @_;
    confess "Need to call with a vrtrack handle" unless $vrtrack;
    confess "The interface has changed, expected vrtrack reference." if $vrtrack->isa('DBI::db');
    
    my $dbh = $vrtrack->{_dbh};
    my $table = $class->_class_to_table;
    
    # prevent adding an object with an existing name, if name supplied. In case of mapstats, the name is void
    if ($name && $class->is_name_in_database($vrtrack, $name, $name)){
        confess "Already a $table entry with value $name";
    }
    
    $vrtrack->transaction_start();
    
    # insert a fake record to obtain a unique id (row_id)
    my $query = qq[INSERT INTO $table SET ${table}_id=0];
    my $sth   = $dbh->prepare($query) or confess qq[The query "$query" failed: $!];
    my $rv    = $sth->execute or confess qq[The query "$query" failed: $!];
    
    # now update the inserted record
    my $next_id = $dbh->last_insert_id(undef, undef, $table, 'row_id') or confess "No last_insert_id? $!";
    
    if ($name) {
        my $hierarchy_name;
	
        my $fieldsref = $class->fields_dispatch();
        if ( exists($fieldsref->{hierarchy_name}) )
        {
            $hierarchy_name = $name;
            $hierarchy_name =~ s/\W+/_/g;
        }
	
        $name = qq[name='$name', ];
        if ($hierarchy_name) {
            $name .= qq[hierarchy_name='$hierarchy_name', ];
        }
    }
    $name ||= '';
    
    $query = qq[UPDATE $table SET ${table}_id=$next_id, $name changed=now(), latest=true WHERE row_id=$next_id];
    $sth   = $dbh->prepare($query) or confess qq[The query "$query" failed: $!];
    $sth->execute or confess qq[The query "$query" failed: $!];
    
    $vrtrack->transaction_commit();
    
    return $class->new($vrtrack, $next_id);
}


=head2 is_name_in_database

  Arg [1]    : name
  Arg [2]    : hierarchy name
  Example    : if (VRTrack::Core_obj->is_name_in_database($vrtrack, $name, $hname)
  Description: Class method. Checks to see if a name or hierarchy name is
               already used in the database table.
  Returntype : boolean

=cut

sub is_name_in_database {
    my ($class, $vrtrack, $name, $hname) = @_;
    confess "Need to call with a vrtrack handle, name, hierarchy name" unless ($vrtrack && $name && $hname);
    if ($vrtrack->isa('DBI::db')) {
	confess "The interface has changed, expected vrtrack reference.\n";
    }
    
    my $table = $class->_class_to_table;
    
    my $dbh = $vrtrack->{_dbh};
    my $sql = qq[select ${table}_id from $table where latest=true and (name = ? or hierarchy_name = ?)];
    my $sth = $dbh->prepare($sql);
    
    my $already_used = 0;
    if ($sth->execute($name, $hname)) {
        my $data = $sth->fetchrow_hashref;
        if ($data) {
            $already_used = 1;
        }
    }
    else {
        confess(sprintf('Cannot retrieve project by $name: %s', $DBI::errstr));
    }
    
    return $already_used;
}


=head2 update

  Arg [1]    : None
  Example    : $obj->update();
  Description: Update a object whose properties you have changed.  If properties
               haven't changed (i.e. dirty flag is unset) do nothing.  
	       Changes the changed datestamp to now() on the mysql server (i.e.
	       you don't have to set changed yourself, and indeed if you do,
	       it will be overridden).
  Returntype : boolean

=cut

sub update {
    my $self = shift;
    
    my $table = $self->_class_to_table;
    
    my $success = 0;
    
    if ($self->dirty) {
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
	    $dbh->do ($updsql, undef, $self->id);
	    $dbh->do ($addsql, undef, map {$_->()} @$fieldsref{@fields});
	    $row_id = $dbh->{'mysql_insertid'};
	    $self->{vrtrack}->transaction_commit();
	};
	
        if ($@) {
            cluck "Transaction failed, rolling back. Error was:\n$@\n";
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
        my $table = $class->_class_to_table;
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
        cluck "Transaction failed, rolling back. Error was:\n$@\n";
        $self->{vrtrack}->transaction_rollback();
    }
    else {
        $success = 1;
    }

    return $success;
}


=head2 changed

  Arg [1]    : changed (optional)
  Example    : my $changed = $obj->changed();
               $obj->changed('2010-01-04 10:49:10');
  Description: Get/Set for project changed
  Returntype : string

=cut

sub changed {
    my $self = shift;
    return $self->_get_set('changed', 'string', @_);
}


=head2 descendants

  Arg [1]    : none
  Example    : my $desc_objs = $obj->descendants();
  Description: Returns a ref to an array of all objects that are descendants of
               this object
  Returntype : arrayref of objects

=cut

sub descendants {
    my $self = shift;
    my @desc;
    foreach my $child_method ($self->_get_child_methods) {
	foreach my $obj (@{$self->$child_method}){
	    push @desc, $obj;
	    if ($obj->can('descendants')) {
		push @desc, @{$obj->descendants};
	    }
	}
    }
    return \@desc;
}


=head2 _get_child_methods

  Arg [1]    : none
  Example    : my @methods = $obj->_get_child_methods();
  Description: For internal use only by inheritors of Core_obj. Gets the methods
               that return a reference to an array of child objects. For use by
	       descendants().
  Returntype : array of method names

=cut

sub _get_child_methods {
    return;
}


=head2 note_id

  Arg [1]    : note_id (optional)
  Example    : my $note_id = $obj->note_id();
               $obj->note_id(104);
  Description: Get/Set for internal ID of the text note on this object
  Returntype : integer

=cut

sub note_id {
    my $self = shift;
    return $self->_get_set('note_id', 'number', @_);
}


=head2 is_latest

  Arg [1]    : boolean for is_latest status
  Example    : $obj->is_latest(1);
  Description: Get/Set for object being the latest
  Returntype : boolean

=cut

sub is_latest {
    my $self = shift;
    return $self->_get_set('is_latest', 'boolean', @_);
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
    my $self = shift;
    return $self->_get_set('row_id', 'number', @_);
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
    my $row = $dbh->selectrow_hashref("SHOW COLUMNS FROM $table LIKE ?", undef, $column);
    my $type = lc($row->{Type});
    unless ($type =~ /^enum/){
        confess "$table:$column is not of type enum";
    }
    $type =~ s/^enum\('//;
    $type =~ s/'\)$//;
    my @vals = split /','/, $type;
    return \@vals;
}


sub _check_status_value {
    my ($self, $type, $value) = @_;
    
    if (defined $value) {
	my $class = ref($self);
	confess "Could not determine the class name [$class]." unless $class=~/([^:]+)$/;
	my $table = lc($1);
	
        my %allowed = map {$_ => 1} @{$self->list_enum_vals($table, $type)};
        unless ($allowed{lc($value)}){
            confess "'$value' is not a defined $type";
        }
    }
}


=head2 vrtrack

  Arg [1]    : vrtrack (optional)
  Example    : my $vrtrack = $obj->vrtrack();
               $obj->vrtrack($vrtrack);
  Description: Get/Set for vrtrack object.  NB you probably really shouldn't be setting vrtrack from outside this object unless you know what you're doing.
  Returntype : integer

=cut

sub vrtrack {
    my ($self, $vrtrack) = @_;
    if (defined $vrtrack and $vrtrack != $self->{vrtrack}) {
        $self->{_dbh} = $vrtrack->{_dbh};
        $self->{vrtrack} = $vrtrack;
    }
    return $self->{vrtrack};
}


=head2 allowed_processed_flags

  Example    : my %flags = VRTrack::Core_obj->allowed_processed_flags();
  Description: Get a hash of allowed processed flags and their values.
  Returntype : Hash

=cut

sub allowed_processed_flags {
    my %flags = (import => 1,
		 qc => 2,
		 mapped => 4,
		 stored => 8,
		 deleted => 16,
		 swapped => 32,
		 altered_fastq => 64);
    return %flags;
}

1;
