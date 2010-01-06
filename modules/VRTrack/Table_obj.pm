package VRTrack::Table_obj; 

=head1 NAME

VRTrack::Base_obj - Sequence Tracking Table_obj object

=head1 SYNOPSIS

=head1 DESCRIPTION

This is the superclass of all VRTrack modules that represent a table in the
database, providing general methods.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;
use Carp qw(cluck confess);
use Scalar::Util qw(looks_like_number);


=head2 new

  Arg [1]    : vrtrack handle
  Arg [2]    : obj id
  Example    : my $obj= $class->new($vrtrack, $id)
  Description: Returns object by id
  Returntype : $class object

=cut

sub new {
    my ($class, $vrtrack, $id) = @_;
    
    my $self = {};
    bless $self, ref($class) || $class;
    
    if ($vrtrack) {
	confess "The interface has changed, expected vrtrack reference." unless $vrtrack->isa('VRTrack::VRTrack');
	confess "Need to call with a vrtrack reference and id" unless ($vrtrack && $id);
	
	$self->{_dbh} = $vrtrack->{_dbh};
	$self->{vrtrack} = $vrtrack;
	
	$self->_initialize($id);
    }
    
    return $self;
}

sub _initialize {
    my ($self, $id) = @_;
    
    my $table = $self->_class_to_table;
    
    my $fields = $self->fields_dispatch;
    my $sql = qq[select ].(join ", ", keys %$fields).qq[ from $table ];
    $sql .= qq[where ${table}_id = ?];
    
    my $sth = $self->{_dbh}->prepare($sql);
    
    if ($sth->execute($id)){
        my $data = $sth->fetchrow_hashref;
	unless ($data){
	    return undef;
	}
        foreach (keys %$fields){
            $fields->{$_}->($data->{$_});
        }
	
	$self->dirty(0); # unset the dirty flag
    }
    else {
	confess(sprintf('Cannot retrieve $table: %s', $DBI::errstr));
	$DBI::errstr;
    }
}


# database tables for objects are named the same as the class
sub _class_to_table {
    my $thing = shift;
    my $class = ref($thing) || $thing;
    confess "Could not determine the class name [$class]." unless $class =~ /([^:]+)$/; 
    return lc($1);
}


=head2 fields_dispatch

  Arg [1]    : none
  Example    : my $fieldsref = $obj->fields_dispatch();
  Description: Returns hash dispatch table keyed on database field.
               Used internally for new and update methods. Should be extended
	       by inheriting classes to add their unique table columns and
	       methods.
  Returntype : hash ref

=cut

sub fields_dispatch {
    return {};
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
    my ($class, $vrtrack, $field, $value) = @_;
    confess "Need to call with a vrtrack handle, field name, field value" unless ($vrtrack && $field && defined $value);
    confess "The interface has changed, expected vrtrack reference." if $vrtrack->isa('DBI::db');
    
    my $dbh = $vrtrack->{_dbh};
    my $table = $class->_class_to_table;

    # check field exists
    my $colnames = $dbh->selectcol_arrayref(qq[select column_name from information_schema.columns where table_name='$table']);
    my %cols = map { $_ => 1 } @$colnames;
    unless (exists($cols{lc($field)})) {
        confess "No such column $field in $table table\n";
    }

    my $id = $class->_get_id_by_field_value($dbh, $table, $field, $value);
    return unless $id;
    return $class->new($vrtrack, $id);
}

sub _get_id_by_field_value {
    my ($self, $dbh, $table, $field, $value) = @_;
    
    my $sql = qq[select ${table}_id from $table where $field = ?];
    my $sth = $dbh->prepare($sql);
    my $id;
    if ($sth->execute($value)) {
        my $data = $sth->fetchall_arrayref({}); # return array of hashes
        unless (@$data) {
            return undef;
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
  Arg [2]    : hash of column => value pairs, excluding the autoincrement column
  Example    : my $obj = $class->create($vrtrack, name => 'unique_name')
  Description: Class method. Creates new object in the database.
  Returntype : VRTrack::Table_obj inheriting object

=cut

sub create {
    my ($class, $vrtrack, %data) = @_;
    confess "Need to call with a vrtrack handle" unless $vrtrack;
    confess "The interface has changed, expected vrtrack reference." if $vrtrack->isa('DBI::db');
    
    my $dbh = $vrtrack->{_dbh};
    my $table = $class->_class_to_table;
    my $auto_increment = $table.'_id';
    delete $data{$auto_increment};
    
    my ($qms, @cols, @values);
    while (my ($col, $val) = each %data) {
	next unless defined $val;
	push(@cols, $col);
	push(@values, $val);
	$qms .= ',?';
    }
    confess "at least one column => value pair other than the autoincrement is required" unless @cols;
    my $cols = join(',', @cols);
    
    my $sql = qq[INSERT INTO $table ($auto_increment, $cols) VALUES (NULL$qms)];
    
    my $sth = $dbh->prepare($sql);
    my $id;
    if ($sth->execute(@values)) {
        $id = $dbh->{'mysql_insertid'};
    }
    else {
        confess( sprintf('DB load insert failed: %s', $DBI::errstr));
    }
    
    return $class->new($vrtrack, $id);
}


=head2 update

  Arg [1]    : None
  Example    : $study->update();
  Description: Update a study whose properties you have changed.
               If properties haven't changed (i.e. dirty flag is unset) do
	       nothing. Unsets the dirty flag on success.
  Returntype : boolean

=cut

sub update {
    my $self = shift;
    
    my $success = 0;
    
    if ($self->dirty) {
	my $dbh = $self->{_dbh};
	my $table = $self->_class_to_table;
	my $key = $table.'_id';
	my $fieldsref = $self->fields_dispatch;
        my @fields = grep {!/^$key$/} keys %$fieldsref;
	
	my $save_re = $dbh->{RaiseError};
	my $save_pe = $dbh->{PrintError};
	$dbh->{RaiseError} = 1; # raise exception if an error occurs
	$dbh->{PrintError} = 0; # don't print an error message
	
	my $updsql = qq[UPDATE $table SET ];
	my @sql_bits;
	foreach my $field (@fields) {
	    push(@sql_bits, "$field=?");
	}
	$updsql .= join(', ', @sql_bits).qq[ WHERE $key = ? ];
	
	my @values = map {$_->()} @$fieldsref{@fields};
	eval {
	    $dbh->do($updsql, undef, @values, $self->id);
	};
	
	if ($@) {
	    cluck "update failed; error was:\n$@\n\noriginal sql was:\n$updsql\n";
	}
	else {
	    $success = 1;
	    $self->dirty(0);
	}
	
	# restore attributes to original state
	$dbh->{PrintError} = $save_pe;
	$dbh->{RaiseError} = $save_re;
    }
    
    return $success;
}


=head2 _get_set_child_object

  Arg [1]    : name of method used to get existing child object
  Arg [2]    : class of the child
  Arg [3]    : identifier of child to find/create (optional; can be a list if
               Arg[1] takes a list)
  Example    : my $child = $obj->_get_set_child_object('get_study_by_acc', 'VRTrack::Study', 'SRP000031');
  Description: For internal use only by inheritors of Core_obj. Gets/sets a
               child object on this object. Does not create child objects if
	       they don't already exist in the db.
  Returntype : instance of child

=cut

sub _get_set_child_object {
    my ($self, $existing_method, $child_class, @child_identifiers) = @_;
    
    my $child_storage_key = $child_class->_class_to_table;
    my $child_id_method = $child_storage_key.'_id';
    
    if (@child_identifiers) {
        # get existing child by name
        my $obj = $self->$existing_method(@child_identifiers);
	
        if ($obj) {
            # Have we actually changed?
            if ($self->$child_id_method != $obj->id) {
                $self->$child_id_method($obj->id);
                $self->dirty(1);
            }
            $self->{$child_storage_key} = $obj;
        }
        else {
            return undef; # explicitly return nothing.
        }
    }
    elsif (defined $self->{$child_storage_key}) {
        # check that our current id matches the object we've cached
	my $obj = $self->{$child_storage_key};
	if ($self->$child_id_method != $obj->id) {
	    undef $self->{$child_storage_key};
	}
    }
    
    if (! defined $self->{$child_storage_key}) {
        # lazy-load child from database
        if ($self->$child_id_method) {
            my $obj = $child_class->new($self->{vrtrack}, $self->$child_id_method);
            $self->{$child_storage_key} = $obj;
        }
    }
    
    return $self->{$child_storage_key};
}

=head2 _create_child_object

  Arg [1]    : name of method used to get existing child object
  Arg [2]    : class of the child
  Arg [3]    : identifier of child to create (can supply list of ids as needed
               by Arg [1])
  Example    : my $child = $obj->_create_child_object('get_study_by_acc', 'VRTrack::Study', 'SRP000031');
  Description: For internal use only by inheritors of Core_obj. Creates a
               child object on this object. Returns undef if the child already
	       exists in the db. For use only on child types where only one
	       child is allowed.
  Returntype : instance of child

=cut

sub _create_child_object {
    my ($self, $existing_method, $child_class, @child_identifiers) = @_;
    
    my $child_storage_key = $child_class->_class_to_table;
    my $child_id_storage_key = $child_storage_key.'_id';
    
    my $obj = $self->$existing_method(@child_identifiers);
    if ($obj) {
        cluck "$child_class (@child_identifiers) is already present in the database\n";
        return undef;
    }
    else {
	$obj = $child_class->create($self->{vrtrack}, @child_identifiers);
	
        # populate caches
        $self->{$child_id_storage_key} = $obj->id;
        $self->{$child_storage_key} = $obj;
        $self->dirty(1);
    }
    
    return $self->{$child_storage_key};
}


=head2 _get_child_objects

  Arg [1]    : class of the child
  Example    : my $children = $obj->_get_child_objects('VRTrack::Sample');
  Description: For internal use only by inheritors of Core_obj. Gets the child
               objects of the given class that have us as a parent.
  Returntype : array ref of children

=cut

sub _get_child_objects {
    my ($self, $child_class) = @_;
    
    my $table = $child_class->_class_to_table;
    my $children_storage_key = $table.'s';
    my $ids_method = $table.'_ids';
    
    unless ($self->{$children_storage_key}){
        my @objs;
        foreach my $id (@{$self->$ids_method()}){
            my $obj = $child_class->new($self->{vrtrack}, $id);
            push @objs, $obj;
        }
        $self->{$children_storage_key} = \@objs;
    }
    
    return $self->{$children_storage_key};
}


=head2 _get_child_ids

  Arg [1]    : class of the child
  Example    : my $children_ids = $obj->_get_child_ids('VRTrack::Sample');
  Description: For internal use only by inheritors of Core_obj. Gets the child
               ids of the given class that have us as a parent. 
  Returntype : array ref of children ids

=cut

sub _get_child_ids {
    my ($self, $child_class) = @_;
    
    my $table = $child_class->_class_to_table;
    my $child_table_column = $table.'_id';
    my $children_ids_storage_key = $child_table_column.'s';
    my $child_parent_id_column = $self->_class_to_table.'_id';
    
    unless ($self->{$children_ids_storage_key}) {
	my $history_sql = '';
	if ($child_class->isa('VRTrack::Core_obj')) {
	    $history_sql = ' and latest=true';
	}
	
        my $sql = qq[select distinct($child_table_column) from $table where $child_parent_id_column=?$history_sql];
        my @objs;
        my $sth = $self->{_dbh}->prepare($sql);
	
        if ($sth->execute($self->id)) {
            foreach(@{$sth->fetchall_arrayref()}){
                push @objs, $_->[0];
            }
        }
        else {
            confess(sprintf('Cannot retrieve child instances: %s', $DBI::errstr));
        }
	
        $self->{$children_ids_storage_key} = \@objs;
    }
    
    return $self->{$children_ids_storage_key};
}


=head2 _add_child_object

  Arg [1]    : name of method used to get existing child object (optional)
  Arg [2]    : class of the child
  Arg [3]    : identifier of child to create (can supply list of ids as needed
               by Arg [1]; optional if Arg[1] is undef)
  Example    : my $child = $obj->_add_child_object('new_by_name_project', 'VRTrack::Sample', $sample_name, $project_id);
  Description: For internal use only by inheritors of Core_obj. Adds a
               child object to this object. Returns undef if the child already
	       exists in the db. For use only on child types where many children
	       are stored at once.
  Returntype : instance of child

=cut

sub _add_child_object {
    my ($self, $existing_method, $child_class, @child_identifiers) = @_;
    
    my $table = $child_class->_class_to_table;
    my $children_storage_key = $table.'s';
    my $children_ids_storage_key = $table.'_ids';
    my $child_parent_id_method = $self->_class_to_table.'_id';
    
    if ($existing_method) {
	@child_identifiers > 0 || confess "Must call with a child identifier";
	
	my $obj = $child_class->$existing_method($self->{vrtrack}, @child_identifiers);
	if ($obj) {
	    cluck "$child_class (@child_identifiers) is already present in the database\n";
	    return undef;
	}
    }
    
    my $obj = $child_class->create($self->{vrtrack}, @child_identifiers);
    if ($obj) {
        $obj->$child_parent_id_method($self->id);
	if ($obj->can('hierarchy_name')) {
	    my $hierarchy_name = $child_identifiers[0];
	    $hierarchy_name =~ s/\W+/_/g;
	    $obj->hierarchy_name($hierarchy_name);
	}
        $obj->update;
    }
    delete $self->{$children_ids_storage_key};
    delete $self->{$children_storage_key};
    
    return $obj;
}


=head2 _get_child_by_field_value

  Arg [1]    : name of method used to get all of the desired children
  Arg [2]    : field (name of method on child that returns a value to match against)
  Arg [3]    : value (to search for amongst children)
  Example    : my $child = $obj->_get_child_by_field_value('samples', 'name', 'mysample_name');
  Description: For internal use only by inheritors of Core_obj. Gets one of the
               descendant objects by field value. Dies if more than one child
	       matches the value. Only for use on child types where many
	       children are stored at once.
  Returntype : instance of child

=cut

sub _get_child_by_field_value {
    my ($self, $get_children_method, $field, $value) = @_;
    $value || confess "Must call with a child identifier";
    
    my @match = grep {$_->$field eq $value} @{$self->$get_children_method};
    
    if (@match == 1) {
        return $match[0];
    }
    elsif (@match > 1) {
	confess "More than one child with $field $value";
    }

    return;
}


=head2 dirty

  Arg [1]    : boolean for dirty status
  Example    : $obj->dirty(1);
  Description: Get/Set for object properties having been altered.
  Returntype : boolean

=cut

sub dirty {
    # (dirty can't use _get_set because _get_set calls dirty)
    my ($self, $dirty) = @_;
    if (defined $dirty) {
	$self->{_dirty} = $dirty ? 1 : 0;
    }
    return $self->{_dirty};
}


=head2 _get_set

  Arg [1]    : a self hash key to get values from or store them to
  Arg [2]    : expected type of data: 'number', 'boolean' or 'string'
  Arg [3]    : a value to set (optional)
  Example    : my $note_id = $obj->_get_set('note_id', 'number');
               $obj->_get_set('note_id', 'number', 104);
  Description: Internal generic get/setter for modules inheriting from us;
               not for end-user use.
  Returntype : int or string

=cut

sub _get_set {
    my ($self, $key, $type, $value) = @_;
    
    if (defined $value) {
	if ($type eq 'number') {
	    confess "value to set must be a number" unless looks_like_number($value);
	}
	elsif ($type eq 'boolean') {
	    $value = $value ? 1 : 0;
	}
	
	if (! defined $self->{$key} || "$self->{$key}" ne "$value") {
	    $self->{$key} = $value;
	    $self->dirty(1);
	}
    }
    
    return $self->{$key};
}


=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $obj->id();
               $obj->id('104');
  Description: Get/Set for ID of a VRTrack object
  Returntype : Internal ID integer

=cut

sub id {
    my $self = shift;
    return $self->_get_set('id', 'number', @_);
}


=head2 vrtrack

  Arg [1]    : vrtrack (optional)
  Example    : my $vrtrack = $obj->vrtrack();
               $obj->vrtrack($vrtrack);
  Description: Get/Set for vrtrack object.  NB you probably really shouldn't be
               setting vrtrack from outside this object unless you know what
	       you're doing.
  Returntype : VRTrack::VRTrack

=cut

sub vrtrack {
    my ($self, $vrtrack) = @_;
    
    if (defined $vrtrack && (! defined $self->{vrtrack} || "$vrtrack" ne "$self->{vrtrack}")){
        $self->{_dbh} = $vrtrack->{_dbh};
        $self->{vrtrack} = $vrtrack;
    }
    
    return $self->{vrtrack};
}

1;
