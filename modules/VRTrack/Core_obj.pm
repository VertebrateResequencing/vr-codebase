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
no warnings 'uninitialized';


=head2 new

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : obj id
  Example    : my $obj= $class->new($dbh, $id)
  Description: Returns core objects by id
  Returntype : $class object

=cut

sub new {
    my ($class, $dbh, $id) = @_;
    die "Need to call with a db handle and id" unless ($dbh && $id);
    my $self = {};
    bless ($self, $class);
    $self->{_dbh} = $dbh;
    
    # database tables for core objects are named the same as the class,
    $class =~/VRTrack::(\w+)$/;
    my $table = lc($1);
    $table or die "Unrecognised classname $class\n";

    my $fieldsref = $self->fields_dispatch;
    my $sql = qq[select ].(join ", ",keys %$fieldsref).qq[ from $table ];
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

	$self->dirty(0);    # unset the dirty flag
    }
    else{
	die(sprintf('Cannot retrieve $table: %s', $DBI::errstr));
    }

    return $self;
}


=head2 new_by_field_value

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : field name
  Arg [3]    : field value
  Example    : my $obj = $class->new_by_field_value($dbh, 'name',$name)
  Description: Class method. Returns latest $class object by field name and value.  If no such value is in the database, returns undef.
               Dies if there is more than one matching record.
  Returntype : $class object

=cut

sub new_by_field_value {
    my ($class,$dbh, $field, $value) = @_;
    die "Need to call with a db handle, field name, field value" unless ($dbh && $field && defined $value);
    
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
    return $class->new($dbh, $id);
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
	my $save_re = $dbh->{RaiseError};
	my $save_pe = $dbh->{PrintError};
	my $save_ac = $dbh->{AutoCommit};
	$dbh->{RaiseError} = 1; # raise exception if an error occurs
	$dbh->{PrintError} = 0; # don't print an error message
	$dbh->{AutoCommit} = 0; # disable auto-commit

        my $fieldsref = $self->fields_dispatch;
        my @fields = grep {!/^changed$/ && !/^latest$/} keys %$fieldsref;
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
	    $dbh->commit ( );
	};

	if ($@) {
	    warn "Transaction failed, rolling back. Error was:\n$@\n";
	    # roll back within eval to prevent rollback
	    # failure from terminating the script
	    eval { $dbh->rollback ( ); };
	}
	else {
	    $success = 1;
	}

	# restore attributes to original state
	$dbh->{AutoCommit} = $save_ac;
	$dbh->{PrintError} = $save_pe;
	$dbh->{RaiseError} = $save_re;

    }

    return $success;
}


=head2 delete

  Arg [1]    : none
  Example    : $obj->delete();
  Description: Deletes the current object and all its descendant objects from the database.  R
  Returntype : 1 if successful, otherwise undef.

=cut

sub delete {
    my ($self) = @_;
    my $success;
    my $objs_to_delete = [$self];
    push @$objs_to_delete, @{$self->descendants};
    my %lanes_from_table;
    foreach (@$objs_to_delete){
        my $class=ref($_);
        $class =~/VRTrack::(\w+)$/;
        my $table = lc($1);
        $table or die "Unrecognised classname $class\n";
        push @{$lanes_from_table{$table}}, $_->id;
    }

    # now delete them all, one table at a time
    my $dbh = $self->{_dbh};
    my $save_re = $dbh->{RaiseError};
    my $save_pe = $dbh->{PrintError};
    my $save_ac = $dbh->{AutoCommit};
    $dbh->{RaiseError} = 1; # raise exception if an error occurs
    $dbh->{PrintError} = 0; # don't print an error message
    $dbh->{AutoCommit} = 0; # disable auto-commit

    eval {
        foreach my $table (keys %lanes_from_table){
            my $lanes = join ",", @{$lanes_from_table{$table}};
            my $delsql = qq[delete from $table WHERE ${table}_id in ($lanes)];
            $dbh->do ($delsql);
        }
        $dbh->commit ( );
    };

    if ($@) {
        warn "Transaction failed, rolling back. Error was:\n$@\n";
        # roll back within eval to prevent rollback
        # failure from terminating the script
        eval { $dbh->rollback ( ); };
    }
    else {
        $success = 1;
    }

    # restore attributes to original state
    $dbh->{AutoCommit} = $save_ac;
    $dbh->{PrintError} = $save_pe;
    $dbh->{RaiseError} = $save_re;

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



1;
