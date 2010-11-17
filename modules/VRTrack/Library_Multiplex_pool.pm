package VRTrack::Library_Multiplex_pool; 

=head1 NAME

VRTrack::Library_Multiplex_library  - Sequence Tracking Library_Multiplex_library object

=head1 SYNOPSIS
    my $lib = VRTrack::Library_Multiplex_library->new($vrtrack, $library_id,$multiplex_library);
    
=head1 DESCRIPTION

An object describing the relationship b/w library and multiplexed_library.

=head1 AUTHOR

rn2@sanger.ac.uk (author)

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

use base qw(VRTrack::Table_obj);


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
               library_multiplex_pool_id        => sub { $self->id(@_)},
               multiplex_pool_id         => sub { $self->multiplex_pool_id(@_)},
               library_id              => sub { $self->library_id(@_)});
    return \%fields;
}


###############################################################################
# Class methods
###############################################################################

=head2 new

  Arg [1]    : vrtrack handle
  Arg [2]    : obj id
  Example    : my $obj= $class->new($vrtrack, $id)
  Description: Returns object by id
  Returntype : $class object or undef if no such object id

=cut

=head2 new_by_field_value

  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : field name
  Arg [3]    : field value
  Example    : my $obj = $class->new_by_field_value($vrtrack, 'name',$name)
  Description: Class method. Returns latest $class object by field name and value.  If no such value is in the database, returns undef.
               Dies if there is more than one matching record.
  Returntype : $class object

=cut

=head2 create
  
  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : ssid
  Example    : my $file = VRTrack::Library_Multiplex_library->create($vrtrack, $ssid)
  Description: Class method.  Creates new Library_Multiplex Library object in the database.
  Returntype : VRTrack::Library_Multiplex_library object
   
=cut

sub create{
my ($class, $vrtrack,$library_id,$multiplex_pool_id) = @_;
    my $obj=$class->SUPER::create($vrtrack,('library_id'=>$library_id,'multiplex_pool_id'=>$multiplex_pool_id));
    return $obj;
}


###############################################################################
# Object methods
###############################################################################

=head2 dirty

  Arg [1]    : boolean for dirty status
  Example    : $obj->dirty(1);
  Description: Get/Set for object properties having been altered.
  Returntype : boolean

=cut

=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $obj->id();
               $obj->id('104');
  Description: Get/Set for ID of a VRTrack object
  Returntype : Internal ID integer

=cut

=head2 vrtrack

  Arg [1]    : vrtrack (optional)
  Example    : my $vrtrack = $obj->vrtrack();
               $obj->vrtrack($vrtrack);
  Description: Get/Set for vrtrack object.  NB you probably really shouldn't be
               setting vrtrack from outside this object unless you know what
	       you're doing.
  Returntype : VRTrack::VRTrack

=cut


=head2 new_by_library_id_multiplex_pool_id
  
  Arg [1]    : library_id 
  Arg [2]    : multiplex_pool_id
  Example    : my $vlib_multiplex_pool=VRTrack::Library_Multiplex_pool->new_by_library_id_multiplex_pool_id($vrtrack,$vlib->id,$vmultiplex_pool->id);
  Description: Class method.  Returns a library_multiplex_pool record associated with this library_id and the multiplex_pool_id.
  Returntype : VRTrack::Library_Multiplex_pool object
   
=cut

sub new_by_library_id_multiplex_pool_id{
    my ($class, $vrtrack, $library_id,$multiplex_pool_id) = @_;
    confess "Need to call with a vrtrack handle, library_id, multiplex_pool_id" unless ($vrtrack && $library_id && $multiplex_pool_id);
    confess "The interface has changed, expected vrtrack reference." if $vrtrack->isa('DBI::db');
    
    my $dbh = $vrtrack->{_dbh};
    my $table = $class->_class_to_table;

    my $id = $class->_get_id_by_library_id_multiplex_pool_id($dbh,$library_id,$multiplex_pool_id);
    return unless $id;
    return $class->new($vrtrack, $id);
}    
    


=head2 _get_id_by_library_id_multiplex_pool_id
  
  Arg [1]    : library_id
  Arg [2]    : multiplex_pool_id
  Example    : my $id = $class->_get_id_by_library_id_multiplex_pool_id($dbh,$library_id,$multiplex_pool_id);
  Description: Class method.  Returns the library_multiplex_pool_id associated with the library_id and multiplex_pool_id
  Returntype : library_multiplex_pool_id
   
=cut

sub _get_id_by_library_id_multiplex_pool_id{
    my ($self, $dbh,$library_id,$multiplex_pool_id) = @_;
    
    my $table='library_multiplex_pool';
    my $sql = qq[select library_multiplex_pool_id from $table where library_id=$library_id and multiplex_pool_id=$multiplex_pool_id];
    my $sth = $dbh->prepare($sql);
    my $id;
    if ($sth->execute) {
        my $data = $sth->fetchall_arrayref({}); # return array of hashes
        unless (@$data) {
            return;
        }
        if (scalar @$data > 1) {
            confess " more than 1 records returned for $table\n";
        }
        $id = $data->[0]{"${table}_id"};
    }
    else {
        confess "Cannot retrieve $table by $library_id,$multiplex_pool_id: ".$DBI::errstr;
    }
    
    return $id;
}    

=head2 library_id

  Arg [1]    : library_id (optional)
  Example    : my $library_id = $lib_multiplex_lib->library_id();
               $lib_multiplex_lib->library_id(104);
  Description: Get/Set for ID of a library
  Returntype : Internal ID integer

=cut

sub library_id {
    my $self = shift;
    return $self->_get_set('library_id', 'number', @_);
}


=head2  multiplex_pool_id

  Arg [1]    : multiplex_pool_id (optional)
  Example    : my $multiplex_pool_id = $lib_multiplex_lib->multiplex_pool_id();
               $lib_multiplex_lib->multiplex_pool_id(104);
  Description: Get/Set for ID of a library
  Returntype : Internal ID integer

=cut

sub multiplex_pool_id {
    my $self = shift;
    return $self->_get_set('multiplex_pool_id', 'number', @_);
}

1;
