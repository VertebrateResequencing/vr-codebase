package VRTrack::File; 

=head1 NAME

VRTrack::File - Sequence Tracking File object

=head1 SYNOPSIS
    my $file= VRTrack::File->new($vrtrack, $file_id);

    my $id = $file->id();

=head1 DESCRIPTION

An object describing the tracked properties of a file.

=head1 AUTHOR

jws@sanger.ac.uk (author)

=head1 METHODS

=cut

use strict;
use warnings;
use Carp qw(cluck confess);

use base qw(VRTrack::Core_obj
            VRTrack::Hierarchy_obj
	    VRTrack::Named_obj);


=head2 fields_dispatch

  Arg [1]    : none
  Example    : my $fieldsref = $file->fields_dispatch();
  Description: Returns hashref dispatch table keyed on database field
               Used internally for new and update methods
  Returntype : hashref

=cut

sub fields_dispatch {
    my $self = shift;
    
    my %fields = %{$self->SUPER::fields_dispatch()};
    %fields = (%fields,
               file_id           => sub { $self->id(@_)},
               lane_id           => sub { $self->lane_id(@_)},
               name              => sub { $self->name(@_)},
               hierarchy_name    => sub { $self->hierarchy_name(@_)},
               processed         => sub { $self->processed(@_)},
               type              => sub { $self->type(@_)},
               readlen           => sub { $self->read_len(@_)},
               raw_reads         => sub { $self->raw_reads(@_)},
               raw_bases         => sub { $self->raw_bases(@_)},
               mean_q            => sub { $self->mean_q(@_)},
               md5               => sub { $self->md5(@_)});
    
    return \%fields;
}


###############################################################################
# Class methods
###############################################################################

=head2 new_by_name

  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : file name
  Example    : my $file = VRTrack::File->new_by_name($vrtrack, $name)
  Description: Class method. Returns latest File object by name.  If no such name is in the database, returns undef.  Dies if multiple names match.
  Returntype : VRTrack::File object

=cut


=head2 new_by_hierarchy_name

  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : file hierarchy_name
  Example    : my $file = VRTrack::File->new_by_hierarchy_name($vrtrack, $hierarchy_name)
  Description: Class method. Returns latest File object by hierarchy_name.  If no such hierarchy_name is in the database, returns undef.  Dies if multiple hierarchy_names match.
  Returntype : VRTrack::File object

=cut


=head2 create
  
  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : file name
  Example    : my $file = VRTrack::File->create($vrtrack, $name)
  Description: Class method.  Creates new File object in the database.
  Returntype : VRTrack::File object
   
=cut


=head2 is_name_in_database

  Arg [1]    : file name
  Arg [2]    : hierarchy name
  Example    : if(VRTrack::File->is_name_in_database($vrtrack, $name,$hname)
  Description: Class method. Checks to see if a name or hierarchy name is already used in the file table.
  Returntype : boolean

=cut


###############################################################################
# Object methods
###############################################################################

=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $lib->id();
               $lib->id(104);
  Description: Get/Set for internal db ID of a lane
  Returntype : integer

=cut


=head2 hierarchy_name

  Arg [1]    : directory name (optional)
  Example    : my $hname = $file->hierarchy_name();
  Description: Get/set file hierarchy name.  This is the directory name (without path) that the file will be named in a file hierarchy.
  Returntype : string

=cut


=head2 name

  Arg [1]    : name (optional)
  Example    : my $name = $file->name();
	       $file->name('1111_s_3_1.fastq');
  Description: Get/Set for name of a file
  Returntype : string

=cut


=head2 md5

  Arg [1]    : md5 (optional)
  Example    : my $md5 = $file->md5();
	       $file->md5('1027759a77eb562fcab253c9f01d7661');
  Description: Get/Set for md5 of a file
  Returntype : string

=cut

sub md5 {
    my $self = shift;
    return $self->_get_set('md5', 'string', @_);
}


=head2 lane_id

  Arg [1]    : lane_id (optional)
  Example    : my $lane_id = $file->lane_id();
	       $file->lane_id('104');
  Description: Get/Set for ID of a file
  Returntype : SequenceScape ID (usu. integer)

=cut

sub lane_id {
    my $self = shift;
    return $self->_get_set('lane_id', 'number', @_);
}


=head2 changed

  Arg [1]    : changed (optional)
  Example    : my $changed = $file->changed();
	       $file->changed('20090101102300');
  Description: Get/Set for file changed
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

  Arg [1]    : One of flags listed in Core_obj::allowed_processed_flags();
  Arg [2]    : processed: 0 or 1 (optional)
  Example    : my $processed = $file->is_processed('import');
               $file->is_processed('import',1);
  Description: Get/Set for file processed
  Returntype : 1 or 0

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


=head2 type

  Arg [1]    : type (optional) [0,1,2,3]
  Example    : my $type = $file->type();
	       $file->type(1);
  Description: Get/Set for file type - 0 is single-end, 1 fwd, 2 rev
                3 is for _s_ files, i.e. fwd & rev concatenated on one
                line
  Returntype : integer

=cut

sub type {
    my $self = shift;
    return $self->_get_set('type', 'number', @_);
}


=head2 raw_reads

  Arg [1]    : number of sequences in file
  Example    : my $num_seqs = $file->raw_reads();
	       $file->raw_reads(1_000_000);
  Description: Get/Set for number of sequences in a fasta file
  Returntype : integer

=cut

sub raw_reads {
    my $self = shift;
    return $self->_get_set('seq_count', 'number', @_);
}


=head2 raw_bases

  Arg [1]    : number of basepairs in file
  Example    : my $num_bp = $file->raw_bases();
	       $file->raw_bases(1_000_000);
  Description: Get/Set for total number of unfiltered basepairs in file
  Returntype : integer

=cut

sub raw_bases {
    my $self = shift;
    return $self->_get_set('raw_bases', 'number', @_);
}


=head2 mean_q

  Arg [1]    : mean quality score for bases in file
  Example    : my $mean_q = $file->mean_q();
	       $file->mean_q(27);
  Description: Get/Set for mean quality score for bases in file
  Returntype : float

=cut

sub mean_q {
    my $self = shift;
    return $self->_get_set('mean_q', 'number', @_);
}


=head2 read_len

  Arg [1]    : read length
  Example    : my $readlen = $file->read_len();
	       $file->read_len(54);
  Description: Get/Set for unclipped read length in file
  Returntype : integer

=cut

sub read_len {
    my $self = shift;
    return $self->_get_set('readlen', 'number', @_);
}

1;
