package Sfind::Library; 
=head1 NAME

Sfind::Library - Sequence Tracking Library object

=head1 SYNOPSIS
    my $lib = Sfind::Library->new({dbh => $dbh, id => $library_id});

    #get arrayref of requests on a library
    my $reqs = $library->requests();

    my $id = $library->id();
    my $name = $library->name();


=head1 DESCRIPTION

An object describing the tracked properties of a library.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=head2 new

  Arg [1]    : hashref: dbh => database handle to seqtracking database
                        id  => library id
  Example    : my $lib = Sfind::Library->new({dbh=>$dbh, id=>$id};)
  Description: Returns Library object by library_id
  Returntype : Sfind::Library object

=cut

use Moose;
use namespace::autoclean;
use Sfind::Types qw(MysqlDateTime);
use Sfind::Seq_Request;
use Sfind::Lane;

# Populate the parameters from the database
around BUILDARGS => sub {
    my $orig  = shift;
    my $class = shift;
    
    my $argref = $class->$orig(@_);

    die "Need to call with a library id" unless $argref->{id};
    my $sql = qq[select * from requests where internal_id = ? and is_current=1];
    my $id_ref = $argref->{dbh}->selectrow_hashref($sql, undef, ($argref->{id}));
    if ($id_ref){
        foreach my $field(keys %$id_ref){
            $argref->{$field} = $id_ref->{$field};
        }
    };
    return $argref;
};
sub new {
    my ($class,$dbh, $id) = @_;
    die "Need to call with a db handle and id" unless ($dbh && $id);
    my $self = {};
    bless ($self, $class);
    $self->{_dbh} = $dbh;

    $self->id($id);
    
    # name
    my $sql = qq[select name from assets where asset_id=?];
    my $id_ref = $self->{_dbh}->selectrow_hashref($sql, undef, ($id));
    unless ($id_ref){
        die "No name for asset $id";
    }
    my $name = $id_ref->{name};
    $name =~ s/\s+$//;  # trim trailing whitespace

    # 2010-10-22 problem with warehouse - item_information is missing.
    # So, let's put in a hack to at least create a library if there is no name
    # hopefully it will get auto-fixed when the item_information comes back
    $name ||= $id;

    $self->name($name);

    # prep_status
    $sql = qq[select state from requests_new where target_asset_id=?];
    $id_ref = $self->{_dbh}->selectrow_hashref($sql, undef, ($id));
    if ($id_ref){
        $self->prep_status($id_ref->{state});
    }


    # Are we a tagged library?
    # If so, we should have two parents - a sample, and a taginstance.
    # Check for taginstance
    $sql = qq[  select count(*) as count
                from assets 
                join asset_links on (assets.asset_id=asset_links.parent_id)
                where assets.asset_type='TagInstance'
                and asset_links.child_id=?];

    $id_ref = $self->{_dbh}->selectrow_hashref($sql, undef, ($id));
    
    if ($id_ref->{count}){
        # This is a tagged library.
        $self->is_tagged(1);

        # retrieve multiplex assets, tag sequence, tag id and tag group
        # (same tag id can be in multiple tag groups)
        $sql = qq[select distinct asset_id,tag_id,tag_group_id,sequence from  multiplex_information where library_asset_id=?];
 	my $sth = $self->{_dbh}->prepare($sql);
	$sth->execute($id);
        my @mplex_assets;
        while (my $mplex = $sth->fetchrow_hashref) {
            push @mplex_assets, $mplex->{'asset_id'};
            if ($self->tag_id){
                unless ($self->tag_id == $mplex->{'tag_id'}){
                    die "Library $id has multiple different tag ids";
                }
            }
            else {
                $self->tag_id($mplex->{'tag_id'});
                $self->tag_group_id($mplex->{'tag_group_id'});
                $self->tag_sequence($mplex->{'sequence'});
            }
        }
        $self->multiplex_pool_asset_ids(\@mplex_assets);
    }


    # populate _data for things like type and fragment_size from the library
    # creation request id in the property_information table:
    # get request id
    $sql = qq[select request_id from requests_new where target_asset_id=?];
    warn "LIB: $id";
    my $req_ref = $self->{_dbh}->selectrow_hashref($sql, undef, ($id));
    my $req_id = $req_ref->{'request_id'};
    if ($req_id){
        # get property information
        $sql = qq[select `key`, value from property_information where obj_type='request' and obj_id = ?];
        my $sth = $self->{_dbh}->prepare($sql);
        $sth->execute($req_id);
        my ($key, $value);
        $sth->bind_columns ( \$key, \$value);
        while ($sth->fetch) {
            $self->{'_data'}{$key} = $value;
        }
    }

    # if the above hasn't returned anything (i.e. for some pulldown libraries
    # ) fall back to getting the data from the asset_information table
    unless ($self->{'_data'}){
        $sql = qq[select param, value from asset_information where asset_id = ?];
        my $sth = $self->{_dbh}->prepare($sql);
        $sth->execute($id);
        my ($key, $value);
        $sth->bind_columns ( \$key, \$value);
        while ($sth->fetch) {
            $self->{'_data'}{$key} = $value;
        }
    }
    

    return $self;
}



=head2 seq_requests

  Arg [1]    : None
  Example    : my $seq_requests = $library->seq_requests();
  Description: Returns a ref to an array of the seq_request objects that are associated with this library.
  Returntype : ref to array of Sfind::Seq_Request objects

=cut

sub seq_requests {
    my ($self) = @_;
    unless ($self->{'requests'}){
	my @seq_requests;
    	foreach my $id (@{$self->seq_request_ids()}){
	    my $obj = $self->get_seq_request_by_id($id);
	    push @seq_requests, $obj if $obj->status;
	}
	$self->{'seq_requests'} = \@seq_requests;
    }

    return $self->{'seq_requests'};
}


=head2 seq_request_ids

  Arg [1]    : None
  Example    : my $seq_request_ids = $library->seq_request_ids();
  Description: Returns a ref to an array of the seq_request IDs that are associated with this library, including those on multiplexes that this library is in
  Returntype : ref to array of integer seq request IDs

=cut

sub seq_request_ids {
    my ($self) = @_;
    unless ($self->{'seq_request_ids'}){
        
        # seq requests can be on this library or a multiplex that this library
        # is in.
        my @asset_ids = ($self->id);

        if ($self->is_tagged){
            push @asset_ids,@{$self->multiplex_pool_asset_ids};
        }
            

        my $sql = qq[select distinct request_id 
                from 
                requests_new where asset_id=?
                and type like '%sequencing'];

        my @seq_requests;
        my $sth = $self->{_dbh}->prepare($sql);
            foreach my $asset_id(@asset_ids){
                $sth->execute($asset_id);
                foreach(@{$sth->fetchall_arrayref()}){
                    push @seq_requests, $_->[0];
                }
            }
        @seq_requests = sort {$a <=> $b} @seq_requests;

        $self->{'seq_request_ids'} = \@seq_requests;
    }
 
    return $self->{'seq_request_ids'};
}


=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $lib->id();
	       $lib->id('104');
  Description: Get/Set for ID of a library
  Returntype : SequenceScape ID (usu. integer)

=cut

sub id {
    my ($self,$id) = @_;
    if (defined $id and $id != $self->{'id'}){
	$self->{'id'} = $id;
    }
    return $self->{'id'};
}


=head2 sample_id

  Arg [1]    : sample_id (optional)
  Example    : my $sample_id = $lib->sample_id();
	       $lib->sample_id('104');
  Description: Get/Set for sample ID of a library
  Returntype : SequenceScape ID (usu. integer)

=cut

sub sample_id {
    my ($self,$sample_id) = @_;
    if (defined $sample_id and $sample_id != $self->{'sample_id'}){
	$self->{'sample_id'} = $sample_id;
    }
    unless ($self->{'sample_id'}){
	$self->_load_sample_study_id();
    }
    return $self->{'sample_id'};
}


=head2 study_id

  Arg [1]    : study_id (optional)
  Example    : my $study_id = $lib->study_id();
	       $lib->study_id('104');
  Description: Get/Set for study ID of a library
  Returntype : SequenceScape ID (usu. integer)

=cut

sub study_id {
    my ($self,$study_id) = @_;
    if (defined $study_id and $study_id ne $self->{'study_id'}){
	$self->{'study_id'} = $study_id;
    }
    unless ($self->{'study_id'}){
	$self->_load_sample_study_id();
    }
    return $self->{'study_id'};
}


=head2 name

  Arg [1]    : name (optional)
  Example    : my $name = $lib->name();
	       $lib->name('NA12006-CEU-1 1');
  Description: Get/Set for library name
  Returntype : string

=cut

sub name {
    my ($self,$name) = @_;
    if (defined $name and $name ne $self->{'name'}){
	$self->{'name'} = $name;
    }
    return $self->{'name'};
}


=head2 prep_status

  Arg [1]    : prep_status (optional)
  Example    : my $prep_status = $lib->prep_status();
	       $lib->prep_status('passed');
  Description: Get/Set for library prep_status
  Returntype : string

=cut

sub prep_status {
    my ($self,$prep_status) = @_;
    if (defined $prep_status and $prep_status ne $self->{'prep_status'}){
	$self->{'prep_status'} = $prep_status;
    }
    return $self->{'prep_status'};
}


=head2 get_seq_request_by_id

  Arg [1]    : request id from sequencescape
  Example    : my $seqrequest = $lib->get_seq_request_by_id(7447);
  Description: retrieve seq_request object by sequencescape id
  Returntype : Sfind::Seq_Request object

=cut

sub get_seq_request_by_id {
    my ($self, $id) = @_;
    my $obj = Sfind::Seq_Request->new($self->{_dbh},$id);
    return $obj;
}

=head2 get_lane_by_id

  Arg [1]    : lane id from sequencescape
  Example    : my $lane = $lib->get_lane_by_id(7447);
  Description: retrieve NPG lane object by id
  Returntype : Sfind::Lane object

=cut

sub get_lane_by_id {
    my ($self, $id) = @_;
    my $obj = Sfind::Lane->new($self->{_dbh},$id);
    return $obj;
}

=head2 fragment_size

  Arg [1]    : none
  Example    : my ($frag_from, $frag_to) = @{$lib->fragment_size};
  Description: Fetches requested library fragment size, as a pair of from and to sizes
  Returntype : arrayref of [from, to] sizes or undefs

=cut

sub fragment_size {
    my ($self, $id) = @_;
    my $from = lc($self->{'_data'}{'fragment_size_required_from'});
    my $to = lc($self->{'_data'}{'fragment_size_required_to'});
    # convert kb to 1000s
    foreach ($from,$to){
        if ($_ =~ /kb/){
            $_ =~ s/\s*kb//;
            $_ = $_ * 1000;
        }
    }
    ($from, $to) = ($to, $from) if $from > $to; # swap if necessary
    return [$from, $to];
}


=head2 type

  Arg [1]    : none
  Example    : my $lib_type = $lib->type();
  Description: returns type of library, if in database.  e.g. 'Standard', 'No PCR'
  Returntype : string or undef

=cut

sub type {
    my ($self, $id) = @_;
    return $self->{'_data'}{'library_type'};
}


=head2 multiplex_pool_asset_ids

  Arg [1]    : multiplex_tube_asset_ids
  Description: Get/Set for multiplex_tube_asset_ids
  Returntype : reference to a array

=cut

sub multiplex_pool_asset_ids{
    my ($self,$multiplex_pool_asset_ids) = @_;
    if (defined $multiplex_pool_asset_ids and $multiplex_pool_asset_ids ne $self->{'multiplex_pool_asset_ids'}){
	$self->{'multiplex_pool_asset_ids'} = $multiplex_pool_asset_ids;
    }
    return $self->{'multiplex_pool_asset_ids'} || [];
}

=head2 is_tagged

  Arg [1]    : boolean for whether library is tagged
  Example    : my $tag = $lib->is_tagged();
  Description: Get/Set for whether Library is tagged.  The tag is a short sequence tag added to each Library molecule so that it can be sequenced in a multiplex and the resulting mix resolved by the tag sequence.
  Returntype : boolean

=cut

sub is_tagged {
    my ($self,$is_tagged) = @_;
    if (defined $is_tagged and $is_tagged ne $self->{'is_tagged'}){
	$self->{'is_tagged'} = $is_tagged ? 1 : 0;
    }
    return $self->{'is_tagged'};
}


=head2 tag_id

  Arg [1]    : tag_id (optional)
  Example    : my $tag_id = $lib->tag_id();
	       $lib->tag_id('104');
  Description: Get/Set for tag ID of a library
  Returntype : int

=cut

sub tag_id {
    my ($self,$tag_id) = @_;
    if (defined $tag_id and $tag_id != $self->{'tag_id'}){
	$self->{'tag_id'} = $tag_id;
    }
    return $self->{'tag_id'};
}


=head2 tag_group_id

  Arg [1]    : tag_group_id (optional)
  Example    : my $tag_group_id = $lib->tag_group_id();
	       $lib->tag_group_id('104');
  Description: Get/Set for tag group ID.  Can have the same tag id in multiple tag groups.
  Returntype : int

=cut

sub tag_group_id {
    my ($self,$tag_group_id) = @_;
    if (defined $tag_group_id and $tag_group_id != $self->{'tag_group_id'}){
	$self->{'tag_group_id'} = $tag_group_id;
    }
    return $self->{'tag_group_id'};
}


=head2 tag_sequence

  Arg [1]    : tag_sequence (optional)
  Example    : my $tag_sequence = $lib->tag_sequence();
	       $lib->tag_sequence('ACTGATCGT');
  Description: Get/Set for tag sequence.
  Returntype : string

=cut

sub tag_sequence {
    my ($self,$tag_sequence) = @_;
    if (defined $tag_sequence and $tag_sequence ne $self->{'tag_sequence'}){
	$self->{'tag_sequence'} = $tag_sequence;
    }
    return $self->{'tag_sequence'};
}


# Internal function to populate sample_id and study_id
sub _load_sample_study_id {
    my ($self) = @_;

    my $sql = qq[select sample_id, study_id from requests_new where target_asset_id = ?];
    my $id_ref = $self->{_dbh}->selectrow_hashref($sql, undef, ($self->id));
    if($id_ref){
	$self->sample_id($id_ref->{sample_id});
	$self->study_id($id_ref->{study_id});
    }
}
1;
