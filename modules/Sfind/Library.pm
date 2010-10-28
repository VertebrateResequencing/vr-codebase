package Sfind::Library; 
=head1 NAME

Sfind::Library - Sequence Tracking Library object

=head1 SYNOPSIS
    my $lib = Sfind::Library->new($dbh, $library_id);

    #get arrayref of requests on a library
    my $reqs = $library->requests();

    #get arrayref of lanes from a library
    my $lanes = $library->lanes();
    
    Note that while the hierarchy should go library->request->lane
    the link between request and lane is often broken, so the two
    sets can't always be correlated.  This is why this object also
    provides direct access to lanes.
    
    my $id = $library->id();
    my $name = $library->name();
    my $sequenced_bp = $library->sequenced_bases();

=head1 DESCRIPTION

An object describing the tracked properties of a library.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;
no warnings 'uninitialized';
use Sfind::Request;
use Sfind::Lane;

=head2 new

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : library id
  Example    : my $lib = Sfind::Sfind->new($dbh, $id)
  Description: Returns Library object by library_id
  Returntype : Sfind::Library object

=cut

sub new {
    my ($class,$dbh, $id) = @_;
    die "Need to call with a db handle and id" unless ($dbh && $id);
    my $self = {};
    bless ($self, $class);
    $self->{_dbh} = $dbh;

    # get item_information values
    my $sql = qq[select item_name, param, value from item_information where item_id = ?];
    my $sth = $self->{_dbh}->prepare($sql);
    $sth->execute($id);
    my ($name, $param, $value);
    $sth->bind_columns (\$name, \$param, \$value);
    while ($sth->fetch) {
	$self->{'_data'}{$param} = $value;
    }

    $name =~ s/\s+$//;  # trim trailing whitespace
    $self->id($id);
    # 2010-10-22 problem with warehouse - item_information is missing.
    # So, let's put in a hack to at least create a library if there is no name
    # hopefully it will get auto-fixed when the item_information comes back
    $name ||= $id;
    $self->name($name);

    # get library creation status
    $sql = qq[select state from requests where item_id=? and type = 'Library creation'];
    my $id_ref = $self->{_dbh}->selectrow_hashref($sql, undef, ($id));
    $self->prep_status($id_ref->{state});
    return $self;
}



=head2 requests

  Arg [1]    : None
  Example    : my $requests = $library->requests();
  Description: Returns a ref to an array of the request objects that are associated with this library.
  Returntype : ref to array of Sfind::Request objects

=cut

sub requests {
    my ($self) = @_;
    unless ($self->{'requests'}){
	my @requests;
    	foreach my $id (@{$self->request_ids()}){
	    my $obj = $self->get_request_by_id($id);
	    push @requests, $obj if $obj->status;
	}
	$self->{'requests'} = \@requests;
    }

    return $self->{'requests'};
}


=head2 request_ids

  Arg [1]    : None
  Example    : my $request_ids = $library->request_ids();
  Description: Returns a ref to an array of the request IDs that are associated with this library
  Returntype : ref to array of integer request IDs

=cut

sub request_ids {
    my ($self) = @_;
    unless ($self->{'request_ids'}){
	my $sql = qq[select distinct(request_id) from requests where item_id=? and type in ('single ended sequencing','paired end sequencing')];
	my @requests;
	my $sth = $self->{_dbh}->prepare($sql);

	$sth->execute($self->id);
	foreach(@{$sth->fetchall_arrayref()}){
	    push @requests, $_->[0];
	}
	@requests = sort {$a <=> $b} @requests;

	$self->{'request_ids'} = \@requests;
    }
 
    return $self->{'request_ids'};
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
    if (defined $id and $id ne $self->{'id'}){
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
    if (defined $sample_id and $sample_id ne $self->{'sample_id'}){
	$self->{'sample_id'} = $sample_id;
    }
    unless ($self->{'sample_id'}){
	$self->_load_sample_project_id();
    }
    return $self->{'sample_id'};
}


=head2 project_id

  Arg [1]    : project_id (optional)
  Example    : my $project_id = $lib->project_id();
	       $lib->project_id('104');
  Description: Get/Set for project ID of a library
  Returntype : SequenceScape ID (usu. integer)

=cut

sub project_id {
    my ($self,$project_id) = @_;
    if (defined $project_id and $project_id ne $self->{'project_id'}){
	$self->{'project_id'} = $project_id;
    }
    unless ($self->{'project_id'}){
	$self->_load_sample_project_id();
    }
    return $self->{'project_id'};
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


=head2 get_request_by_id

  Arg [1]    : request id from sequencescape
  Example    : my $request = $lib->get_request_by_id(7447);
  Description: retrieve request object by sequencescape id
  Returntype : Sfind::Request object

=cut

sub get_request_by_id {
    my ($self, $id) = @_;
    my $obj = Sfind::Request->new($self->{_dbh},$id);
    return $obj;
}


=head2 lanes

  Arg [1]    : None
  Example    : my $lanes = $library->lanes();
  Description: Returns a ref to an array of the lanes objects that have been generated from this library.
  Returntype : ref to array of Sfind::Lane objects

=cut

sub lanes {
    my ($self) = @_;
    unless ($self->{'lanes'}){
	my @lanes;
    	foreach my $id (@{$self->lane_ids()}){
	    push @lanes, $self->get_lane_by_id($id);
	}
	$self->{'lanes'} = \@lanes;
    }
    return $self->{'lanes'};
}


=head2 lane_ids

  Arg [1]    : None
  Example    : my $lane_ids = $library->lane_ids();
  Description: Returns a ref to an array of the lane IDs that are associated with this library.
  Returntype : ref to array of integer lane IDs

=cut

sub lane_ids {
    my ($self) = @_;
    unless ($self->{'lane_ids'}){
	my $sql = qq[select id_npg_information from npg_information n, library l where l.item_id=? and l.batch_id =n.batch_id and l.position = n.position and (n.id_run_pair=0 or n.id_run_pair is null);];
	my @lanes;
	my $sth = $self->{_dbh}->prepare($sql);

	$sth->execute($self->id);
	foreach(@{$sth->fetchall_arrayref()}){
	    push @lanes, $_->[0];
	}
	@lanes = sort {$a <=> $b} @lanes;

	$self->{'lane_ids'} = \@lanes;
    }
 
    return $self->{'lane_ids'};
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


=head2 sequenced_bases

  Arg [1]    : none
  Example    : my $tot_bp = $lib->sequenced_bases();
  Description: the total number of sequenced bases on this library.
		This is the sum of the bases from the fastq files associated
		with this library in NPG.
  Returntype : integer

=cut

sub sequenced_bases {
    my ($self, $id) = @_;
    unless ($self->{'seq_bases'}){
	$self->{'seq_bases'} = 0;
	foreach my $lane(@{$self->lanes}){
	    $self->{'seq_bases'} += $lane->basepairs;
	}
    }
    return $self->{'seq_bases'};
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


# Internal function to populate sample_id and project_id
sub _load_sample_project_id {
    my ($self) = @_;

    my $sql = qq[select sample_id, project_id from requests where item_id = ? and type='Library creation'];
    my $id_ref = $self->{_dbh}->selectrow_hashref($sql, undef, ($self->id));
    if($id_ref){
	$self->sample_id($id_ref->{sample_id});
	$self->project_id($id_ref->{project_id});
    }
}
1;
