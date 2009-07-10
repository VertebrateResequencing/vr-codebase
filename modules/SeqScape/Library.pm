package SeqScape::Library;
=head1 NAME

SeqScape::Library - object for SequenceScape Library web API

=head1 SYNOPSIS
    my $lib = SeqScape::Library->new($lib_id);

    #get arrayref of request objects in this library
    my $reqs = $lib->requests();
    
    my $id = $lib->id();
    my $name = $lib->name();

=head1 DESCRIPTION

An object to fetch/contain/parse the XML from the SequenceScape web API that
describes a library.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut


use strict;
use SeqScape::XMLObj;
use SeqScape::Request;
use Data::Dumper;
our @ISA = qw(SeqScape::XMLObj);
our $PE_SEQ = "paired end sequencing";
our $LIB_CREATE = "library creation";


=head2 url_from_id

  Arg [1]    : Library id, e.g. 140
  Example    : my $url = $library->url_from_id($id);
  Description: Returns a full URL for a library xml request given a library id
  Returntype : string

=cut

sub url_from_id {
    my ($self, $id) = @_;
    my $root_url = $self->{'root_url'};
    my $url = "$root_url/items/$id.xml";
    return $url;
}


=head2 status

  Arg [1]    : None
  Example    : my $status = $library->status();
  Description: Returns the status description of this library, e.g. "passed"
  Exceptions : Dies if can't find status
  Returntype : string

=cut

sub status {
    my ($self) = @_;
    unless ($self->{'status'}){
	my $state;
	foreach my $req (@{$self->requests}){
	    next unless $req->type eq $LIB_CREATE;
	    $state = $req->status;
	}
	#die "Can't find status" unless $state;
	$state ||= "unknown";
	$self->{'status'} = $state;
    }

    return $self->{'status'};

}


=head2 requests

  Arg [1]    : None
  Example    : my $requests = $library->requests();
  Description: Returns a ref to an array of the request objects that are associated with this library
  Returntype : ref to array of SeqScape::Request objects

=cut

sub requests {
    my ($self) = @_;

    unless ($self->{'requests'}){
	my @requests;
    	foreach my $req_id (@{$self->request_ids()}){
	    my $req = SeqScape::Request->new($req_id);
	    push @requests, $req;
	}
	$self->{'requests'} = \@requests;
    }

    return $self->{'requests'};
}


=head2 seq_requests

  Arg [1]    : None
  Example    : my $seq_requests = $library->seq_requests();
  Description: Returns a ref to an array of the sequencing request objects that are associated with this library
  Returntype : ref to array of SeqScape::Request objects

=cut

sub seq_requests {
    my ($self) = @_;
    my @seq_req = grep {$_->type eq $PE_SEQ} @{$self->requests()};
    return \@seq_req;
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
	my $xml = $self->{xml};
	my @requests;
	foreach my $request(@{$xml->{requests}[0]{'request'}}){
	    push @requests, $request->{'id'}
	}
	$self->{'request_ids'} = \@requests;
    }

    return $self->{'request_ids'};
}
1;
