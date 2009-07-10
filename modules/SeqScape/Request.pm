package SeqScape::Request;
=head1 NAME

SeqScape::Request - object for SequenceScape Request web API

=head1 SYNOPSIS
    my $proj = SeqScape::Request->new($request_id);

    my $id	= $request->id();
    my $name	= $request->name();
    my $status	= $request->status();

=head1 DESCRIPTION

An object to fetch/contain/parse the XML from the SequenceScape web API that
describes a request.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut

use strict;
use SeqScape::XMLObj;
our @ISA = qw(SeqScape::XMLObj);


=head2 url_from_id

  Arg [1]    : Request id, e.g. 140
  Example    : my $url = $request->url_from_id($id);
  Description: Returns a full URL for a request xml request given a request id
  Returntype : string

=cut

sub url_from_id {
    my ($self, $id) = @_;
    my $root_url = $self->{'root_url'};
    my $url = "$root_url/requests/$id.xml";
    return $url;
}


=head2 type

  Arg [1]    : None
  Example    : my $type = $request->type();
  Description: Returns the type of request, e.g. "paired end sequencing"
  Returntype : lowercase string
  Note	     : This is going to change in the web api to something more robust

=cut

sub type {
    my ($self) = @_;
    return lc($self->{xml}{'template'}[0]{content}); 
}


=head2 status

  Arg [1]    : None
  Example    : my $status = $request->status();
  Description: Returns a status description (e.g. 'passed') for this request
  Returntype : string

=cut

sub status {
    my ($self) = @_;
    return $self->{xml}{'state'}[0];
}


=head2 read_length

  Arg [1]    : None
  Example    : my $len = $request->read_length();
  Description: Returns integer read length, e.g. 37 of one end
  Returntype : integer

=cut

sub read_length {
    my ($self) = @_;
    return $self->{xml}{'read_length'}[0];
}


=head2 seq_run

  Arg [1]    : None
  Example    : my $npg_run = $request->seq_run();
  Description: when a sequencing request is complete, the NPG id of the run
	       should be available.  seq_run gives the forward run id, as used
	       by dfind to locate the run output.
	       Returns undef if the run is incomplete or the id isn't available.
  Returntype : integer

=cut

sub seq_run {
    my ($self) = @_;
    my $descs = $self->{xml}{'descriptors'}[0]{'descriptor'};
    my $fwdrun;
    foreach (@$descs){
	if ($_->{name}[0] eq 'Forward run ID'){
	    $fwdrun = $_->{value}[0];
	}
    }
    return $fwdrun;
}


=head2 seq_lane

  Arg [1]    : None
  Example    : my $npg_lane = $request->seq_lane();
  Description: when a sequencing request is complete, the lane number of the run
	       should be available.  See seq_lane for the run id.
	       Returns undef if the run is incomplete or the id isn't available.
  Returntype : integer

=cut

sub seq_lane {
    my ($self) = @_;
    my $descs = $self->{xml}{'descriptors'}[0]{'descriptor'};
    my $lane;
    foreach (@$descs){
	if ($_->{name}[0] eq 'Lane'){
	    $lane = $_->{value}[0];
	}
    }
    return $lane;
}


=head2 last_updated

  Arg [1]    : None
  Example    : my $ss_update = $request->last_updated();
  Description: Gives the last updated date from SequenceScape
  Returntype : Date string

=cut

sub last_updated {
    my ($self) = @_;
    return $self->{xml}{'updated_at'}[0];
}

1;

