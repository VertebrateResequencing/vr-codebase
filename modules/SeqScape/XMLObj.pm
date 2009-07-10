package SeqScape::XMLObj;
=head1 NAME

SeqScape::XMLObj - general XML object for SequenceScape web API

=head1 SYNOPSIS
    # Not used directly, is subclassed by SeqScape::Project, Sample, etc

=head1 DESCRIPTION

An object to fetch/contain/parse the XML from the SequenceScape web API

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut


use strict;
use warnings;
#no warnings 'uninitialized';
use XML::Simple;
use WWW::Mechanize;
use Data::Dumper;

=head2 new

  Arg [1]    : ID of object to retrieve from web API
  Example    : my $obj = $obj->new($obj_id)
  Description: Returns object requested.
  Returntype : object requested

=cut

sub new {
    my ($class, $id) = @_;

    die "Need to call with an ID" unless $id;

    my $self = {};
    bless ($self, $class);

    $self->{'root_url'} = qq(http://psd-production.internal.sanger.ac.uk:6600);
    $self->{'id'} = $id;

    my $url = $self->url_from_id($id);

    my $xml = $self->xml_data($url);

    $self->{xml} = $xml;

    return $self;
}

=head2 id

  Arg [1]    : None
  Example    : my $id = $obj->id();
  Description: Returns SequenceScape ID of object
  Returntype : SequenceScape ID (usu. integer)

=cut

sub id {
    my $self = shift;
    return $self->{'id'};
}


=head2 xml_data

  Arg [1]    : URL to XML resource
  Example    : my $xml = $obj->xml_data($url);
  Description: Fetches the XML from a SequenceScape (or other) URL
  Returntype : XML::Simple parsed XML string

=cut

sub xml_data {
  my ($self, $url) = @_;
  my $browser = WWW::Mechanize->new();
  $browser->add_header('Accept' => 'application/xml');

  my $obj;
  eval{
    # download the json page:
    #print "Getting xml $url\n";
    $browser->get( $url );
    my $content = $browser->content();
    my $xml = XML::Simple->new(KeyAttr => [], ForceArray=>1);
    $obj = $xml->XMLin($content);
  };
  if($@){
    die "[[XML ERROR]] $url: $@\n";
  }
  $obj->{url} = $url;
  return $obj;
}


=head2 name

  Arg [1]    : None
  Example    : my $name = $obj->name();
  Description: Returns the SequenceScape name of this object
  Returntype : string

=cut

sub name {
    my ($self) = @_;
    return $self->{xml}{name}[0];
}

1;
