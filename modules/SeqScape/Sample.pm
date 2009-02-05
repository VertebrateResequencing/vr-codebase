package SeqScape::Sample;
=head1 NAME

SeqScape::Sample - object for SequenceScape Sample web API

=head1 SYNOPSIS
    my $proj = SeqScape::Sample->new($sample_id);

    #get arrayref of library objects in a sample
    my $libs = $sample->libraries();
    
    my $id = $sample->id();
    my $name = $sample->name();

=head1 DESCRIPTION

An object to fetch/contain/parse the XML from the SequenceScape web API that
describes a sample.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut


use SeqScape::XMLObj;
use SeqScape::Library;
our @ISA = qw(SeqScape::XMLObj);


=head2 url_from_id

  Arg [1]    : Sample id, e.g. 140
  Example    : my $url = $sample->url_from_id($id);
  Description: Returns a full URL for a sample xml request given a sample id
  Returntype : string

=cut

sub url_from_id {
    my ($self, $id) = @_;
    my $root_url = $self->{'root_url'};
    my $url = "$root_url/samples/$id.xml";
    return $url;
}


=head2 libraries

  Arg [1]    : None
  Example    : my $libraries = $sample->libraries();
  Description: Returns a ref to an array of the library objects that are associated with this sample
  Returntype : ref to array of SeqScape::Library objects

=cut

sub libraries {
    my ($self) = @_;

    unless ($self->{'libraries'}){
	my @libraries;
    	foreach my $id (@{$self->library_ids()}){
	    my $obj = SeqScape::Library->new($id);
	    push @libraries, $obj;
	}
	$self->{'libraries'} = \@libraries;
    }

    return $self->{'libraries'};
}


=head2 library_ids

  Arg [1]    : None
  Example    : my $library_ids = $sample->library_ids();
  Description: Returns a ref to an array of the library IDs that are associated with this sample
  Returntype : ref to array of integer library IDs

=cut

sub library_ids {
    my ($self) = @_;

    unless ($self->{'library_ids'}){
	my @libs;

	foreach my $lib (@{$self->{xml}{'items'}[0]{'item'}}){
	    push @libs, $lib->{'id'};
	}

	$self->{'library_ids'} = \@libs;
    }

    return $self->{'library_ids'};
}

1;
