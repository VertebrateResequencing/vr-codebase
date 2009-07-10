package SeqScape::Project;
=head1 NAME

SeqScape::Project - object for SequenceScape Project web API

=head1 SYNOPSIS
    my $proj = SeqScape::Project->new($project_id);

    #get arrayref of sample objects in a project
    my $samples = $project->samples();
    
    my $id = $project->id();
    my $name = $project->name();

=head1 DESCRIPTION

An object to fetch/contain/parse the XML from the SequenceScape web API that
describes a project.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut

use strict;
use SeqScape::XMLObj;
use SeqScape::Sample;
our @ISA = qw(SeqScape::XMLObj);


=head2 url_from_id

  Arg [1]    : Project id, e.g. 140
  Example    : my $url = $project->url_from_id($id);
  Description: Returns a full URL for a project xml request given a project id
  Returntype : string

=cut

sub url_from_id {
    my ($self, $id) = @_;
    my $root_url = $self->{'root_url'};
    my $url = "$root_url/projects/$id.xml";
    return $url;
}

=head2 samples

  Arg [1]    : None
  Example    : my $samples = $project->samples();
  Description: Returns a ref to an array of the sample objects that are associated with this project
  Returntype : ref to array of SeqScape::Sample objects

=cut

sub samples {
    my ($self) = @_;

    unless ($self->{'samples'}){
	my @samples;
    	foreach my $id (@{$self->sample_ids()}){
	    my $obj = SeqScape::Sample->new($id);
	    push @samples, $obj;
	}
	$self->{'samples'} = \@samples;
    }

    return $self->{'samples'};
}


=head2 sample_ids

  Arg [1]    : None
  Example    : my $sample_ids = $project->sample_ids();
  Description: Returns a ref to an array of the sample IDs that are associated with this project
  Returntype : ref to array of integer sample IDs

=cut

sub sample_ids {
    my ($self) = @_;

    unless ($self->{'sample_ids'}){
	my $root_url = $self->{'root_url'};
	my $id = $self->id;
	my $url = "$root_url/projects/samples/$id";
	my $xml = $self->xml_data($url);
	my @samples;
	foreach my $sample(@{$xml->{sample}}){
	    push @samples, $sample->{'id'}[0]{content};
	}
	$self->{'sample_ids'} = \@samples;
    }

    return $self->{'sample_ids'};
}

1;
