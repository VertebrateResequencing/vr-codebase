package Sfind::BAM; 
=head1 NAME

Sfind::BAM - Sequence Tracking BAM object

=head1 SYNOPSIS
    my $file = Sfind::BAM->new($bam_irods_location);

    my $name = $file->name();
    my $readlen = $file->readlen();
    my $md5 = $file->md5();
    my $basepairs = $file->basepairs();
    my $reads = $file->reads();

=head1 DESCRIPTION

An object describing the tracked properties of a bam file.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=head2 new

  Arg [1]    : bam irods location (e.g. /seq/5322/5322_8.bam)
  Example    : my $file= Sfind::BAM->new($irods_loc)
  Description: Returns BAM object by irods location.
  Returntype : Sfind::BAM object

=head2 location

  Arg [1]    : None
  Example    : my $file_location = $bam->location();
  Description: BAM irods location, e.g. /seq/5322/5322_8.bam

  Returntype : string

=head2 name

  Arg [1]    : None
  Example    : my $file_name = $bam->name();
  Description: BAM file name. e.g 5322_8.bam
  Returntype : string

=head2 md5

  Arg [1]    : none
  Example    : my $md5 = $bam->md5();
  Description: Get bam md5 from irods
  Returntype : string

=cut

use Moose;
use namespace::autoclean;
use File::Basename;
use VertRes::Wrapper::iRODS; # for finding bam files

has 'location'=> (
    is          => 'ro',
    isa         => 'Str',
    required    => 1,
);

has 'name'=> (
    is          => 'ro',
    isa         => 'Str',
    required    => 1,
);

has 'md5'=> (
    is          => 'ro',
    isa         => 'Str',
    lazy        => 1,
    builder     => '_get_md5',
);

around BUILDARGS => sub {
    my $orig  = shift;
    my $class = shift;
    
    my $argref;
    if ( @_ == 1 && !ref $_[0] ) {
          $argref = $class->$orig( location => $_[0] );
    }
    else {
          $argref = $class->$orig(@_);
    }

    my($filename, $path, $suffix) = fileparse($argref->{location}, qr/\.[^.]*/);
    $argref->{name} = $filename.$suffix;

    return $argref;
};
    

###############################################################################
# BUILDERS
###############################################################################

sub _get_md5 {
    my ($self) = @_;
    my $irods = VertRes::Wrapper::iRODS->new();
    my $md5 = $irods->get_file_md5($self->location);
    return $md5;
}

# TODO - add metadata methods
1;
