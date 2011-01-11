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

=cut

use strict;
use warnings;
no warnings 'uninitialized';
use File::Basename;
use VertRes::Wrapper::iRODS; # for finding bam files

=head2 new

  Arg [1]    : bam irods location (e.g. /seq/5322/5322_8.bam)
  Example    : my $file= Sfind::BAM->new($irods_loc)
  Description: Returns BAM object by irods location.
  Returntype : Sfind::BAM object

=cut

sub new {
    my ($class,$file_loc) = @_;
    die "Need to call with a file location" unless $file_loc;
    # TODO change to check irods for file
    # die "No such file: $file_loc" unless -f $file_loc;
    my $self = {};
    bless ($self, $class);
    
    my($filename, $path, $suffix) = fileparse($file_loc, qr/\.[^.]*/);
    #my $bamcheck = "$path/$filename.bamcheck";
    #if (-f $bamcheck){
    #    $self->bamcheck_location($bamcheck);
    #}
    #else {
    #    die "No bamcheck file $bamcheck";
    #}
    $self->location($file_loc);
    $self->name($filename.$suffix);
    return $self;
}


=head2 name

  Arg [1]    : None
  Example    : my $file_name = $bam->name();
  Description: BAM file name. e.g 5322_8.bam
  Returntype : string

=cut

sub name {
    my ($self,$name) = @_;
    if (defined $name and $name ne $self->{'name'}){
	$self->{'name'} = $name;
    }
    return $self->{'name'};
}


=head2 location

  Arg [1]    : None
  Example    : my $file_location = $bam->location();
  Description: BAM irods location, e.g. /seq/5322/5322_8.bam

  Returntype : string

=cut

sub location {
    my ($self,$location) = @_;
    if (defined $location and $location ne $self->{'location'}){
	$self->{'location'} = $location;
    }
    return $self->{'location'};
}


=head2 bamcheck_location

  Arg [1]    : None
  Example    : my $file_bamcheck_location = $bam->bamcheck_location();
  Description: BAM file bamcheck_location, e.g. /seq/5322/5322_8.bam)

  Returntype : string

=cut

sub bamcheck_location {
    my ($self,$bamcheck_location) = @_;
    if (defined $bamcheck_location and $bamcheck_location ne $self->{'bamcheck_location'}){
	$self->{'bamcheck_location'} = $bamcheck_location;
    }
    return $self->{'bamcheck_location'};
}


=head2 reads

  Arg [1]    : none
  Example    : my $reads = $bam->reads();
  Description: get total number of reads in resulting bam
  Returntype : integer count of reads

=cut

sub reads {
    my ($self) = @_;
    unless ($self->{'reads'}){
	$self->_get_bamcheck_stats;
    }
    return $self->{'reads'};
}

=head2 mean_quality

  Arg [1]    : none
  Example    : my $mean_quality = $bam->mean_quality();
  Description: get mean_quality (from bamcheck) in resulting bam
  Returntype : mean_quality to 1dp

=cut

sub mean_quality {
    my ($self) = @_;
    unless ($self->{'mean_quality'}){
	$self->_get_bamcheck_stats;
    }
    return $self->{'mean_quality'};
}


=head2 basepairs

  Arg [1]    : none
  Example    : my $basepairs = $bam->basepairs();
  Description: get total number of basepairs in resulting bam
  Returntype : integer count of basepairs

=cut

sub basepairs {
    my ($self) = @_;
    unless ($self->{'basepairs'}){
	$self->_get_bamcheck_stats;
    }
    return $self->{'basepairs'};
}


=head2 read_len

  Arg [1]    : none
  Example    : my $read_len = $bam->read_len();
  Description: Get bam read_len, i.e. number of cycles
  Returntype : integer

=cut

sub read_len {
    my ($self) = @_;
    unless ($self->{'read_len'}){
	$self->_get_bamcheck_stats;
    }
    return $self->{'read_len'};
}


=head2 md5

  Arg [1]    : none
  Example    : my $md5 = $bam->md5();
  Description: Get bam md5 from irods
  Returntype : string

=cut

sub md5 {
    my ($self) = @_;
    unless ($self->{'md5'}){
        my $irods = VertRes::Wrapper::iRODS->new();
	$self->{'md5'} = $irods->get_file_md5($self->location);
    }
    return $self->{'md5'};
}


# Internal function to populate basepairs, reads, and mean_quality from
# bamcheck
sub _get_bamcheck_stats {
    my ($self) = @_;
    my $fqc = $self->bamcheck_location;
    unless ($fqc){
        # Currently don't have a way to make bamcheck, and they aren't in irods
        # so bail out undef.
        return undef;
    }

    my $read_tot = 0;
    my $bp_tot = 0;
    my $mean_q = 0;
    my $readlen = 0;
    my $q_tot = 0;
    my $q_n = 0;

    open(my $FQC, $fqc) or die "Can't open $fqc to retrieve readcount: $!\n";
    my $header = <$FQC>;
    chomp $header;
    my $totals;
    while(<$FQC>){
        if (/^Total/){
            chomp;
            $totals = $_;
            last;
        }
    }
    close $FQC;

    my ($reads, $seqfoo, $bp, $tfoo,$lfoo,$cycles,$foo) = split " ", $header;
    $read_tot += $reads;
    $bp_tot += $bp;
    $readlen = $cycles;

    my ($label, $atot, $ctot, $gtot, $ttot, $ntot, @vals) = split " ", $totals;	
    my $AQ = pop @vals; # NB, this is NOT average quality
    # the numbers in @vals are the counts (in millions I think) of bases
    # at the quality of the subscript (i.e. from 0..x)
    # want average of these
    for (my $i=0; $i < scalar @vals; ++$i){
        $q_n += $vals[$i]; 
        $q_tot += $vals[$i] * $i;
    }

    if ($q_n){
	$mean_q = $q_tot/$q_n;
    }

    $self->{'reads'} = $read_tot;
    $self->{'basepairs'} = $bp_tot;
    $self->{'read_len'} = $readlen;
    $self->{'mean_quality'} = sprintf("%.1f",$mean_q);
}

1;
