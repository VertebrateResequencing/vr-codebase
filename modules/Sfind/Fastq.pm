package Sfind::Fastq; 
=head1 NAME

Sfind::Fastq - Sequence Tracking Fastq object

=head1 SYNOPSIS
    my $file = Sfind::Fastq->new($fastq_fuse_location);

    my $name = $lane->file();

=head1 DESCRIPTION

An object describing the tracked properties of a fastq file.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;
no warnings 'uninitialized';
use File::Basename;

=head2 new

  Arg [1]    : fastq fuse location (e.g. /fuse/mpsafs/runs/1378/1378_s_1.fastq
  )
  Example    : my $file= Sfind::Fastq->new($file_loc)
  Description: Returns Fastq object by file location.  Assumes a fastqcheck file exists alongside the fastq.
  Returntype : Sfind::Fastq object

=cut

sub new {
    my ($class,$file_loc) = @_;
    die "Need to call with a file location" unless $file_loc;
    die "No such file: $file_loc" unless -f $file_loc;
    my $self = {};
    bless ($self, $class);
    
    my($filename, $path, $suffix) = fileparse($file_loc, qr/\.[^.]*/);
    my $fastqcheck = "$path/$filename.fastqcheck";
    if (-f $fastqcheck){
        $self->fastqcheck_location($fastqcheck);
    }
    else {
        die "No fastqcheck file $fastqcheck";
    }
    $self->location($file_loc);
    $self->name($filename.$suffix);
    return $self;
}


=head2 name

  Arg [1]    : None
  Example    : my $file_name = $fastq->name();
  Description: Fastq file name. e.g 1378_s_1.fastq
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
  Example    : my $file_location = $fastq->location();
  Description: Fastq file location, e.g. /fuse/mpsafs/runs/1378/1378_s_1.fastq

  Returntype : string

=cut

sub location {
    my ($self,$location) = @_;
    if (defined $location and $location ne $self->{'location'}){
	$self->{'location'} = $location;
    }
    return $self->{'location'};
}


=head2 fastqcheck_location

  Arg [1]    : None
  Example    : my $file_fastqcheck_location = $fastq->fastqcheck_location();
  Description: Fastq file fastqcheck_location, e.g. /fuse/mpsafs/runs/1378/1378_s_1.fastqcheck

  Returntype : string

=cut

sub fastqcheck_location {
    my ($self,$fastqcheck_location) = @_;
    if (defined $fastqcheck_location and $fastqcheck_location ne $self->{'fastqcheck_location'}){
	$self->{'fastqcheck_location'} = $fastqcheck_location;
    }
    return $self->{'fastqcheck_location'};
}


=head2 reads

  Arg [1]    : none
  Example    : my $reads = $fastq->reads();
  Description: get total number of reads in resulting fastq
  Returntype : integer count of reads

=cut

sub reads {
    my ($self) = @_;
    unless ($self->{'reads'}){
	$self->_get_fastqcheck_stats;
    }
    return $self->{'reads'};
}

=head2 mean_quality

  Arg [1]    : none
  Example    : my $mean_quality = $fastq->mean_quality();
  Description: get mean_quality (from fastqcheck) in resulting fastq
  Returntype : mean_quality to 1dp

=cut

sub mean_quality {
    my ($self) = @_;
    unless ($self->{'mean_quality'}){
	$self->_get_fastqcheck_stats;
    }
    return $self->{'mean_quality'};
}


=head2 basepairs

  Arg [1]    : none
  Example    : my $basepairs = $fastq->basepairs();
  Description: get total number of basepairs in resulting fastq
  Returntype : integer count of basepairs

=cut

sub basepairs {
    my ($self) = @_;
    unless ($self->{'basepairs'}){
	$self->_get_fastqcheck_stats;
    }
    return $self->{'basepairs'};
}


=head2 read_len

  Arg [1]    : none
  Example    : my $read_len = $fastq->read_len();
  Description: Get fastq read_len, i.e. number of cycles
  Returntype : integer

=cut

sub read_len {
    my ($self) = @_;
    unless ($self->{'read_len'}){
	$self->_get_fastqcheck_stats;
    }
    return $self->{'read_len'};
}


=head2 md5

  Arg [1]    : none
  Example    : my $md5 = $fastq->md5();
  Description: Get fastq md5 (uses mpsa_download)
  Returntype : string

=cut

sub md5 {
    my ($self) = @_;
    unless ($self->{'md5'}){
	$self->{'md5'} = $self->_get_mpsa_md5();
    }
    return $self->{'md5'};
}


# Internal function to populate basepairs, reads, and mean_quality from
# fastqcheck
sub _get_fastqcheck_stats {
    my ($self) = @_;
    my $fqc = $self->fastqcheck_location;

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

# Internal function to get md5 for an mpsa file
sub _get_mpsa_md5 {
    my ($self) = @_;
    my $name = $self->name;
    open(my $MD5, "/software/solexa/bin/mpsa_download -m -f $name|") or die "Can't run mpsa_download to get md5: $!\n";
    my $line = <$MD5>;
    chomp $line;
    my ($md5, $file) = split / +/, $line;
    return $md5;
}
1;
