package Sfind::Fastq; 
=head1 NAME

Sfind::Fastq - Sequence Tracking Fastq object

=head1 SYNOPSIS
    my $file = Sfind::Fastq->new($fastq_fuse_location);

    my $readlen = $file->readlen();
    my $md5 = $file->md5();
    my $basepairs = $file->basepairs();
    my $reads = $file->reads();
    my $name = $file->name();

=head1 DESCRIPTION

An object describing the tracked properties of a fastq file.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=head2 new

  Arg [1]    : fastq fuse location (e.g. /fuse/mpsafs/runs/1378/1378_s_1.fastq
  )
  Example    : my $file= Sfind::Fastq->new($file_loc)
  Description: Returns Fastq object by file location.  Assumes a fastqcheck file exists alongside the fastq.
  Returntype : Sfind::Fastq object


=head2 name

  Arg [1]    : None
  Example    : my $file_name = $fastq->name();
  Description: Fastq file name. e.g 1378_s_1.fastq
  Returntype : string


=head2 location

  Arg [1]    : None
  Example    : my $file_location = $fastq->location();
  Description: Fastq file location, e.g. /fuse/mpsafs/runs/1378/1378_s_1.fastq
  Returntype : string


=head2 fastqcheck_location

  Arg [1]    : None
  Example    : my $file_fastqcheck_location = $fastq->fastqcheck_location();
  Description: Fastq file fastqcheck_location, e.g. /fuse/mpsafs/runs/1378/1378_s_1.fastqcheck
  Returntype : string


=head2 md5

  Arg [1]    : none
  Example    : my $md5 = $fastq->md5();
  Description: Get fastq md5 (uses mpsa_download)
  Returntype : string


=head2 reads

  Arg [1]    : none
  Example    : my $reads = $fastq->reads();
  Description: get total number of reads in resulting fastq
  Returntype : integer count of reads


=head2 mean_quality

  Arg [1]    : none
  Example    : my $mean_quality = $fastq->mean_quality();
  Description: get mean_quality (from fastqcheck) in resulting fastq
  Returntype : mean_quality to 1dp


=head2 basepairs

  Arg [1]    : none
  Example    : my $basepairs = $fastq->basepairs();
  Description: get total number of basepairs in resulting fastq
  Returntype : integer count of basepairs


=head2 read_len

  Arg [1]    : none
  Example    : my $read_len = $fastq->read_len();
  Description: Get fastq read_len, i.e. number of cycles
  Returntype : integer


=cut

use Moose;
use namespace::autoclean;
use File::Basename;

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

has 'fastqcheck_location'=> (
    is          => 'ro',
    isa         => 'Str',
    required    => 0,
);

has 'md5'=> (
    is          => 'ro',
    isa         => 'Str',
    lazy        => 1,
    builder     => '_get_md5',
);

has 'reads' => (
    is          => 'ro',
    isa         => 'Int',
    lazy        => 1,
    builder     => '_get_reads',
);

has 'basepairs' => (
    is          => 'ro',
    isa         => 'Int',
    lazy        => 1,
    builder     => '_get_basepairs',
);

has 'mean_quality' => (
    is          => 'ro',
    isa         => 'Num',
    lazy        => 1,
    builder     => '_get_mean_quality',
);

has 'read_len' => (
    is          => 'ro',
    isa         => 'Int',
    lazy        => 1,
    builder     => '_get_read_len',
);

has '_fastqcheck_stats' => (
    is          => 'ro',
    isa         => 'HashRef',
    lazy        => 1,
    builder     => '_get_fastqcheck_stats',
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
    my $fastqcheck = "$path/$filename.fastqcheck";
    # there may not be fastqcheck files for _nonhuman.fastq files
    if (-f $fastqcheck || $filename=~/nonhuman/){
	$argref->{fastqcheck_location} = $fastqcheck;
    }
    else {
	die "No fastqcheck file $fastqcheck";
    }

    return $argref;
};


###############################################################################
# BUILDERS
###############################################################################

sub _get_reads {
    my ($self) = @_;
    return $self->_fastqcheck_stats->{reads};
}


sub _get_mean_quality {
    my ($self) = @_;
    return $self->_fastqcheck_stats->{mean_quality};
}


sub _get_basepairs {
    my ($self) = @_;
    return $self->_fastqcheck_stats->{basepairs};
}


sub _get_read_len {
    my ($self) = @_;
    return $self->_fastqcheck_stats->{read_len};
}


# Internal function to populate basepairs, reads, and mean_quality from
# fastqcheck
sub _get_fastqcheck_stats {
    my ($self) = @_;
    my $fqc = $self->fastqcheck_location;

    my %stats = (reads => 0,
                 basepairs => 0,
                 read_len => 0,
                 mean_quality => 0
                 );

    if (-f $fqc){
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

        my $read_tot = 0;
        my $bp_tot = 0;
        my $mean_q = 0;
        my $readlen = 0;
        my $q_tot = 0;
        my $q_n = 0;
        
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

        $stats{'reads'} = $read_tot;
        $stats{'basepairs'} = $bp_tot;
        $stats{'read_len'} = $readlen;
        $stats{'mean_quality'} = sprintf("%.1f",$mean_q);
    }
    return \%stats;
}


sub _get_md5 {
    my ($self) = @_;
    my $name = $self->name;
    open(my $MD5, "/software/solexa/bin/mpsa_download -m -f $name|") or die "Can't run mpsa_download to get md5: $!\n";
    my $line = <$MD5>;
    chomp $line;
    my ($md5, $file) = split / +/, $line;
    return $md5;
}

1;
