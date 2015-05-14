package Sfind::Lane; 
=head1 NAME

Sfind::Lane - Sequence Tracking Lane object

=head1 SYNOPSIS
    my $lane= Sfind::Lane->new({dbh => $dbh, id => $lane_id});

    my $id = $lane->id();

=head1 DESCRIPTION

An object describing the tracked properties of a lane.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=head2 new

  Arg [1]    : hashref: dbh => database handle to seqtracking database
                        id  => id_npg_information 
  Example    : my $lane = Sfind::Lane->new({dbh=>$dbh, id=>$id};)
  Description: Returns Lane object by lane_id
  Returntype : Sfind::Lane object


=head2 id

  Arg [1]    : None
  Example    : my $id = $lane->id();
  Description: Returns ID of a lane
  Returntype : SequenceScape ID (usu. integer)


=head2 created

  Arg [1]    : None
  Example    : my $created = $lane->created();
  Description: Returns created timestamp.  This is the 'run complete' timestamp from NPG
  Returntype : DateTime object


=head2 batch_id

  Arg [1]    : None
  Example    : my $batch_id = $lane->batch_id();
  Description: Returns batch ID of a lane
  Returntype : SequenceScape ID integer


=head2 run_name

  Arg [1]    : None
  Example    : my $run_name = $lane->run_name();
  Description: Returns run name of lane.  This is the name of the forward
		run for a paired-end lane
  Returntype : run name


=head2 name

  Arg [1]    : None
  Example    : my $lane_name = $lane->name();
  Description: 'name' of lane.  This is the concatenation of run_name with run_lane, e.g. '1234_1'
  Returntype : common name for lane


=head2 is_cancelled

  Arg [1]    : None
  Example    : my $cancelled = $lane->is_cancelled();
  Description: Returns whether lane is cancelled during run
  Returntype : boolean


=head2 is_paired

  Arg [1]    : None
  Example    : my $paired = $lane->is_paired();
  Description: Returns whether lane is paired-end read
  Returntype : boolean


=head2 pair_run_name

  Arg [1]    : None
  Example    : my $pair_run_name = $lane->pair_run_name();
  Description: Returns run name of reverse run for a paired-end run, if it has one.  Only old runs had separate run names for each direction.
  Returntype : run name


=head2 read_len

  Arg [1]    : None
  Example    : my $read_len = $lane->read_len();
  Description: Returns lane read_len
  Returntype : integer


=head2 run_lane

  Arg [1]    : None
  Example    : my $run_lane = $lane->run_lane();
  Description: Returns run_lane of Lane
  Returntype : integer


=head2 fastq

  Arg [1]    : None
  Example    : my $fastq = $lane->fastq();
  Description: Returns a ref to an array of the fastq objects that have been generated from this lane.
  Returntype : ref to array of Sfind::Fastq objects


=head2 fastq_filenames

  Arg [1]    : None
  Example    : my $fastq_files = $lane->fastq_filenames();
  Description: fetch array ref of mpsa fastq files for this lane
  Returntype : array ref of mpsa fuse file locations


=head2 get_fastq_by_filename

  Arg [1]    : None
  Example    : my $fastqfile = $lane->get_fastq_by_filename('/fuse/mpsafs/runs/1378/1378_s_1.fastq);
  Description: retrieve Sfind::Fastq object by file location 
  Returntype : Sfind::Fastq object


=head2 bam

  Arg [1]    : None
  Example    : my $bam = $lane->bam();
  Description: Returns a ref to an array of the bam objects that have been generated from this lane.
  Returntype : ref to array of Sfind::BAM objects


=head2 bam_filenames

  Arg [1]    : None
  Example    : my $bam_files = $lane->bam_filenames();
  Description: fetch array ref of irods bam files for this lane
  Returntype : array ref of irods file locations


=head2 get_bam_by_filename

  Arg [1]    : None
  Example    : my $bamfile = $lane->get_bam_by_filename('/seq/5322/5322_1.bam');
  Description: retrieve Sfind::BAM object by irods location
  Returntype : Sfind::BAM object


=head2 npg_qc

  Arg [1]    : None
  Example    : my $npg_qc = $lane->npg_qc();
  Description: get npg manual qc ('pass','fail','pending') for this lane.
  Returntype : npg_qc status string


=head2 mean_quality

  Arg [1]    : None
  Example    : my $mean_quality = $lane->mean_quality();
  Description: get mean_quality (from fastqcheck) in resulting fastq.  Note that this is the mean of the mean_qualities in the fastqcheck.
  Returntype : mean_quality to 1dp


=head2 reads

  Arg [1]    : None
  Example    : my $reads = $lane->reads();
  Description: get total number of reads in resulting fastq
  Returntype : integer count of reads


=head2 basepairs

  Arg [1]    : None
  Example    : my $basepairs = $lane->basepairs();
  Description: get total number of basepairs in resulting fastq
  Returntype : integer count of basepairs

=cut

use Moose;
use namespace::autoclean;
use Sfind::Types qw(MaybeMysqlDateTime);
use Sfind::Fastq;
use Sfind::BAM;
use VertRes::Wrapper::iRODS; # for finding bam files

has '_dbh'  => (
    is          => 'ro',
    isa         => 'DBI::db',
    required    => 1,
    init_arg    => 'dbh',
);

has 'id'    => (
    is          => 'ro',
    isa         => 'Int',
    required    => 1,
);

has 'run_name'=> (
    is          => 'ro',
    isa         => 'Str',
    required    => 1,
    init_arg    => 'id_run',
);

has 'pair_run_name'=> (
    is          => 'ro',
    isa         => 'Str',
    required    => 0,
    init_arg    => 'pair_id_run',
);

has 'run_lane'=> (
    is          => 'ro',
    isa         => 'Str',
    required    => 1,
    init_arg    => 'position',
);

has 'name'    => (
    is          => 'ro',
    isa         => 'Str',
    lazy        => 1,
    default     => sub {
            my $self = shift;
            return join "_", ($self->run_name, $self->run_lane);
                        },
);

has 'is_cancelled'  => (
    is          => 'ro',
    isa         => 'Bool',
    required    => 0,
    init_arg    => 'cancelled',
);

has 'is_paired'  => (
    is          => 'ro',
    isa         => 'Bool',
    required    => 0,
    init_arg    => 'paired_read',
);

has 'read_len'    => (
    is          => 'ro',
    isa         => 'Int',
    required    => 0,
);

has 'batch_id'    => (
    is          => 'ro',
    isa         => 'Int',
    required    => 0,
);

has 'npg_qc'    => (
    is          => 'ro',
    isa         => 'Str',
    required    => 1,
);

has 'created' => (
    is          => 'ro',
    isa         => MaybeMysqlDateTime,
    coerce      => 1,   # accept mysql dates
    init_arg    => 'run_complete',
);

has 'fastq'=> (
    is          => 'ro',
    isa         => 'ArrayRef[Sfind::Fastq]',
    lazy        => 1,
    builder     => '_get_fastq',
);

has 'fastq_filenames'=> (
    is          => 'ro',
    isa         => 'ArrayRef[Str]',
    lazy        => 1,
    builder     => '_get_fastq_filenames',
);

has 'bam'=> (
    is          => 'ro',
    isa         => 'ArrayRef[Sfind::BAM]',
    lazy        => 1,
    builder     => '_get_bam',
);

has 'bam_filenames'=> (
    is          => 'ro',
    isa         => 'ArrayRef[Str]',
    lazy        => 1,
    builder     => '_get_bam_filenames',
);


# Populate the parameters from the database
around BUILDARGS => sub {
    my $orig  = shift;
    my $class = shift;
    
    my $argref = $class->$orig(@_);

    die "Need to call with a npg_information id" unless $argref->{id};
    # id_run_pair = 0 means this is the first of any pair
    # which is what the srf, fastq & fastqcheck files are named for
    my $sql = qq[select * from npg_information 
                 where id_npg_information = ? 
                 and (id_run_pair=0 or id_run_pair is null)
                 ];
    my $row_ref = $argref->{dbh}->selectrow_hashref($sql, undef, ($argref->{id}));
    if ($row_ref){
        foreach my $field(keys %$row_ref){
            $argref->{$field} = $row_ref->{$field};
        }
    };
    # hacks

    # if paired, and only one runfolder, then cycles is total cycles (i.e.
    # fwd+rev) rather than one end.  Need to divide by two.
    if($row_ref->{paired_read} &! $row_ref->{has_two_runfolders}){
        $argref->{read_len} = int($row_ref->{cycles}/2);
    }
    else {
        $argref->{read_len} = $row_ref->{cycles};
    }

    # NPG manual qc
    
    #uninitialise error in .err files
    $argref->{npg_qc} = 'pending';
    if ($row_ref->{manual_qc}) {
	    if($row_ref->{manual_qc} == 1){
    	    $argref->{npg_qc} = 'pass';
    	}
    	elsif ($row_ref->{manual_qc} == 0){
        	$argref->{npg_qc} = 'fail';
    	}
    }



    # populate pair-end information if it exists
    $sql = qq[select id_run, run_complete 
              from npg_information 
              where batch_id = ? 
              and position = ? 
              and id_run_pair = ?];

    $row_ref = $argref->{dbh}->selectrow_hashref($sql, undef, ($argref->{batch_id},$argref->{position},$argref->{id_run}));
    if (keys %$row_ref){
        $argref->{pair_id_run} = $row_ref->{id_run};
        # the 'lane' was finished when the second end completed
        $argref->{run_complete} = $row_ref->{run_complete}; 
    }

    return $argref;
};


###############################################################################
# BUILDERS
###############################################################################


sub _get_fastq {
    my ($self) = @_;
    my @fastq;
    foreach my $name (@{$self->fastq_filenames()}){
        push @fastq, $self->get_fastq_by_filename($name);
    }
    return \@fastq;
}

sub _get_fastq_filenames {
    my ($self) = @_;
    my @fastq;
    return \@fastq;
}


sub get_fastq_by_filename {
    my ($self, $filename) = @_;
    my $obj = Sfind::Fastq->new($filename);
    return $obj;
}


sub _get_bam {
    my ($self) = @_;
    my @bam;
    foreach my $name (@{$self->bam_filenames()}){
        push @bam, $self->get_bam_by_filename($name);
    }
    return \@bam;
}


sub _get_bam_filenames {
    my ($self) = @_;
    my $run = $self->run_name;
    my $lane = $self->run_lane;
    my $irods = VertRes::Wrapper::iRODS->new();
    my @bam = @{$irods->find_files_by_run_lane($run,$lane)};
    return \@bam;
}


sub get_bam_by_filename {
    my ($self, $filename) = @_;
    my $obj = Sfind::BAM->new($filename);
    return $obj;
}


sub mean_quality {
    my ($self) = @_;
    unless ($self->{'mean_quality'}){
	$self->_get_fastq_stats;
    }
    return $self->{'mean_quality'};
}




sub reads {
    my ($self) = @_;
    unless ($self->{'reads'}){
	$self->_get_fastq_stats;
    }
    return $self->{'reads'};
}



sub basepairs {
    my ($self) = @_;
    unless ($self->{'basepairs'}){
	$self->_get_fastq_stats;
    }
    return $self->{'basepairs'};
}


# Internal function to populate basepairs, reads, and mean_quality from
# fastqcheck
# TODO change to get BAM stats instead if appropriate
# 2011-04-17 - Might pull this entirely, and the associated stat methods.
# No way of getting this for bam yet, and with MPSA going, this is obsolete.
sub _get_fastq_stats {
    my ($self) = @_;
}

1;
