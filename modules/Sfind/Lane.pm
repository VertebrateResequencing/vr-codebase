package Sfind::Lane; 
=head1 NAME

Sfind::Lane - Sequence Tracking Lane object

=head1 SYNOPSIS
    my $lane= Sfind::Lane->new($dbh, $lane_id);

    my $id = $lane->id();

=head1 DESCRIPTION

An object describing the tracked properties of a lane.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;
no warnings 'uninitialized';
BEGIN { push(@INC,qw(/software/solexa/lib/site_perl/5.8.8)) }
use npg::api::run_lane;
use Sfind::Fastq;

=head2 new

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : lane id
  Example    : my $lane= Sfind::Sfind->new($dbh, $id)
  Description: Returns Lane object by lane_id
  Returntype : Sfind::Lane object

=cut

sub new {
    my ($class,$dbh, $id) = @_;
    die "Need to call with a db handle and id" unless ($dbh && $id);
    my $self = {};
    bless ($self, $class);
    $self->{_dbh} = $dbh;

    # id_run_pair = 0 means this is the first of any pair
    # which is what the srf, fastq & fastqcheck files are named for
    my $sql = qq[select batch_id, id_run, position, run_complete, cycles, paired_read, has_two_runfolders from npg_information where id_npg_information = ? and (id_run_pair=0 or id_run_pair is null)];
    my $sth = $self->{_dbh}->prepare($sql);

    $sth->execute($id);
    my $data = $sth->fetchall_arrayref()->[0];	# only one row for each id
    unless ($data){
	return undef;
    }
    $self->id($id);
    $self->batch_id($data->[0]);
    $self->run_name($data->[1]);
    $self->run_lane($data->[2]);
    $self->created($data->[3]);
    $self->is_paired($data->[5]);
    if ($data->[5] &! $data->[6]){
        # if paired, and only one runfolder, then cycles is total cycles (i.e.
        # fwd+rev) rather than one end.  Need to divide by two.
        $self->read_len(int($data->[4]/2));
    }
    else {
        $self->read_len($data->[4]);
    }
    $self->_load_pair_info();
    return $self;
}


=head2 new_by_run_lane

  Arg [1]    : database handle to seqtracking database
  Arg [2]    : lane id
  Example    : my $lane= Sfind::Sfind->new_by_run_lane($dbh, '2064', '3')
  Description: Returns Lane object from a run and lane
  Returntype : Sfind::Lane object

=cut

sub new_by_run_lane {
    my ($class,$dbh, $run, $lane) = @_;
    die "Need to call with a db handle, run and lane" unless ($dbh && $run && $lane);
    my $self = {};
    bless ($self, $class);
    $self->{_dbh} = $dbh;

    # id_run_pair = 0 means this is the first of any pair
    # which is what the srf, fastq & fastqcheck files are named for
    my $sql = qq[select batch_id, id_run, position, run_complete, cycles, paired_read, id_npg_information from npg_information where id_run = ? and position = ? and id_run_pair=0];
    my $sth = $self->{_dbh}->prepare($sql);

    $sth->execute($run, $lane);
    my $data = $sth->fetchall_arrayref()->[0];	# only one row for each id
    unless ($data){
	return undef;
    }
    $self->batch_id($data->[0]);
    $self->run_name($data->[1]);
    $self->run_lane($data->[2]);
    $self->created($data->[3]);
    $self->read_len($data->[4]);
    $self->is_paired($data->[5]);
    $self->id($data->[6]);

    $self->_load_pair_info();
    return $self;
}


=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $lane->id();
	       $lane->id('104');
  Description: Get/Set for ID of a lane
  Returntype : SequenceScape ID (usu. integer)

=cut

sub id {
    my ($self,$id) = @_;
    if (defined $id and $id ne $self->{'id'}){
	$self->{'id'} = $id;
    }
    return $self->{'id'};
}


=head2 created

  Arg [1]    : created (optional)
  Example    : my $created = $lane->created();
	       $lane->created('2008-10-24 11:40:59');
  Description: Get/Set for created timestamp.  This is the 'run complete' timestamp from NPG
  Returntype : timestamp string

=cut

sub created {
    my ($self,$created) = @_;
    if (defined $created and $created ne $self->{'created'}){
	$self->{'created'} = $created;
    }
    return $self->{'created'};
}


=head2 batch_id

  Arg [1]    : batch_id (optional)
  Example    : my $batch_id = $lane->batch_id();
	       $lane->batch_id('104');
  Description: Get/Set for batch ID of a lane
  Returntype : SequenceScape ID integer

=cut

sub batch_id {
    my ($self,$batch_id) = @_;
    if (defined $batch_id and $batch_id ne $self->{'batch_id'}){
	$self->{'batch_id'} = $batch_id;
    }
    return $self->{'batch_id'};
}



=head2 run_name

  Arg [1]    : run_name (optional)
  Example    : my $run_name = $lane->run_name();
	       $lane->run_name('1234');
  Description: Get/Set for run name of lane.  This is the name of the forward
		run for a paired-end lane
  Returntype : run name

=cut

sub run_name {
    my ($self,$run_name) = @_;
    if (defined $run_name and $run_name ne $self->{'run_name'}){
	$self->{'run_name'} = $run_name;
    }
    return $self->{'run_name'};
}


=head2 name

  Arg [1]    : None
  Example    : my $lane_name = $lane->name();
  Description: 'name' of lane.  This is the concatenation of run_name with run_lane, e.g. '1234_1'
  Returntype : common name for lane

=cut

sub name {
    my ($self) = @_;
    return join "_", ($self->run_name, $self->run_lane);
}


=head2 is_paired

  Arg [1]    : boolean for whether lane is paired (optional)
  Example    : my $paired = $lane->is_paired();
	       $lane->is_paired('1');
  Description: Get/Set for whether lane is paired-end read
  Returntype : boolean

=cut

sub is_paired {
    my ($self,$is_paired) = @_;
    if (defined $is_paired and $is_paired ne $self->{'is_paired'}){
	$self->{'is_paired'} = $is_paired ? 1 : 0;
    }
    return $self->{'is_paired'};
}


=head2 pair_run_name

  Arg [1]    : pair_run_name (optional)
  Example    : my $pair_run_name = $lane->pair_run_name();
	       $lane->pair_run_name('1234');
  Description: Get/Set for run name of reverse run for a paired-end run
  Returntype : run name

=cut

sub pair_run_name {
    my ($self,$pair_run_name) = @_;
    if (defined $pair_run_name and $pair_run_name ne $self->{'pair_run_name'}){
	$self->{'pair_run_name'} = $pair_run_name;
    }
    return $self->{'pair_run_name'};
}



=head2 read_len

  Arg [1]    : read_len (optional)
  Example    : my $read_len = $lane->read_len();
	       $lane->read_len(54);
  Description: Get/Set for lane read_len
  Returntype : integer

=cut

sub read_len {
    my ($self,$read_len) = @_;
    if (defined $read_len and $read_len ne $self->{'read_len'}){
	$self->{'read_len'} = $read_len;
    }
    return $self->{'read_len'};
}


=head2 run_lane

  Arg [1]    : run_lane (optional)
  Example    : my $run_lane = $lane->run_lane();
	       $run_lane->lane(7);
  Description: Get/Set for run_lane of Lane
  Returntype : integer

=cut

sub run_lane {
    my ($self,$run_lane) = @_;
    if (defined $run_lane and $run_lane ne $self->{'run_lane'}){
	$self->{'run_lane'} = $run_lane;
    }
    return $self->{'run_lane'};
}

=head2 fastq

  Arg [1]    : None
  Example    : my $fastq = $lane->fastq();
  Description: Returns a ref to an array of the fastq objects that have been generated from this lane.
  Returntype : ref to array of Sfind::Fastq objects

=cut

sub fastq {
    my ($self) = @_;
    unless ($self->{'fastq'}){
	my @fastq;
    	foreach my $name (@{$self->fastq_filenames()}){
	    push @fastq, $self->get_fastq_by_filename($name);
	}
	$self->{'fastq'} = \@fastq;
    }
    return $self->{'fastq'};
}


=head2 fastq_filenames

  Arg [1]    : none
  Example    : my $fastq_files = $lane->fastq_filenames();
  Description: fetch array ref of mpsa fastq files for this lane
  Returntype : array ref of mpsa fuse file locations

=cut

sub fastq_filenames {
    my ($self) = @_;
    unless ($self->{'fastq_filenames'}){
	my $run = $self->run_name;
	my $lane = $self->run_lane;
	open(my $FASTQ, "/software/solexa/bin/dfind -run $run -lane $lane -filetype fastq |grep fuse |") or die "Can't run dfind to locate fastq: $!\n";
	my @fastq;
	while (my $fq = <$FASTQ>){
	    chomp $fq;
	    push @fastq, $fq;
	}
	close $FASTQ;
	$self->{'fastq_filenames'} = \@fastq;
    }
    return $self->{'fastq_filenames'};
}


=head2 get_fastq_by_filename

  Arg [1]    : fastq filename, including path
  Example    : my $fastqfile = $lane->get_fastq_by_filename('/fuse/mpsafs/runs/1378/1378_s_1.fastq);
  Description: retrieve Sfind::Fastq object by id
  Returntype : Sfind::Fastq object

=cut

sub get_fastq_by_filename {
    my ($self, $filename) = @_;
    my $obj = Sfind::Fastq->new($filename);
    return $obj;
}


=head2 npg_qc

  Arg [1]    : none
  Example    : my $npg_qc = $lane->npg_qc();
  Description: get npg manual qc for this lane.
  Returntype : npg_qc status string

=cut

sub npg_qc {
    my ($self) = @_;
    unless ($self->{'npg_qc'}){
        my $rl=npg::api::run_lane->new({id_run=>$self->run_name,position=>$self->run_lane}); 
        if ($rl){
            $self->{'npg_qc'} = $rl->manual_qc;
        }
        else {
            warn "Can't get NPG run_lane for ",$self->run_name,"_",$self->run_lane,"\n";
        }
    }
    return $self->{'npg_qc'};
}


=head2 mean_quality

  Arg [1]    : none
  Example    : my $mean_quality = $lane->mean_quality();
  Description: get mean_quality (from fastqcheck) in resulting fastq.  Note that this is the mean of the mean_qualities in the fastqcheck.
  Returntype : mean_quality to 1dp

=cut

sub mean_quality {
    my ($self) = @_;
    unless ($self->{'mean_quality'}){
	$self->_get_fastq_stats;
    }
    return $self->{'mean_quality'};
}


=head2 reads

  Arg [1]    : none
  Example    : my $reads = $lane->reads();
  Description: get total number of reads in resulting fastq
  Returntype : integer count of reads

=cut

sub reads {
    my ($self) = @_;
    unless ($self->{'reads'}){
	$self->_get_fastq_stats;
    }
    return $self->{'reads'};
}


=head2 basepairs

  Arg [1]    : none
  Example    : my $basepairs = $lane->basepairs();
  Description: get total number of basepairs in resulting fastq
  Returntype : integer count of basepairs

=cut

sub basepairs {
    my ($self) = @_;
    unless ($self->{'basepairs'}){
	$self->_get_fastq_stats;
    }
    return $self->{'basepairs'};
}


# Internal function to populate basepairs, reads, and mean_quality from
# fastqcheck
sub _get_fastq_stats {
    my ($self) = @_;
    my $fastq = $self->fastq;
    my $read_tot = 0;
    my $bp_tot = 0;
    my $mean_q = 0;
    my $q_tot;
    my $count;
    foreach my $fq(@$fastq){
        $count++;
        $q_tot += $fq->mean_quality;
        $read_tot += $fq->reads;
        $bp_tot += $fq->basepairs;
    }

    if ($count){
	$mean_q = $q_tot/$count;
    }
    $self->{'reads'} = $read_tot;
    $self->{'basepairs'} = $bp_tot;
    $self->{'mean_quality'} = sprintf("%.1f",$mean_q);
}


# Internal function to populate pair-end information if it exists

sub _load_pair_info {
    my ($self) = @_;

    my $sql = qq[select id_run, run_complete, cycles from npg_information where batch_id = ? and position = ? and id_run_pair = ?];
    my $sth = $self->{_dbh}->prepare($sql);

    $sth->execute($self->batch_id, $self->run_lane, $self->run_name);
    my $data = $sth->fetchall_arrayref()->[0];	# only one row for each id
    if ($data){
	$self->pair_run_name($data->[0]);
	$self->created($data->[1]);    # the 'lane' was finished when the
					    # second end completed
    }
}
1;
