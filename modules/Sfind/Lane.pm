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

#use Sfind::Lane;

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
    # which is what the srf,fastq & fastqcheck files are named for
    my $sql = qq[select batch_id, id_run, position, run_complete, cycles from npg_information where id_npg_information = ? and id_run_pair=0];
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
    $self->read_len($data->[4]);

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

  Arg [1]    : none
  Example    : my $fastq = $lane->fastq();
  Description: fetch array ref of mpsa fastq files for this lane
  Returntype : array ref of mpsa fuse file locations

=cut

sub fastq {
    my ($self) = @_;
    unless ($self->{'fastq'}){
	my $run = $self->run_name;
	my $lane = $self->run_lane;
	open(my $FASTQ, "/software/solexa/bin/dfind -run $run -lane $lane -filetype fastq |grep fuse |") or die "Can't run dfind to locate fastq: $!\n";
	my @fastq;
	while (my $fq = <$FASTQ>){
	    chomp $fq;
	    push @fastq, $fq;
	}
	close $FASTQ;
	$self->{'fastq'} = \@fastq;
    }
    return $self->{'fastq'};
}


=head2 mean_quality

  Arg [1]    : none
  Example    : my $mean_quality = $lane->mean_quality();
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


=head2 reads

  Arg [1]    : none
  Example    : my $reads = $lane->reads();
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


=head2 basepairs

  Arg [1]    : none
  Example    : my $basepairs = $lane->basepairs();
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


# Internal function to populate basepairs, reads, and mean_quality from
# fastqcheck
sub _get_fastqcheck_stats {
    my ($self) = @_;
    my $fastq = $self->fastq;
    my $read_tot = 0;
    my $bp_tot = 0;
    my $mean_q = 0;
    my $q_tot = 0;
    my $q_n = 0;
    foreach my $fq(@$fastq){
	my $fqc = $fq."check";	# should really dfind, but this is faster
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
	# TODO check cycles from here too?
	my ($reads, $foo, $bp) = split " ", $header;
	$read_tot += $reads;
	$bp_tot += $bp;

	my ($label, $atot, $ctot, $gtot, $ttot, $ntot, @vals) = split " ", $totals;	
	my $AQ = pop @vals; # NB, this is NOT average quality
	# the numbers in @vals are the counts (in millions I think) of bases
	# at the quality of the subscript (i.e. from 0..x)
	# want average of these
	for (my $i=0; $i < scalar @vals; ++$i){
	    $q_n += $vals[$i]; 
	    $q_tot += $vals[$i] * $i;
	}
    }
    if ($q_n){
	$mean_q = $q_tot/$q_n;
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
