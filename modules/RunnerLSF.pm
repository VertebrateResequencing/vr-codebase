package RunnerLSF;

use strict;
use warnings;
use Carp;
use DateTime;

our $Running = 1;
our $Error   = 2;
our $Unknown = 4;
our $No      = 8;
our $Done    = 16;

sub is_job_array_running
{
    my ($jids_file, $ids, $nmax) = @_;

    my @jobs = ();
    for my $id (@$ids) { push @jobs, {status=>$No, nfailures=>0}; }
    if ( ! -e $jids_file ) { return \@jobs; }

    # The same file with job ids may contain many job arrays runs. Check
    #   the jobs in the reverse order in case there were many failures.
    #   The last success counts, failures may be discarded in such a case.
    #
    open(my $fh, '<', $jids_file) or confess("$jids_file: $!");
    my $path;
    my @jids = ();
    while (my $line=<$fh>)
    {
        if ( !($line=~/^(\d+)\s+([^\t]*)/) ) { confess("Uh, could not parse \"$line\".\n") }
        push @jids, $1;     # LSF array ID
        if ( !defined $path ) { $path = $2; }
        if ( $path ne $2 ) { confess("$path ne $2\n"); }
    }
    for (my $i=@jids-1; $i>=0; $i--)
    {
        my $info = parse_bjobs_l($jids[$i]);
        if ( !defined $info ) { next; }
        for my $job (values %$info)
        {
            check_job($job);
        }
        for (my $i=0; $i<@$ids; $i++)
        {
            my $id = $$ids[$i];
            if ( !exists($$info{$id}) ) { next; }
            if ( $jobs[$i]{status} ne $No ) { next; }   # the job was submitted multiple times and already has a status
            if ( $$info{$id}{status}==$Done ) { $jobs[$i]{status} = $Done; }
            if ( $$info{$id}{status}==$Running ) { $jobs[$i]{status} = $Running; }
        }
    }
    my $ntodo = 0;
    for (my $i=0; $i<@$ids; $i++)
    {
        if ( $jobs[$i]{status} & $Running || $jobs[$i]{status} & $Error ) { $ntodo++; }
        if ( $nmax && $ntodo >= $nmax ) { last; } 
        if ( $jobs[$i]{status} ne $No ) { next; }
        my $info = parse_output($$ids[$i], $path);
        if ( defined $info )
        {
            $jobs[$i]{status} = $$info{status};
            $jobs[$i]{nfailures} = $$info{nfailures};
        }
    }
    return \@jobs;
}

sub parse_bjobs_l
{
    my ($jid) = @_;

    my @lines;
    for (my $i=0; $i<3; $i++)
    {
        @lines = `bjobs -l $jid 2>/dev/null`;
        if ( $? ) { sleep 5; next; }
        if ( !scalar @lines ) { return undef; }
    }

    my %months = qw(Jan 1 Feb 2 Mar 3 Apr 4 May 5 Jun 6 Jul 7 Aug 8 Sep 9 Oct 10 Nov 11 Dec 12);
    my $year = (gmtime())[5] + 1900;

    my $info   = {};
    my $job;
    for (my $i=0; $i<@lines; $i++)
    {
        if ( $lines[$i]=~/^\s*$/ ) { next; }

        my $job_info;
        if ( $lines[$i]=~/^Job <(\d+)/ )
        {
            # Runner's ID is $id, LSF job ID is lsf_id
            my $id = $1;
            my $rest = $';
            my $lsf_id = $id;
            if ( $rest=~/^\[(\d+)/ )
            {
                $lsf_id = "$id\[$1\]";
                $id = $1;
            }

            if ( scalar keys %$job) { $$info{$$job{id}} = $job; }
            $job = { id=>$id, lsf_id=>$lsf_id };

            my $job_info = $lines[$i];
            chomp($job_info);
            $i++; 

            while ( $i<@lines && $lines[$i]=~/^\s{21}?/ )
            {
                $job_info .= $';
                chomp($job_info);
                $i++;
            }
            if ( !($job_info=~/,\s*Status <([^>]+)>/) ) { confess("Could not determine the status: [$job_info]"); }
            $$job{status} = $1;
            if ( !($job_info=~/,\s*Queue <([^>]+)>/) ) { confess("Could not determine the queue: [$job_info]"); }
            $$job{queue} = $1;
        }
        # Tue Mar 19 13:00:35: [685] started on <uk10k-4-1-07>...
        if ( $lines[$i]=~/^\w+\s+(\w+)\s+(\d+) (\d+):(\d+):(\d+):.+ started on </ ) 
        {
            $$job{started} = DateTime->new(month=>$months{$1}, day=>$2, hour=>$3, minute=>$4, year=>$year)->epoch;
        }
        # Tue Mar 19 13:58:23: Resource usage collected...
        if ( $lines[$i]=~/^\w+\s+(\w+)\s+(\d+) (\d+):(\d+):(\d+):\s+Resource usage collected/ ) 
        {
            if ( !exists($$job{started}) ) { confess("No wall time for job $$job{id}??", @lines); }
            my $wall_time = DateTime->new(month=>$months{$1}, day=>$2, hour=>$3, minute=>$4, year=>$year)->epoch - $$job{started};
            if ( !exists($$job{cpu_time}) or $$job{cpu_time} < $wall_time ) { $$job{cpu_time} = $wall_time; }
        }
        if ( $lines[$i]=~/The CPU time used is (\d+) seconds./ ) 
        { 
            if ( !exists($$job{cpu_time}) or $$job{cpu_time} < $1 ) { $$job{cpu_time} = $1; }
        }
    }
    if ( scalar keys %$job) { $$info{$$job{id}} = $job; }
    return $info;
}

sub check_job
{
    my ($job) = @_;
    my $status = { DONE=>$Done, PEND=>$Running, EXIT=>$Error, RUN=>$Running, UNKWN=>$Running, SSUSP=>$Running };
    if ( !exists($$status{$$job{status}}) ) 
    { 
        if ( $$job{status} eq 'ZOMBI' )
        {
            confess(
                    "FIXME: \n" .
                    "Some of the jobs are in ZOMBI state, please see `man bjobs` for possible reasons why this happened.\n" .
                    "Because the job might have been already requeued with a new job ID which this pipeline would be unaware of,\n" .
                    "the pipeline cannot proceed. Please fix manually by deleting the appropriate line in .jobs/*.w.jid\n" .
                    "and remove the appropriate LSF output file(s) .jobs/*.w.<ID>.o\n"
                   );
        }
        confess("Todo: $$job{status}\n"); 
    }
    $$job{status} = $$status{$$job{status}};
    if ( $$job{status}==$Running )
    {
        if ( !exists($$job{cpu_time}) ) { $$job{cpu_time} = 0; }
        my $queue = $$job{queue};
        my %queue_limits = ( basement=>1e9, long=>2880*60, normal=>720*60, small=>30*60 );
        if ( exists($queue_limits{$queue}) && $$job{cpu_time}*1.3 > $queue_limits{$queue} )
        {
            my $bswitch;
            for my $q (sort {$queue_limits{$a} <=> $queue_limits{$b}} keys %queue_limits)
            {
                if ( $$job{cpu_time}*1.3 > $queue_limits{$q} ) { next; }
                warn("Changing queue of the job $$job{lsf_id} from $$job{queue} to $q\n");
                `bswitch $q '$$job{lsf_id}'`;
                $$job{queue} = $q;
                last;
            }
        }
    }
}


sub parse_output
{
    my ($jid,$output) = @_;

    my $fname = "$output.$jid.o";
    if ( !-e $fname ) { return undef; }
    
    my $out = {};
    open(my $fh,'<',$fname) or confess("$fname: $!");
    while (my $line=<$fh>)
    {
        # Subject: Job 822187: <_2215_1_graphs> Done
        if ( $line =~ /^Subject: Job.+\s+(\S+)$/ )
        {
            if ( $1 eq 'Exited' ) { $$out{status} = $Error; $$out{nfailures}++; }
            if ( $1 eq 'Done' ) { $$out{status} = $Done; $$out{nfailures} = 0; }
        }
    }
    close($fh);
    return $out;
}

sub past_limits
{
    my ($jid,$output) = @_; 
    my $fname = "$output.$jid.o";
    if ( ! -e $fname ) { return (); }
    open(my $fh,'<',$fname) or confess("$fname: $!");
    my (%out,$killed,$mem);
    while (my $line=<$fh>)
    {
        if ( $line=~/^TERM_MEMLIMIT:/) { $killed = 1; }
        elsif ( $line=~/^\s+Max Memory\s+:\s+(\S+)\s+(\S+)/) 
        { 
            $mem = $1;
            if ($2 eq 'KB') { $mem /= 1024; }
            elsif ($2 eq 'GB') { $mem *= 1024; }

            if ( !exists($out{memory}) or $out{memory}<$mem )
            {
                $out{memory} = $mem;
                if ( $killed ) { $out{MEMLIMIT} = $mem; }
                else { delete($out{MEMLIMIT}); }
            }
        }
    }
    close($fh);
    return %out;
}


sub run_array
{
    my ($jids_file, $job_name, $opts, $cmd, $ids) = @_;

    if ( !scalar @$ids ) { confess("No IDs given??\n"); }

    # Process the list of IDs. The maximum job name length is 255 characters.
    my @ids = sort { $a<=>$b } @$ids;
    my @bsub_ids;
    my $from = $ids[0];
    my $prev = $from;
    for (my $i=1; $i<@ids; $i++)
    {
        my $id = $ids[$i];
        if ( $id != $prev+1 )
        {
            if ( $prev>$from ) { push @bsub_ids, "$from-$prev"; }
            else { push @bsub_ids, $from; }
            $from = $id;
            $prev = $id;
        }
        $prev = $id;
    }
    if ( $prev>$from ) { push @bsub_ids, "$from-$prev"; }
    else { push @bsub_ids, $from; }
    my $bsub_ids  = join(',', @bsub_ids);
    while ( length($job_name) + length($bsub_ids) > 250 && scalar @bsub_ids ) 
    {
        pop(@bsub_ids);
        $bsub_ids = join(',', @bsub_ids);
    }

  
    $cmd =~ s/{JOB_INDEX}/\$LSB_JOBINDEX/g;
    my $bsub_opts = '';
    if ( exists($$opts{memory}) && $$opts{memory} ) 
    { 
        my $mem = int($$opts{memory});
        if ( $mem )
        {
            $bsub_opts = sprintf " -M%d -R 'select[type==X86_64 && mem>%d] rusage[mem=%d]'", $mem*1000,$mem,$mem; 
        }
    }
    my $bsub_cmd  = qq[bsub -J '${job_name}[$bsub_ids]' -e $job_name.\%I.e -o $job_name.\%I.o $bsub_opts '$cmd'];

    # Submit to LSF
    print STDERR "$bsub_cmd\n";
    my @out = `$bsub_cmd`;
    if ( scalar @out!=1 || !($out[0]=~/^Job <(\d+)> is submitted/) )
    {
        my $cwd = `pwd`;
        confess("Expected different output from bsub. The command was:\n\t$cmd\nThe working directory was:\n\t$cwd\nThe output was:\n", @out);
    }

    # Write down info about the submitted command
    my $jid = $1;
    open(my $jids_fh, '>>', $jids_file) or confess("$jids_file: $!");
    print $jids_fh "$jid\t$job_name\t$bsub_cmd\n";
    close $jids_fh;
}


1;

