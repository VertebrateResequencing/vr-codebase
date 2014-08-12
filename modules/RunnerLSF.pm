package RunnerLSF;

use strict;
use warnings;
use Carp;
use DateTime;

sub new
{
    my ($class,@args) = @_;
    my $self = @args ? {@args} : {};
    bless $self, ref($class) || $class;

    $$self{Running} = 1;
    $$self{Error}   = 2;
    $$self{Zombi}   = 4;
    $$self{No}      = 8;
    $$self{Done}    = 16;
    $$self{Waiting} = 32;

    $$self{lsf_status_codes} = 
    { 
        DONE  => $$self{Done}, 
        PEND  => $$self{Running} | $$self{Waiting}, 
        WAIT  => $$self{Running} | $$self{Waiting}, 
        EXIT  => $$self{Error}, 
        ZOMBI => $$self{Zombi},
        RUN   => $$self{Running}, 
        UNKWN => $$self{Running}, 
        SSUSP => $$self{Running}
    };

    # runtime and queue_limits are in minutes
    $$self{default_limits} = { runtime=>40, memory=>1_000, queue=>'normal' };
    $$self{queue_limits}   = { basement=>1e9, long=>48*60, normal=>12*60, small=>30 };

    $self->_set_lsf_limits_unit();
    $self->_init_zombies();

    return $self;
}

sub job_running
{
    my ($self,$task) = @_;
    return $$task{status} & $$self{Running};
}
sub job_done
{
    my ($self, $task) = @_;
    return $$task{status} & $$self{Done};
}
sub job_failed
{
    my ($self, $task) = @_;
    return $$task{status} & $$self{Error};
}
sub job_nfailures
{
    my ($self, $task) = @_;
    return $$task{nfailures} ? $$task{nfailures} : 0;
}
sub set_max_jobs
{
    my ($self,$nmax) = @_;
    $$self{nmax_jobs} = $nmax;
}
# runtime and queue_limits are in minutes
sub set_limits
{
    my ($self,%limits) = @_;
    $$self{limits} = { %{$$self{default_limits}}, %limits };
    if ( exists($limits{queues}) )
    {
        $$self{queue_limits} = { %{$limits{queues}} };
    }
}
sub clean_jobs
{
    my ($self,$wfile,$ids) = @_;
}
sub kill_job
{
    my ($self,$job) = @_;
    if ( !exists($$job{lsf_id}) ) { return; }
    my $cmd = "bkill -s7 -b '$$job{lsf_id}'";
    warn("$cmd\n");
    `$cmd`;
}

sub _init_zombies
{
    my ($self) = @_;
    if ( !exists($ENV{LSF_ZOMBI_IS_DEAD}) ) { return; }
    my @arrays = split(/\./,$ENV{LSF_ZOMBI_IS_DEAD});
    for my $arr (@arrays)
    {
        if ( !($arr=~/^(\d+)\[(.+)\]$/) ) { confess("Could not parse LSF_ZOMBI_IS_DEAD=$ENV{LSF_ZOMBI_IS_DEAD}\n"); }
        my $id  = $1;
        my $ids = $2;
        my @items = split(/,/,$ids);
        for my $item (@items)
        {
            if ( $item=~/^(\d+)$/ ) { $$self{ignore_zombies}{"${id}[$1]"} = 1; next; }
            my ($from,$to) = split(/-/,$item); 
            for (my $i=$from; $i<=$to; $i++)
            {
                $$self{ignore_zombies}{"${id}[$i]"} = 1;
            }
        }
    }
}

sub get_jobs
{
    my ($self, $jids_file, $ids) = @_;

    # For each input job id create a hash with info: status, number of failuers
    my @jobs_out = ();
    for my $id (@$ids) { push @jobs_out, {status=>$$self{No}, nfailures=>0}; }
    if ( ! -e $jids_file ) { return \@jobs_out; }

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

    # ZOMBI jobs need special care, we cannot be 100% sure that the non-responsive
    # node stopped writing to network disks. Let the user decide if they can be 
    # safely ignored.
    my %zombi_warning = ();

    # Get info from bjobs -l: iterate over LSF array IDs we remember for this task in the jids_file
    for (my $i=@jids-1; $i>=0; $i--)
    {
        my $info = $self->_parse_bjobs_l($jids[$i]);
        if ( !defined $info ) { next; }

        # Check if the time limits are still OK for all running jobs. Switch queues if not
        for my $job_l (values %$info)
        {
            if ( $$job_l{status} eq 'ZOMBI' && !$$self{ignore_zombies}{$$job_l{lsf_id}} ) { $zombi_warning{$$job_l{lsf_id}} = 1; }
            $self->_check_job($job_l,$jids_file);
        }

        # Update status of input jobs present in the bjobs -l listing. Note that the failed jobs
        # get their status from the output files as otherwise we wouldn't know how many times
        # they failed already. 
        for (my $j=0; $j<@$ids; $j++)
        {
            my $id = $$ids[$j];
            if ( !exists($$info{$id}) ) { next; }
            if ( $jobs_out[$j]{status} ne $$self{No} ) { next; }   # the job was submitted multiple times and already has a status

            if ( $$info{$id}{status} & $$self{Done} ) 
            { 
                $jobs_out[$j]{status} = $$self{Done}; 
            }
            elsif ( $$info{$id}{status} & $$self{Running} ) 
            { 
                $jobs_out[$j]{status} = $$self{Running}; 
                $jobs_out[$j]{lsf_id} = $$info{$id}{lsf_id};
            }
            elsif ( $$info{$id}{status} & $$self{Zombi} )
            {
                # Set as failed if ZOMBI should be ignored, otherwise say it's running.
                my $lsf_id = $$info{$id}{lsf_id};
                if ( $$self{ignore_zombies}{$lsf_id} ) { $jobs_out[$j]{status} = $$self{Error}; }
                else { $jobs_out[$j]{status} = $$self{Running}; }
                $jobs_out[$j]{lsf_id} = $lsf_id;
            }
        }
    }

    if ( scalar keys %zombi_warning )
    {
        my %arrays = ();
        for my $lsf_id (keys %zombi_warning)
        {
            if ( $lsf_id =~ s/^(\d+)\[(\d+)\]$// ) { push @{$arrays{$1}},$2; }
        }
        my @id_strings = ();
        for my $id (keys %arrays)
        {
            while ( @{$arrays{$id}} )
            {
                push @id_strings, "${id}[". $self->_create_bsub_ids_string($id,$arrays{$id}) . "]";
            }
        }
        warn(
            "\n\n----\n\n" .
            "WARNING:  Some jobs were found in ZOMBI state and are still considered as\n" .
            "  running by the pipeline. To ignore these jobs, set the environment variable\n" .
            "       LSF_ZOMBI_IS_DEAD='" .join('.',@id_strings). "'\n" .
            "  and restart the pipeline. See also \"$jids_file\".\n" .
            "----\n\n"
            );
    }

    # For jobs which are not present in the bjobs -l listing we get the info from output files
    my $ntodo = 0;
    for (my $i=0; $i<@$ids; $i++)
    {
        if ( $jobs_out[$i]{status} & $$self{Running} || $jobs_out[$i]{status} & $$self{Error} ) { $ntodo++; }
        if ( $$self{nmax_jobs} && $ntodo >= $$self{nmax_jobs} ) { last; } 
        if ( $jobs_out[$i]{status} ne $$self{No} ) { next; }
        my $info = $self->_parse_output($$ids[$i], $path);
        if ( defined $info )
        {
            $jobs_out[$i]{status} = $$info{status};
            $jobs_out[$i]{nfailures} = $$info{nfailures};
        }
    }
    return \@jobs_out;
}

sub _parse_bjobs_l
{
    my ($self,$jid) = @_;

    my @lines;
    for (my $i=0; $i<3; $i++)
    {
        @lines = `bjobs -l $jid 2>/dev/null`;
        if ( $? ) { sleep 5; next; }
        if ( !scalar @lines ) { return undef; }
    }

    my %months = qw(Jan 1 Feb 2 Mar 3 Apr 4 May 5 Jun 6 Jul 7 Aug 8 Sep 9 Oct 10 Nov 11 Dec 12);
    my $year = (gmtime())[5] + 1900;

    my $info = {};
    my $job;
    for (my $i=0; $i<@lines; $i++)
    {
        if ( $lines[$i]=~/^\s*$/ ) { next; }

        my $job_info;
        if ( $lines[$i]=~/^Job <(\d+)(.*)$/ )
        {
            # Runner's ID is $id, LSF job ID is lsf_id
            my $id = $1;
            my $rest = $2;
            my $lsf_id = $id;
            if ( $rest=~/^\[(\d+)/ )
            {
                $lsf_id = "$id\[$1\]";
                $id = $1;
            }

            if ( scalar keys %$job) { $$info{$$job{id}} = $job; }
            $job = { id=>$id, lsf_id=>$lsf_id, cpus=>1 };

            my $job_info = $lines[$i];
            chomp($job_info);
            $i++; 

            while ( $i<@lines && $lines[$i]=~/^\s{21}?(.*)$/ )
            {
                $job_info .= $1;
                chomp($job_info);
                $i++;
            }
            if ( !($job_info=~/,\s*Status <([^>]+)>/) ) { confess("Could not determine the status: [$job_info]"); }
            $$job{status} = $1;
            if ( !($job_info=~/,\s*Queue <([^>]+)>/) ) { confess("Could not determine the queue: [$job_info]"); }
            $$job{queue} = $1;
            if ( !($job_info=~/,\s*Command <([^>]+)>/) ) { confess("Could not determine the command: [$job_info]"); }
            $$job{command} = $1;
        }
        
        # Collect also checkpoint data for LSFCR to avoid code duplication: checkpoint directory, memory, status
        # Wed Mar 19 10:14:17: Submitted from host <vr-2-2-02>...
        if ( $lines[$i]=~/^\w+\s+\w+\s+\d+ \d+:\d+:\d+:\s*Submitted from/ )
        {
            my $job_info = $lines[$i];
            chomp($job_info);
            $i++;

            while ( $i<@lines && $lines[$i]=~/^\s{21}?(.*)$/ )
            {
                $job_info .= $1;
                chomp($job_info);
                $i++;
            }
            if ( $job_info=~/,\s*Checkpoint directory <([^>]+)>/ ) { $$job{chkpnt_dir} = $1; }
            if ( $job_info=~/\srusage\[mem=(\d+)/ ) 
            { 
                $$job{mem_usage} = $1; 
                if ( $$self{lsf_limits_unit} eq 'kB' ) { $$job{mem_usage} /= 1000.0; }
            }
        }
        # elsif ( $lines[$i]=~/^\w+\s+\w+\s+\d+ \d+:\d+:\d+:\s*Completed <exit>; TERM_CHKPNT/ )
        # {
        #     $$job{status} = 'EXIT';
        # }

        # Tue Mar 19 13:00:35: [685] started on <uk10k-4-1-07>...
        # Tue Dec 24 13:12:00: [1] started on 8 Hosts/Processors <8*vr-1-1-05>...
        elsif ( $lines[$i]=~/^\w+\s+(\w+)\s+(\d+) (\d+):(\d+):(\d+):.+ started on/ ) 
        {
            $$job{started} = DateTime->new(month=>$months{$1}, day=>$2, hour=>$3, minute=>$4, year=>$year)->epoch;
        }
        elsif ( $lines[$i]=~/^\w+\s+(\w+)\s+(\d+) (\d+):(\d+):(\d+):.+ dispatched to/ )    # associated with underrun status
        {
            $$job{started} = DateTime->new(month=>$months{$1}, day=>$2, hour=>$3, minute=>$4, year=>$year)->epoch;
        }
        # Tue Mar 19 13:58:23: Resource usage collected...
        elsif ( $lines[$i]=~/^\w+\s+(\w+)\s+(\d+) (\d+):(\d+):(\d+):\s+Resource usage collected/ ) 
        {
            if ( !exists($$job{started}) ) { confess("No wall time for job $$job{id}??", @lines); }
            my $wall_time = DateTime->new(month=>$months{$1}, day=>$2, hour=>$3, minute=>$4, year=>$year)->epoch - $$job{started};
            if ( !exists($$job{cpu_time}) or $$job{cpu_time} < $wall_time ) { $$job{cpu_time} = $wall_time; }
        }
        if ( $lines[$i]=~/The CPU time used is (\d+) seconds./ ) 
        { 
            if ( !exists($$job{cpu_time}) or $$job{cpu_time} < $1 ) { $$job{cpu_time} = $1; }
        }
        if ( $lines[$i]=~/started on (\d+) Hosts\/Processors/ ) 
        { 
            $$job{cpus} = $1;
        }
        if ( $lines[$i]=~/Exited with exit code (\d+)\./ ) 
        { 
            $$job{exit_code} = $1;
        }
    }
    if ( scalar keys %$job) 
    { 
        if ( $$job{command}=~/^cr_restart/ && exists($$job{exit_code}) && $$job{exit_code} eq '16' )
        {
            # temporary failure (e.g. pid in use) of cr_restart, ignore this failure
            return $info;
        }
        $$info{$$job{id}} = $job; 
    }
    return $info;
}

sub _check_job
{
    my ($self,$job,$jids_file) = @_;
    my $status = $$self{lsf_status_codes};
    if ( !exists($$status{$$job{status}}) ) 
    { 
        confess("Todo: $$job{status} $$job{lsf_id}\n"); 
    }
    $$job{status} = $$status{$$job{status}};
    if ( $$job{status}==$$self{Running} )
    {
        if ( !exists($$job{cpu_time}) ) { $$job{cpu_time} = 0; }

        # Estimate how long it might take before we are called again + plus 5 minutes to be safe, and
        # bswitch to a longer queue if necessary.

        my $wakeup_interval = $$self{limits}{wakeup_interval} ? $$self{limits}{wakeup_interval} + 300 : 300;
        my $time_mins = ($$job{cpu_time} / $$job{cpus} + $wakeup_interval) / 60.;
        my $new_queue = $self->_get_queue($time_mins);
        my $cur_queue = $$job{queue};
        if ( defined $new_queue && $new_queue ne $cur_queue && $$self{queue_limits}{$new_queue} > $$self{queue_limits}{$cur_queue} )
        {
            warn("Switching job $$job{lsf_id} from queue $cur_queue to $new_queue\n");
            `bswitch $new_queue '$$job{lsf_id}'`;
            if ( $? ) { warn("Could not switch queues: $$job{lsf_id}"); }
            else { $$job{queue} = $new_queue; }
        }
    }
}

# time is in minutes
sub _get_queue
{
    my ($self,$time) = @_;
    my $queue = exists($$self{limits}{queue}) ? $$self{limits}{queue} : $$self{default_limits}{queue};
    if ( $$self{queue_limits}{$queue} >= $time ) { return $queue; }
    for my $q (sort {$$self{queue_limits}{$a} <=> $$self{queue_limits}{$b}} keys %{$$self{queue_limits}})
    {
        if ( $time > $$self{queue_limits}{$q} ) { next; }
        return $q;
    }
    return undef;
}

sub _parse_output
{
    my ($self,$jid,$output) = @_;

    my $fname = "$output.$jid.o";
    if ( !-e $fname ) { return undef; }
    
    # if the output file is empty, assume the job is running
    my $out = { status=>$$self{Running} };

    # collect command lines and exit status to detect non-critical
    # cr_restart exits
    my @attempts = ();

    open(my $fh,'<',$fname) or confess("$fname: $!");
    while (my $line=<$fh>)
    {
        # Subject: Job 822187: <_2215_1_graphs> Done
        if ( $line =~ /^Subject: Job.+\s+(\S+)$/ )
        {
            if ( $1 eq 'Exited' ) { $$out{status} = $$self{Error}; $$out{nfailures}++; }
            if ( $1 eq 'Done' ) { $$out{status} = $$self{Done}; $$out{nfailures} = 0; }
        }
        if ( $line =~ /^# LSBATCH:/ )
        {
            $line = <$fh>;
            my $cmd = substr($line,0,10);
            push @attempts, { cmd=>$cmd };
            next;
        }
        if ( $line =~ /^Exited with exit code (\d+)\./ )
        {
            if ( !scalar @attempts or exists($attempts[-1]{exit}) ) { warn("Uh, unable to parse $output.$jid.o\n"); next; }
            $attempts[-1]{exit} = $1;
        }
        # Do not count checkpoint and owner kills as a failure. 
        if ( $line =~ /^TERM_CHKPNT/ && $$out{nfailures} ) { $$out{nfailures}--; }
        if ( $line =~ /^TERM_OWNER/ && $$out{nfailures} ) { $$out{nfailures}--; }
    }
    close($fh);
    for (my $i=0; $i<@attempts; $i++)
    {
        # cr_restart exited with a non-critical error
        if ( $attempts[$i]{cmd} eq 'cr_restart' && exists($attempts[$i]{exit}) && $attempts[$i]{exit} eq '16' )
        {
            $$out{nfailures}--;
        }
    }
    return $out;
}

sub past_limits
{
    my ($self,$jid,$output) = @_; 
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

sub _set_lsf_limits_unit
{
    my ($self) = @_;
    if ( exists($$self{lsf_limits_unit}) ) { return; }

    for (my $i=2; $i<15; $i++)
    {
        my @units = grep { /LSF_UNIT_FOR_LIMITS/ } `lsadmin showconf lim 2>/dev/null`;
        if ( $? ) 
        { 
            # lasdmin may be temporarily unavailable and return confusing errors:
            # "Bad host name" or "ls_gethostinfo(): A socket operation has failed: Address already in use"
            print STDERR "lsadmin failed, trying again in $i sec...\n";
            sleep $i; 
            next; 
        }
        if ( @units && $units[0]=~/\s+MB$/ ) { $$self{lsf_limits_unit} = 'MB'; }
        else { $$self{lsf_limits_unit} = 'kB'; }
        return $$self{lsf_limits_unit};
    }
    confess("lsadmin showconf lim failed repeatedly");
}

sub _create_bsub_opts_string
{
    my ($self) = @_;

    # Set bsub options. By default request 1GB of memory, the queues require mem to be set explicitly
    my $bsub_opts = '';
    my $mem    = $$self{limits}{memory} ? int($$self{limits}{memory}) : $$self{default_limits}{memory};
    my $lmem   = $$self{lsf_limits_unit} eq 'kB' ? $mem*1000 : $mem;

    my $runtime = $$self{limits}{runtime} ? $$self{limits}{runtime} : $$self{default_limits}{runtime};
    my $queue   = $self->_get_queue($runtime);
    if ( !defined $queue ) { $queue = $$self{default_limits}{queue}; }

    $bsub_opts  = sprintf " -M%d -R 'select[type==X86_64 && mem>%d] rusage[mem=%d]'", $lmem,$mem,$mem; 
    $bsub_opts .= " -q $queue";
    if ( defined($$self{limits}{cpus}) ) 
    {
        $bsub_opts .= " -n $$self{limits}{cpus} -R 'span[hosts=1]'";
    }
    return $bsub_opts;
}

sub _create_bsub_ids_string
{
    my ($self,$job_name,$ids) = @_;

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
    my @skipped_bsub_ids;
    while ( length($job_name) + length($bsub_ids) > 250 && scalar @bsub_ids ) 
    {
        push @skipped_bsub_ids, pop(@bsub_ids);
        $bsub_ids = join(',', @bsub_ids);
    }
    @$ids = ();
    foreach my $bsub_id (@skipped_bsub_ids)
    {
        if ($bsub_id =~ m/(\d+)-(\d+)/) { push @$ids, ($1..$2); }
        else { push @$ids, $bsub_id; }
    }
    return $bsub_ids;
}

sub _bsub_command
{
    my ($self,$jids_file,$job_name,$bsub_cmd,$cmd) = @_;

    print STDERR "$bsub_cmd\n";
    my @out = `$bsub_cmd`;
    if ( scalar @out!=1 || !($out[0]=~/^Job <(\d+)> is submitted/) )
    {
        my $cwd = `pwd`;
        confess("Expected different output from bsub. The command was:\n\t$cmd\nThe bsub_command was:\n\t$bsub_cmd\nThe working directory was:\n\t$cwd\nThe output was:\n", @out);
    }

    # Write down info about the submitted command
    my $jid = $1;
    open(my $jids_fh, '>>', $jids_file) or confess("$jids_file: $!");
    print $jids_fh "$jid\t$job_name\t$bsub_cmd\n";
    close $jids_fh;
}

sub run_jobs
{
    my ($self,$jids_file,$job_name,$cmd,$ids) = @_;

    if ( !scalar @$ids ) { confess("No IDs given??\n"); }

    $cmd =~ s/{JOB_INDEX}/\$LSB_JOBINDEX/g;
    my $bsub_opts = $self->_create_bsub_opts_string();

    my @ids = @$ids;
    while ( @ids )
    {
        my $bsub_ids = $self->_create_bsub_ids_string($job_name,\@ids);

        # Do not allow the system to requeue jobs automatically, we would loose track of the job ID: -rn 
        my $bsub_cmd  = qq[bsub -rn -J '${job_name}[$bsub_ids]' -e $job_name.\%I.e -o $job_name.\%I.o $bsub_opts '$cmd'];

        # Submit to LSF
        $self->_bsub_command($jids_file,$job_name,$bsub_cmd,$cmd);
    }
}

1;

