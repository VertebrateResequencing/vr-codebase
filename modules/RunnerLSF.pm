package RunnerLSF;

use strict;
use warnings;
use Carp;
use DateTime;
use File::Spec;
use Data::Dumper;

our $Running = 1;
our $Error   = 2;
our $Unknown = 4;
our $No      = 8;
our $Done    = 16;

sub is_job_array_running
{
    my ($jids_file, $ids, $nmax) = @_;

    my @jobs = ();
    for my $id (@$ids) { push @jobs, {status=>$No, nfailures=>0, last_checkpoint_lsf_id=>""}; }
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
    close $fh or confess("$jids_file: $!");
    for (my $i=@jids-1; $i>=0; $i--)
    {
        my $info = parse_bjobs_l($jids[$i]);
        if ( !defined $info ) { next; }
        for my $job (values %$info)
        {
            check_job($job, $jids_file);
        }
        for (my $j=0; $j<@$ids; $j++)
        {
            my $id = $$ids[$j];
            if ( !exists($$info{$id}) ) { next; }
	    if ( $jobs[$j]{last_checkpoint_lsf_id} eq "" && !exists($$info{$id}{chkpnt_failed}) ) 
	    {
		#$jobs[$j]{last_checkpoint_lsf_id} = $$info{$id}{lsf_id};
                $jobs[$j]{last_checkpoint_lsf_id} = $$info{$id}{term_chkpnt};
		$jobs[$j]{queue} = $$info{$id}{queue};
		$jobs[$j]{mem} = $$info{$id}{mem};
		$jobs[$j]{checkpoint_dir} = $$info{$id}{checkpoint_dir};
		$jobs[$j]{name} = $$info{$id}{name};
		$jobs[$j]{chkpnt_failed} = $$info{$id}{chkpnt_failed};
	    }
            if ( $jobs[$j]{status} ne $No ) { next; }   # the job was submitted multiple times and already has a status
            if ( $$info{$id}{status}==$Done ) { $jobs[$j]{status} = $Done; }
            if ( $$info{$id}{status}==$Running ) { $jobs[$j]{status} = $Running; }
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
#	    $jobs[$i]{last_checkpoint_lsf_id} = ""; 
#	    $jobs[$i]{queue} = ""; 
#	    $jobs[$i]{mem} = ""; 
#	    $jobs[$i]{checkpoint_dir} = ""; 
        }
    }
    return \@jobs;
}

sub parse_bjobs_l
{
    my ($jid) = @_;

    my @job_output_sections;
    # FIXME: NO IDEA WHY THIS LOOP IS HERE
#    for (my $i=0; $i<3; $i++)
#    {
    my $output = `bjobs -l $jid 2>/dev/null`;
    if ( $? ) { sleep 5; next; }
    # roll up all continuation lines, signified by beginning with 21 spaces 
    # into a single (long) line beginning with the line before
    $output =~ s/\r//gm;
    $output =~ s/\n\s{21}//gm;
    # split bjobs -l output into separate entries for each job / job array element
    @job_output_sections = split /^[-]+[[:space:]]*$/m, $output; 
    if ( !scalar @job_output_sections ) { return undef; }
#    }

    my $info = {};
    foreach my $job_section (@job_output_sections)
    {
	my $job = parse_bjobs_l_section($job_section);
	if ( scalar keys %$job ) { $$info{$$job{id}} = $job; }
    }

    return $info;
}

sub parse_bjobs_l_section 
{
    my ($jobinfo) = @_;

    my $quoted_entry_capture_regex = '<(.*?)>(,|;|$)[[:space:]]*';
    my %months = qw(Jan 1 Feb 2 Mar 3 Apr 4 May 5 Jun 6 Jul 7 Aug 8 Sep 9 Oct 10 Nov 11 Dec 12);
    my $year = (gmtime())[5] + 1900;

    my $job = {};
    foreach my $line (split /\n/, $jobinfo)
    {
        if ( $line =~ /^\s*$/ ) { next; }
	
	# "Job" line (continuation lines have already been collapsed)
        if ( $line =~ /^Job <(\d+)(.*)$/ )
        {
	    if ( !($line =~ s/Job $quoted_entry_capture_regex//) ) { confess("Could not determine the job id: [$line]"); }
	    $$job{lsf_id} = $1;

	    if ( !($line =~ s/Job Name $quoted_entry_capture_regex//) ) { confess("Could not determine the job name: [$line]"); }
	    $$job{full_name} = $1;

            # Runner's ID is $id, LSF job ID is lsf_id
            $$job{id} = $$job{lsf_id};
            if ( $$job{lsf_id} =~ /\[(\d+)\]/ )
            {
                $$job{id} = $1;
            } 
	    # if there was no index in the lsf id, get it from the job name
	    # (when jobs are brestart-ed from a checkpoint, they are no longer actually job arrays - thanks LSF)
	    elsif ( $$job{full_name} =~ /\[(\d+)\]/ )
	    {
		$$job{id} = $1;
	    }
	    
	    if( $$job{full_name} =~ /^(.*)\[.*?$/ ) 
	    {
		# job name only has first part of name (without job array brackets)
		$$job{name} = $1;
	    }

            if ( !($line =~ s/Status $quoted_entry_capture_regex//) ) { confess("Could not determine the status: [$line]"); }
            $$job{status} = $1;
	    
            if ( !($line =~ s/Queue $quoted_entry_capture_regex//) ) { confess("Could not determine the queue: [$line]"); }
            $$job{queue} = $1;
        }

	# "Submitted" line (again, already including continuation lines, if any)
        if ( $line =~ /^\w+\s+(\w+)\s+(\d+) (\d+):(\d+):(\d+): [Ss]ubmitted/ )
        {
            if ($line =~ s/Checkpoint directory $quoted_entry_capture_regex//)
	    {
		my $checkpoint_dir = $1;
                my $chkpnterr_file = $checkpoint_dir."/echkpnt.err";
		my @dirs = File::Spec->splitdir($checkpoint_dir);
		#pop @dirs;
		$$job{checkpoint_dir} = File::Spec->catdir( @dirs );
                if ( -s $chkpnterr_file && $$job{status} eq "RUN")
                {
                    $$job{chkpnt_failed} = 1; 
                    warn( "Checkpoint error log found\n\n");
                }
                else
                {
	           $$job{term_chkpnt} = $$job{lsf_id};
                }
	    }

	    if ($line =~ s/Requested Resources $quoted_entry_capture_regex//)
	    {
		$$job{resources} = $1;

		if ( $$job{resources} =~ m/rusage\[[^]]*?mem=(\d+)[^]]*?\]/ ) 
		{
		    $$job{mem} = $1;
		}
	    }
	    else
	    {
		confess("Could not find Requested Resources: [$line]");
	    }
        }

        # Tue Mar 19 13:00:35: [685] started on <uk10k-4-1-07>...
        # Tue Dec 24 13:12:00: [1] started on 8 Hosts/Processors <8*vr-1-1-05>...
        if ( $line =~ /^\w+\s+(\w+)\s+(\d+) (\d+):(\d+):(\d+):.+[Ss]tarted on/ ) 
        {
            $$job{started} = DateTime->new(month=>$months{$1}, day=>$2, hour=>$3, minute=>$4, year=>$year)->epoch;
        }

        # Tue Mar 19 13:58:23: Resource usage collected...
        if ( $line =~ /^\w+\s+(\w+)\s+(\d+) (\d+):(\d+):(\d+):\s+Resource usage collected/ ) 
        {
            if ( !exists($$job{started}) ) { confess("No wall time for job $$job{id}??", $line); }
            my $wall_time = DateTime->new(month=>$months{$1}, day=>$2, hour=>$3, minute=>$4, year=>$year)->epoch - $$job{started};
	    # FIXME: this seems broken - CPU time is being forced to be at least as long as wall time???
            if ( !exists($$job{cpu_time}) or $$job{cpu_time} < $wall_time ) { $$job{cpu_time} = $wall_time; }
        }

        # Tue Feb 11 21:38:14: Completed <exit>; TERM_CHKPNT: job killed after checkpointing. Catch last Job_Id checkpointed
	if ( $line =~ /TERM_CHKPNT/ && !$$job{term_chkpnt})
	{
	    if ( $line =~ /Exited with exit code 139/ )
	    {
		$$job{chkpnt_failed} = 1;
	    }
	    $$job{term_chkpnt} = $$job{lsf_id};
	}
        # Job was suspended by the user while pending; Checkpointing failed so fail job so it can be restarted
        if ( $line =~ /Job was suspended by the user while pending\;/ )
        {
            confess(
		"A job was found in a suspended state (supposedly by the user but possibly automatically due to repeated failed attempts to restart from checkpoint).\n"
		"The pipeline cannot proceed. " .
		"If you suspended this job, resume it and retry.\n" .
		"Otherwise, please fix manually by deleting the appropriate line in .jobs/*.w.jid\n" .
		"and removing the appropriate LSF output file(s) .jobs/*.w.<ID>.o\n"
		);
        }

        if ( $line =~ /The CPU time used is (\d+) seconds./ ) 
        { 
	    # FIXME: again, this seems to imply that cpu time will only be set to the actual cpu time reported if it is greater than the previous (wall-time based) setting
            if ( !exists($$job{cpu_time}) or $$job{cpu_time} < $1 ) { $$job{cpu_time} = $1; }
        }
    }
    
    return $job;
}

sub check_job
{
    my ($job) = @_;
    my $status = { DONE=>$Done, PEND=>$Running, EXIT=>$Error, RUN=>$Running, UNKWN=>$Running, SSUSP=>$Running, PSUSP=>$Running };
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
        #my %queue_limits = ( basement=>1e9, long=>2880*60, normal=>720*60, small=>30*60 );
        my %queue_limits = ( basement=>1e9, long=>1*60*60, normal=>1*60*60, small=>30*60 );
        if ( exists($queue_limits{$queue}) && $$job{cpu_time} > $queue_limits{$queue} )
        {
            if ( exists($$job{checkpoint_dir}) )
            {
                # we will just leave checkpointed jobs running till they reach RUNLIMIT
                last;
            }
=head
            {
		# use checkpoint/restart to requeue job
		warn("Checkpointing (and then killing) job $$job{lsf_id} as it is running into runlimit of $$job{queue} queue ($$job{cpu_time} > $queue_limits{$queue})\n");
		# if this section runs again before the checkpoint is complete, no damage is done as `bchkpnt -k` will just 
		# report it is already being checkpointed
		`bchkpnt -k '$$job{lsf_id}'`;
		# once the job terminates, we should find it in term_chkpnt state and attempt to brestart it
            }
=cut
	    else # !exists($$job{checkpoint_dir})
	    {
		# use bswitch instead of checkpoint/restart, for some reason only after 1.3x the cpu time limit has gone
		if ( $$job{cpu_time}*1.3 > $queue_limits{$queue} )
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
    }
}


sub parse_output
{
    my ($jid,$output) = @_;

    my $fname = "$output.$jid.o";
    if ( !-e $fname ) { return undef; }
    
    # if the output file is empty, assume the job is running
    my $out = { status=>$Running };

    open(my $fh,'<',$fname) or confess("$fname: $!");
    while (my $line=<$fh>)
    {
        # Subject: Job 822187: <_2215_1_graphs> Done
        if ( $line =~ /^Subject: Job.+\s+(\S+)$/ )
        {
print "\n\n$line\n\n";
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

our $lsf_limits_unit;
sub get_lsf_limits_unit
{
    if ( defined $lsf_limits_unit ) { return $lsf_limits_unit; }
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
        if ( @units && $units[0]=~/\s+MB$/ ) { $lsf_limits_unit = 'MB'; }
        else { $lsf_limits_unit = 'kB'; }
        return $lsf_limits_unit;
    }
    confess("lsadmin showconf lim failed repeatedly");
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
    my @skipped_bsub_ids;
    while ( length($job_name) + length($bsub_ids) > 250 && scalar @bsub_ids ) 
    {
        push @skipped_bsub_ids, pop(@bsub_ids);
        $bsub_ids = join(',', @bsub_ids);
    }

  
    $cmd =~ s/{JOB_INDEX}/\$LSB_JOBINDEX/g;
    my $bsub_opts = '';
    my $mem;
    if ( exists($$opts{memory}) && $$opts{memory} ) 
    { 
        $mem = int($$opts{memory});
        if ( $mem )
        {
            my $units = get_lsf_limits_unit();
            my $lmem  = $units eq 'kB' ? $mem*1000 : $mem;
            $bsub_opts = sprintf " -M%d -R 'select[type==X86_64 && mem>%d] rusage[mem=%d]'", $lmem,$mem,$mem; 
        }
    }
    if ( !defined($$opts{queue}) ) 
    {            
        if ( !$$opts{_run_with_blcr} && defined($$opts{runtime}) ) 
        { 
            if ( $$opts{runtime} <= 720.0 ) { $$opts{queue} = 'normal'; }
            elsif ( $$opts{runtime} <= 60*24*2 ) { $$opts{queue} = 'long'; }
            else { $$opts{queue} = 'basement'; }
        }
        else 
        { 
            $$opts{queue} = 'normal';
        }
    }
    if ( defined($$opts{queue}) ) 
    {
        $bsub_opts .= " -q $$opts{queue}";
        if ( $$opts{queue} eq 'normal' ){ $$opts{chkpnt_period} = 350; }
        elsif ( $$opts{queue} eq 'long' ){ $$opts{chkpnt_period} = 700; }
        elsif ( $$opts{queue} eq 'basement' ){ $$opts{chkpnt_period} = 1440; }
        else { $$opts{chkpnt_period} = 350; }
    }
    if ( defined($$opts{cpus}) ) 
    {
        $bsub_opts .= " -n $$opts{cpus} -R 'span[hosts=1]'";
    }
    if ( $$opts{_run_with_blcr} ) 
    {
	my $cwd = `pwd`;
	chomp $cwd;
	$jids_file =~ /(.*\/\.jobs\/).*jid/;
	my $temp_dir = $1;
	my $checkpoint_dir = $cwd."/".$temp_dir."checkpoint";
        if ( !$$opts{chkpnt_period} ) { $$opts{chkpnt_period} = 350; }
        #$bsub_opts .= " -k '$checkpoint_dir method=blcr $$opts{chkpnt_period}'";
        $bsub_opts .= " -k '$checkpoint_dir method=blcr 15'";
        $cmd = "cr_run $cmd";
    }

    my $bsub_cmd  = qq[bsub -J '${job_name}[$bsub_ids]' -e $job_name.\%I.e -o $job_name.\%I.o $bsub_opts '$cmd'];

    # Submit to LSF
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

    # If we skipped any subs due to command line limit, call run_array again
    if (@skipped_bsub_ids)
    {
        my @skipped_ids;
        foreach my $bsub_id (@skipped_bsub_ids)
        {
            if ($bsub_id =~ m/(\d+)-(\d+)/) { push @skipped_ids, ($1..$2); }
            else { push @skipped_ids, $bsub_id; }
        }
        run_array($jids_file, $job_name, $opts, $cmd, \@skipped_ids);
    }
}

sub restart_job
{
    my ($jids_file, $job, $opts) = @_;

    my $brestart_opts = "";
    
    if ( defined($$job{queue}) )
    {
	$brestart_opts .= "-q $$job{queue} ";
    }
    
    if ( defined($$job{mem}) )
    {
	my $units = get_lsf_limits_unit();
	my $lmem = $units eq 'kB' ? ($$job{mem}*1000) : $$job{mem};
	$brestart_opts .= "-M $lmem ";
    }
    
    my $brestart_cmd = qq[brestart $brestart_opts  '$$job{checkpoint_dir}' '$$job{last_checkpoint_lsf_id}'];
    
    if ( !defined('$$job{name}') )
    {
	confess("Job name not defined for $$job{last_checkpoint_lsf_id}");
    }
    
    # Submit to LSF
    print STDERR "$brestart_cmd\n";
    my @out = `$brestart_cmd`;
    if ( scalar @out!=1 || !($out[0]=~/^Job <(\d+)> is submitted/) )
    {
	confess("Expected different output from brestart.\nThe brestart_command was:\n\t$brestart_cmd\nThe output was:\n", @out);
    }
    # Checkpoint log is not found or is corrupted. Job not submitted.
    
    # Write down info about the submitted command
    my $jid = $1;
print "Job name set to $$job{name}\n";
    open(my $jids_fh, '>>', $jids_file) or confess("$jids_file: $!");
    print $jids_fh join("\t", $jid, $$job{name} , $brestart_cmd)."\n";
    close $jids_fh;
    
}

1;

