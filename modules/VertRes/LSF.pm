package VertRes::LSF;

use strict;
use warnings;
use Utils;
use File::Basename;
use File::Spec;
use VertRes::Parser::LSF;

our $Running = 1;
our $Error   = 2;
our $Unknown = 4;
our $No      = 8;
our $Done    = 16;
our $Max_jobs_to_check = 10;


=pod

=head1 NAME

VertRes::LSF

=head1 SYNOPSIS

Utilities for manipulating LSF jobs.

=head1 METHODS

=head2 is_job_running

    Arg [1]     : lock file containing LSF job ID's.
    Returntype  : $LSF::Running if at least one of the jobs is still in the queue and/or running.
                  $LSF::No if none of the jobs is present in the queue and there are no records about the job.
                  $LSF::Error if some of the jobs failed.
                  $LSF::Done if at least one job finished OK (so you should check for $LSF:Error first)
                  If the lock file is empty, $LSF::No is returned.
                  If some of the jobs failed while others are running, $LSF::Running|$LSF::Error is returned.


=cut

sub is_job_running
{
    my ($jids_file) = @_;

    if ( ! -e $jids_file ) { return $No }

    my $job_running = $No;

    open(my $fh, '-|', "tail -n ".$Max_jobs_to_check." ".$jids_file) or confess("$jids_file: $!");    
    while (my $jid=<$fh>)
    {
        chomp($jid);
        # For backwards compatibility, parse either a single integer or an integer
        #   followed by \t, path to a LSF output file, and the command.
        if ( !($jid=~/^(\d+)\s*([^\t]*)/) ) { Utils::error("Uh, could not parse \"$jid\".\n") }

        my $status = job_in_queue($1,$2);
        if ( $status == $Error ) { $job_running |= $Error; }
        if ( $status == $Running ) 
        { 
            # Do not exit if one running was found - check all jobs, no one should have the EXIT status.
            $job_running |= $Running; 
        }
        if ( $status == $Done ) { $job_running |= $Done; }
    }
    close($fh);

    if ( $job_running & (~$No) ) { $job_running = $job_running & (~$No) }

    return $job_running;
}


# Return status: $No not in queue, $Running running, $Error aborted, $Done finished.
sub job_in_queue
{
    my ($jid,$lsf_output) = @_;

    my $status = $Unknown;
    if ( defined $jid ) { $status=job_in_bjobs($jid); }
    if ( $status != $Unknown ) { return $status; }

    if ( !$lsf_output ) { return $No; }
    
    # If the job is not in queue, it must be finished, either with the Done or
    #   Exited status.  Parse the output file and see the status. Any failure 
    #   is non-critical. In the worst scenario, the job will be rerun.
    #
    my $fh;
    $status = $No;
    if ( !open($fh,'<',$lsf_output) ) 
    { 
        # This really can happen: The status of the LSF job can be queried before
        #   the file is created.
        return $status; 
    }
    while (my $line=<$fh>)
    {
        # Subject: Job 822187: <_2215_1_graphs> Done
        if ( $line =~ /^Subject: Job $jid: \<[^>]+\> (\S+)/ )
        {
            # It is unlikely that there had been two jobs with the same IDs' 
            #   but the last wins anyway.
            if ( $1 eq 'Exited' ) { $status=$Error; }
            if ( $1 eq 'Done' ) { $status=$Done;  }
            last;
        }
    }
    close($fh);
    return $status;
}


sub job_in_bjobs
{
    my ($jid) = @_;

    # Expecting either
    #   JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
    #   735850  pd3     DONE  normal     sf-2-1-03   sf-1-4-14   mouse100   May 21 09:50
    #
    # on STDOUT or
    #   Job <735851> is not found
    #
    # on STDERR.

    my @out;
    for (my $i=0; $i<3; $i++)
    {
        @out = `bjobs -l $jid 2>/dev/null`;
        if ( $? ) { sleep 5; next; }
        if ( !scalar @out ) { return $Unknown; }
    }

    my $job = parse_bjobs_l(\@out);
    if ( $$job{status} eq 'DONE' ) { return $Done; }
    if ( $$job{status} eq 'PEND' ) { return $Running; }
    if ( $$job{status} eq 'EXIT' ) { return $Error; } 
    if ( $$job{status} eq 'RUN' ) 
    {
        my $queue = $$job{queue};
        my %queue_limits = ( basement=>1e9, long=>2880*60, normal=>720*60, small=>30*60 );
        if ( exists($queue_limits{$queue}) && $$job{cpu_time}*1.3 > $queue_limits{$queue} )
        {
            my $bswitch;
            for my $q (sort {$queue_limits{$a} <=> $queue_limits{$b}} keys %queue_limits)
            {
                if ( $$job{cpu_time}*1.3 > $queue_limits{$q} ) { next; }
                warn("Changing queue of the job $jid from $$job{queue} to $q\n");
                `bswitch $q $jid`;

            }
        }
        return $Running;
    }
    return $Error;
}

sub parse_bjobs_l
{
    my ($lines) = @_;

    my $i=0;
    while ( $i<@$lines && $$lines[$i]=~/^\s*$/ ) { $i++ }
    if ( $i>=@$lines ) { Utils::error("Could not parse bjobs -l output: ", join('',@$lines)); }

    my $job_info = $$lines[$i++]; chomp($job_info);
    while ( $i<@$lines && $$lines[$i]=~/^\s{21}?/ ) 
    { 
        $job_info .= $';
        chomp($job_info);
        $i++;
    }

    if ( !($job_info=~/,\s*Status <([^>]+)>/) ) { Utils::error("Could not determine the status: [$job_info]"); }
    my $status = $1;
    if ( !($job_info=~/,\s*Queue <([^>]+)>/) ) { Utils::error("Could not determine the queue: [$job_info]"); }
    my $queue = $1;

    my $cpu_time = 0;
    while ( $i<@$lines )
    {
        if ( $$lines[$i]=~/The CPU time used is (\d+) seconds./ ) { $cpu_time=$1; last; }
        $i++;
    }
    return { status=>$status, queue=>$queue, cpu_time=>$cpu_time };
}


# Check if the output file exists. If it is not there, return unchanged
#   options. If the job failed due to memory or time, increase the requested
#   memory or change queue.
#
# Performance: stat on non-existent files is very fast. Only failed jobs will slow
#   things down. If it proves to be an issue, DB will be used instead.
#
sub adjust_bsub_options
{
    my ($opts, $output_file) = @_;

    my $mem;
    my $queue;
    
    my $orig_mem = 0;
    my $no_warn = 0;
    if ( !($opts=~/-M/) ) { $mem = -500; $no_warn = 1; }  # if no mem specified, force a reservation of 500MB
    elsif ( $opts =~ /-M(\d+)/ )
    {
        $orig_mem = $1*1e-3;
    }
    elsif ($opts =~ /select\[mem[^\]\d]+(\d+)/) {
        $orig_mem = $1;
    }
    
    # some pipelines make use of VertRes::Pipeline::archive_bsub_files, which
    # moves the output_file to two different possible files, so we need to check
    # all 3 possible output_files
    my $prev = $output_file.'.previous';
    my ($base, $path) = fileparse($output_file);
    my $arch = File::Spec->catfile($path, '.'.$base.'.archive');
    foreach my $o ($output_file, $prev, $arch) {
        if ( -e $o ) 
        {
            my $parser = VertRes::Parser::LSF->new(file=>$o);
            my $n = $parser->nrecords() || 0;
            
            for (my $i=0; $i<$n; $i++)
            {
                my $status = $parser->get('status',$i) || next;
                # Find out the reason of the failure. If the memory was too low, increase it.
                if ($status  eq 'MEMLIMIT')
                {
                    if ( !defined $mem or $mem<$parser->get('memory',$i) ) { $mem=$parser->get('memory',$i); }
                }
                elsif ($status eq 'RUNLIMIT') 
                {
                    $queue = $parser->get('queue',$i);
                    if ( $queue eq 'normal' ) { $queue='long' }
                    else { $queue = 'basement'; }
                }
            }
        }
    }
    
    if (defined $mem) {
      # at some point an attempt to run this failed due to MEMLIMIT
        $mem = calculate_memory_limit($mem);
        
        if ($mem>500000) {
            Utils::error("FIXME: This job cannot be run on the farm, more than 500GB of memory is required.");
        }
        elsif($mem>250000)
        {
                warn("$output_file: changing queue to hugemem\n");   
                $opts =~ s/-q normal/-q hugemem/;
                $opts =~ s/-q long/-q hugemem/;
                if ( !($opts=~/-q/) ) { $opts .= ' -q hugemem'; }
        }
        
        # adjust the command line to include higher memory reservation, but only
        # if the user hasn't since increased the memory manually to higher than
        # our +1000 value
        if ($mem > $orig_mem) {
            # The kind of option line we are trying to produce
            #   -q normal -M3000 -R 'select[type==X86_64 && mem>3000] rusage[mem=3000]'
            warn("$output_file: Increasing memory to $mem\n") unless $no_warn;  # this should be logged in the future
            
            $opts =~ s/-M\d+/-M$mem/;
            $opts =~ s/(select[^]]+mem>)\d+/$1$mem/;
            $opts =~ s/(rusage[^]]+mem=)\d+/$1$mem/;
            
            # The lines above replaced existing values in $opts. If they are not present, add them
            if ( !($opts=~/-M\d+/) ) { $opts .= " -M${mem}" }
            if ( !($opts=~/-R/) ) { $opts .= " -R 'select[mem>$mem] rusage[mem=$mem]'"; }
            else
            {
                if ( !($opts=~/select/) ) { $opts .= " -R 'select[mem>$mem]'"; }
                elsif($opts=~ /select\[mem(.){1,2}[\d]+\]/)
                { 
                  # select[type==X86_64] select[mem>1 ]
                  # If there is more than 1 select statement, dont do anything as the memory has already been modified.
                }
                elsif ( $opts=~/select\[([^]]+)\]/ && !($1=~/mem/) ) { $opts =~ s/select\[/select[mem>$mem &&/ }
                
                # Avoid setting rusage twice
                if(!($opts=~/rusage[^]]+mem=$mem/))
                {
                  if ( !($opts=~/rusage/) ) { $opts .= " -R 'rusage[mem=$mem]'" }
                  elsif ( $opts=~/rusage\[([^]]+)\]/ && !($1=~/mem/) ) { $opts =~ s/rusage\[/rusage[mem=$mem,/ }
                }
            }
        }
    }

    if ( defined $queue )
    {
        warn("$output_file: changing queue to long\n");     # this should be logged in the future

        $opts =~ s/-q normal/-q long/;
        if ( !($opts=~/-q/) ) { $opts .= ' -q long'; }
    }

    return $opts;
}


=head2 past_limits

    Arg [1]     : LSF job name (without the ".o" suffix)
    Description : Find out status and limits of the previous run
    Returntype  : Hash 

=cut

sub past_limits
{
    my ($job_name) = @_; 
    if ( ! -e "$job_name.o" ) { return (); }
    my %out;
    my $parser = VertRes::Parser::LSF->new(file=>"$job_name.o");
    my $n = $parser->nrecords() || 0;
    for (my $i=0; $i<$n; $i++)
    {
        my $status = $parser->get('status',$i) || next;
        my $mem = $parser->get('memory',$i);
        if ( !exists($out{memory}) or $out{memory}<$mem )
        {
            $out{memory} = $mem;
            if ( $status eq 'MEMLIMIT' ) { $out{MEMLIMIT} = $mem; }
            else { delete($out{MEMLIMIT}); }
        }
    }
    return %out;
}


=head2 calculate_memory_limit

    Arg [1]     : memory in mega bytes for previously failed job
    Description : Increases the memory by a set minimum or by a percentage, which ever is greater
    Returntype  : Int

=cut
sub calculate_memory_limit
{
  my ($mem) = @_; 
  my $minimum_memory_increase = 1000;
  my $memory_increase_percentage = 0.3;

  my $updated_memory_limit = $mem*(1+$memory_increase_percentage);
  if($updated_memory_limit < $mem + $minimum_memory_increase)
  {
    $updated_memory_limit =  $mem + $minimum_memory_increase;
  }
  
  return int($updated_memory_limit);
}


=head2 run

    Arg [1]     : lock file where to put the JID (opened in append mode)
    Arg [2]     : before bsub will be called, chdir to the working dir
    Arg [3]     : job name - the output will be redirected to $job.o and $job.e
    Arg [4]     : options
                        append      .. unless explicitly set to zero, the lock file will be opened in the append mode
                        bsub_opts   .. command line options
                        dont_wait   .. if set, don't wait until the job appears in bjobs list
                        memory      .. expected memory requirements (MB) [ignored if bsub_opts present]
                        queue       .. computing farm queue [ignored if bsub_opts present]
                        runtime     .. expected running time (minutes) [ignored if bsub_opts or queue present]
    Arg [5]     : the command to run
    Description : Executes bsub with the given parameters.
    Returntype  : none

=cut

sub run
{
    my ($jids_file,$work_dir,$job_name,$options,$bsub_cmd) = @_;
    
    # The expected output:
    #   Job <771314> is submitted to queue <normal>.
    my $lsf_output_file = "$job_name.o";
    my $lsf_error_file = "$job_name.e";

    my %opts = (%$options);
    if ( !exists($opts{bsub_opts}) )
    {
        if ( !defined($opts{queue}) ) 
        {            
            if( defined($opts{runtime}) ) 
            { 
                if ( $opts{runtime} <= 720.0 ) { $opts{queue} = 'normal'; }
                elsif ( $opts{runtime} <= 60*24*2 ) { $opts{queue} = 'long'; }
                else { $opts{queue} = 'basement'; }
            }
            else 
            { 
                $opts{queue} = 'normal';
            }
        }
        if ( defined($opts{queue}) ) 
        {
            $opts{bsub_opts} = " -q $opts{queue}";
        }
        if ( defined($opts{cpus}) ) 
        {
            $opts{bsub_opts} .= " -n $opts{cpus} -R 'span[hosts=1]'";
        }
        if ( defined($opts{memory}) ) 
        {
            $opts{bsub_opts} .= sprintf " -M%d -R 'select[type==X86_64 && mem>%d] rusage[mem=%d]'", $opts{memory},$opts{memory},$opts{memory};
        }
    }
    if ( !exists($opts{'bsub_opts'}) ) { Utils::error("No 'bsub_opts' given.\n") }
    my $bsub_opts = $opts{'bsub_opts'};

    chomp(my ($cwd) = Utils::CMD("pwd"));
    if ( $work_dir ) { chdir($work_dir) or Utils::error("chdir \"$work_dir\": $!") }

    if($$options{logfile_output_directory})
    {
        my $log_dir = $$options{logfile_output_directory}."/".time().'_'.int( rand(1000));
        if ( !-e $log_dir ) { Utils::CMD("mkdir -p $log_dir"); }
        $lsf_output_file = $log_dir.'/'.$lsf_output_file;
        $lsf_error_file  = $log_dir.'/'.$lsf_error_file;
    }

    # Check if memory or queue should be changed (and change it)
    $bsub_opts = adjust_bsub_options($bsub_opts, $lsf_output_file);
    my $cmd = "bsub -J $job_name -e $lsf_error_file -o $lsf_output_file$bsub_opts '$bsub_cmd'";

    my @out = Utils::CMD($cmd,$options);
    if ( scalar @out!=1 || !($out[0]=~/^Job <(\d+)> is submitted/) ) 
    { 
        Utils::error("Expected different output from bsub. The command was:\n\t$cmd\nThe output was:\n", @out);
    }
    my $jid  = $1;
    my $mode = exists($$options{append}) && !$$options{append} ? '>' : '>>';
    if ( !($lsf_output_file=~m{^/}) ) { $lsf_output_file = "$work_dir/$lsf_output_file"; }
    open(my $jids_fh, $mode, $jids_file) or Utils::error("$jids_file: $!");
    print $jids_fh "$jid\t$lsf_output_file\t$cmd\n";
    close $jids_fh;

    if ( !$$options{dont_wait} )
    {
        # Now wait until the job appears in the queue, it may take few seconds. If another
        #   pipeline was running in the meantime, it would schedule the same job to LSF.
        #   We can safely sleep here, because run_lane protects us by its locking mechanism.
        #
        my $max_wait = 30;
        my $status   = $No;
        while ($max_wait>0)
        {
            $status = job_in_queue($jid,"$lsf_output_file");
            if ( $status!=$No ) { last }
            sleep(1);
            $max_wait-=1;
        }
        if ( $status==$No ) { Utils::error("The job $1 $lsf_output_file still not in queue??\n"); }
    }

    if ( $work_dir ) { chdir($cwd) or Utils::error("chdir \"$cwd\": $!"); }

    return;
}


1;

