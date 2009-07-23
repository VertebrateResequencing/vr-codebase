package LSF;

use strict;
use warnings;
use Utils;

our $No      = 0;
our $Running = 1;
our $Error   = 2;
our $Unknown = 4;


=pod

=head1 NAME

LSF

=head1 SYNOPSIS

Utilities for manipulating LSF jobs.

=head1 METHODS

=head2 is_job_running

    Arg [1]     : lock file containing LSF job ID's.
    Returntype  : $LSF::Running if at least one of the jobs is still in the queue and/or running.
                  $LSF::No if none of the jobs is present in the queue.
                  $LSF::Error if some of the jobs failed.
                  If the lock file is empty, $LSF::No is returned.
                  If some of the jobs failed while others are running, $LSF::Running|$LSF::Error is returned.

=cut

sub is_job_running
{
    my ($jids_file) = @_;

    if ( ! -e $jids_file ) { return $No }

    my $job_running = $No;

    open(my $fh, '<', $jids_file) or Utils::error("$jids_file: $!");
    while (my $jid=<$fh>)
    {
        chomp($jid);
        # For backwards compatibility, parse either a single integer or an integer
        #   followed by \t and path to a LSF output file.
        if ( !($jid=~/^(\d+)\s*(\S*.*)$/) ) { Utils::error("Uh, could not parse \"$jid\".\n") }

        my $status = job_in_queue($1,$2);
        if ( $status == $Error ) { $job_running |= $Error; }
        if ( $status == $Running ) 
        { 
            # Do not exit if one running was found - check all jobs, no one should have the EXIT status.
            $job_running |= $Running; 
        }
    }
    close($fh);

    return $job_running;
}


# Return status: $No not in queue, $Running running, $Error aborted
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
        my $bt = Utils::backtrace();
        print STDERR @$bt,"\nFIXME: $jid .. $lsf_output: $!\n"; 
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
            if ( $1 eq 'Done' ) { $status=$No;  }
            last;
        }
    }
    close($fh);
    return $status;
}


sub job_in_bjobs
{
    my ($jid) = @_;

    # Expexting either
    #   JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
    #   735850  pd3     DONE  normal     sf-2-1-03   sf-1-4-14   mouse100   May 21 09:50
    #
    # on STDOUT or
    #   Job <735851> is not found
    #
    # on STDERR.

    my @out = Utils::CMD("bjobs $jid 2>/dev/null");
    
    if ( ! scalar @out ) { return $Unknown; }
    if ( scalar @out != 2 ) { Utils::error("Expected different output, got: ", @out) }
    if ( $out[1] =~ /^$jid\s+\S+\s+DONE/ ) { return $No; }
    if ( $out[1] =~ /^$jid\s+\S+\s+EXIT/ ) { return $Error; }

    return $Running;
}


=head2 run

    Arg [1]     : lock file where to put the JID (opened in append mode)
    Arg [2]     : before bsub will be called, chdir to the working dir
    Arg [3]     : job name - the output will be redirected to $job.o and $job.e
    Arg [4]     : options, expecting 'bsub_opts' key in the hash with the command line options
    Arg [5]     : the command to run
    Description : Executes bsub with the given parameters.
    Returntype  : none

=cut

sub run
{
    my ($jids_file,$work_dir,$job_name,$options,$bsub_cmd) = @_;

    if ( !exists($$options{'bsub_opts'}) ) { Utils::error("No 'bsub_opts' given.\n") }
    my $bsub_opts = $$options{'bsub_opts'};

    chomp(my ($cwd) = Utils::CMD("pwd"));
    if ( $work_dir ) { chdir($work_dir) or Utils::error("chdir \"$work_dir\": $!") }

    # The expected output:
    #   Job <771314> is submitted to queue <normal>.
    my $lsf_output_file = "$job_name.o";
    my $cmd = "bsub -J $job_name -e $job_name.e -o $lsf_output_file $bsub_opts '$bsub_cmd'";

    my @out = Utils::CMD($cmd,$options);
    if ( scalar @out!=1 || !($out[0]=~/^Job <(\d+)> is submitted/) ) 
    { 
        Utils::error("Expected different output from bsub. The command was:\n\t$cmd\nThe output was:\n", @out);
    }
    open(my $jids_fh, '>>', $jids_file) or Utils::error("$jids_file: $!");
    print $jids_fh "$1\t$work_dir/$lsf_output_file\n";
    close $jids_fh;

    if ( $work_dir ) { chdir($cwd) or Utils::error("chdir \"$cwd\": $!"); }

    return;
}


1;

