#
# Copyright (c) 2013, 2014 Genome Research Ltd. 
# 
# Authors: Petr Danecek <pd3@sanger> 
#          Allan Daly <ad7@sanger.ac.uk>
#          Joshua Randall <jr17@sanger.ac.uk>
#
# This program is free software: you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the Free Software 
# Foundation; either version 3 of the License, or (at your option) any later 
# version. 
#
# This program is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
# details. 
#
# You should have received a copy of the GNU General Public License along with 
# this program. If not, see <http://www.gnu.org/licenses/>.
#

package RunnerLSF;

use strict;
use warnings;
use Carp;
use DateTime;
use File::Spec;
use Data::Dumper;

# TODO: all of these should be farm options that can be overridden in config
my %queue_limits_seconds = ( basement=>1e9, long=>2880*60, normal=>720*60, small=>30*60 );
my $chkpnt_minimum_period_minutes = 1;
my $default_memlimit_mb = 100;

# estimated (minimum) time to checkpoint: ~2.5 minutes per GB of RAM on lustre (scratch113), so 0.15 s/MB
# but it is better to underestimate than overestimate - theoretical maximum is ~500MB/s, so 0.002s/MB
my $min_chkpnt_time_s_per_mb = 0.002;

# maximum allowed time to checkpoint
# lustre can sustain ~10+ GB/s, so 0.0001 s/MB
# loaded by ~10000 jobs, that would be 1 s/MB
my $max_chkpnt_time_s_per_mb = 1;

my %months = qw(Jan 1 Feb 2 Mar 3 Apr 4 May 5 Jun 6 Jul 7 Aug 8 Sep 9 Oct 10 Nov 11 Dec 12);

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
    close $fh or confess("$jids_file: $!");

    my $jidlines_processed = 0;
  JIDLINE: for (my $i=@jids-1; $i>=0; $i--)
  {
      $jidlines_processed++;
      my $info = parse_bjobs_l($jids[$i]);
      if ( !defined $info ) { next JIDLINE; }
      for my $job (values %$info)
      {
	  check_job($job);
      }
    ELEMENT_ID: for (my $j=0; $j<@$ids; $j++)
    {
	my $id = $$ids[$j];

	if ( !exists($$info{$id}) ) { next ELEMENT_ID; }
	# if we don't yet have_chkpnt for this job element, and we do have_chkpnt from this info
	if ( !exists($jobs[$j]{have_chkpnt}) && exists($$info{$id}{have_chkpnt}) ) 
	{
	    # copy all info needed for checkpoint restart from this (the latest) have_chkpnt info
	    $jobs[$j]{have_chkpnt} = $$info{$id}{have_chkpnt};
	    $jobs[$j]{last_chkpnt_lsf_id} = $$info{$id}{lsf_id};
	    $jobs[$j]{queue} = $$info{$id}{queue};
	    $jobs[$j]{mem_mb} = $$info{$id}{mem_mb};
	    $jobs[$j]{chkpnt_dir} = $$info{$id}{chkpnt_dir};
	    $jobs[$j]{name} = $$info{$id}{name};
	    $jobs[$j]{id} = $$info{$id}{id};
	}
	if ( exists($$info{$id}{idle_factor}) )
	{
	    $jobs[$j]{idle_factor} = $$info{$id}{idle_factor};
	}
	if ( exists($$info{$id}{runtime_limit_seconds}) )
	{
	    $jobs[$j]{runtime_limit_seconds} = $$info{$id}{runtime_limit_seconds};
	}
	if ( exists($$info{$id}{cpus}) )
	{
	    $jobs[$j]{cpus} = $$info{$id}{cpus};
	}
	if ( $jobs[$j]{status} ne $No ) { next ELEMENT_ID; }   # the job was submitted multiple times and already has a status
	if ( $$info{$id}{status}==$Done ) { $jobs[$j]{status} = $Done; }
	if ( $$info{$id}{status}==$Running ) { $jobs[$j]{status} = $Running; }
	
    }

      # assume we don't need to process any more jid file lines
      # we will set this to 1 below if an element is still missing needed information
      my $continue_processing_jidlines = 0;
    CHECK_ELEMENTS: for (my $j=0; $j<@$ids; $j++)
    {
	# check if we still need more information for this element
	if ( ( $jobs[$j]{status} != $Done ) &&
	     ( $jobs[$j]{status} != $Running ) &&
	     ( !exists($jobs[$j]{have_chkpnt}) ) )
	{
	    # job is neither Done nor Running and does not have checkpoint information, keep searching through jid file entries
	    $continue_processing_jidlines = 1;
	}
    }	
      # if we don't need to continue processing jid file entries, don't
      if ( ! $continue_processing_jidlines ) 
      { 
	  last JIDLINE; 
      }
    }
    
    if ($jidlines_processed < @jids) 
    {
	my $unprocessed_lines = scalar(@jids) - $jidlines_processed;
	warn(
	    "Only needed to process $jidlines_processed jid file entries out of ".scalar(@jids)."\n" .
	    "You may be able to remove the first $unprocessed_lines lines from the jid file\n"
	    );
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

    my @job_output_sections;

    # loop to retry bjobs up to 3 times if it fails
  BJOBS: for (my $i=0; $i<3; $i++)
  {
      my $output = `bjobs -l $jid 2>/dev/null`;
      if ( $? ) { sleep 5; next BJOBS; }
      # quash carriage returns
      $output =~ s/\r//gm;
      # roll up all continuation lines, signified by beginning with 21 spaces 
      # into a single (long) line beginning with the line before
      $output =~ s/\n\s{21}//gm;
      # roll up *LIMIT or *TIME line with the following line
      $output =~ s/^( [A-Z]+LIMIT| [A-Z]+TIME)\s*\n\s*(.*?)$/$1 $2/gm;
      # split bjobs -l output into separate entries for each job / job array element
      @job_output_sections = split /^[-]+[[:space:]]*$/m, $output; 
      if ( !scalar @job_output_sections ) 
      { 
	  return undef; 
      }
      else 
      {
	  # have some bjobs output
	  last BJOBS;
      }
    }

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
		my $chkpnt_dir = $1;
		my @dirs = File::Spec->splitdir($chkpnt_dir);
		$$job{chkpnt_dir} = File::Spec->catdir( @dirs );
                my $chkpnt_context_file = $$job{chkpnt_dir}."/jobstate.context";
                if ( -e $chkpnt_context_file && -s $chkpnt_context_file )
                {
		    $$job{have_chkpnt} = 1;
                }
	    }

	    if ($line =~ s/Requested Resources $quoted_entry_capture_regex//)
	    {
		$$job{resources} = $1;

		if ( $$job{resources} =~ m/rusage\[[^]]*?mem=(\d+)[^]]*?\]/ ) 
		{
		    $$job{mem_mb} = $1;
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

        # Tue Dec 24 13:12:00: [1] started on 8 Hosts/Processors <8*vr-1-1-05>...
        if ( $line =~ /[Ss]tarted on (\d+) / ) 
        {
            $$job{cpus} = $1;
        }

        # RUNLIMIT 60.0 min of bc-22-4-06
        # RUNTIME 60.0 min of bc-22-4-06
	if ( $line =~ /RUN(LIMIT|TIME) ([0-9.]+) min/ )
	{
	    $$job{runtime_limit_seconds} = $2 * 60;
	}
	
        # Tue Mar 19 13:58:23: Resource usage collected...
        if ( $line =~ /^\w+\s+(\w+)\s+(\d+) (\d+):(\d+):(\d+):\s+Resource usage collected/ ) 
        {
            if ( !exists($$job{started}) ) { confess("No wall time for job $$job{id}??", $line); }
	    $$job{cur_time} = DateTime->new(month=>$months{$1}, day=>$2, hour=>$3, minute=>$4, year=>$year)->epoch;
            my $wall_time = $$job{cur_time} - $$job{started};
            if ( !exists($$job{wall_time}) ) { $$job{wall_time} = $wall_time; }
        }

        # Tue Feb 11 21:38:14: Completed <exit>; TERM_CHKPNT: job killed after checkpointing. Catch last Job_Id checkpointed
	if ( $line =~ /TERM_CHKPNT/ && !$$job{term_chkpnt})
	{
	    if ( $line =~ /exit code 143/ )
	    {
		$$job{term_chkpnt} = $$job{lsf_id};
	    }
	    elsif ( $line =~ /exit code 139/ )
	    {
		warn(
		    "\tc  job $$job{id} ($$job{lsf_id}) appears to have died while checkpointing"
		    );
	    }
	    elsif ( $line =~ /exit code (\d+)/ ) 
	    {
		warn(
		    "\tc  job $$job{id} ($$job{lsf_id}) had unexpected error status while checkpointing (exit code $1)\n"
		    );
	    }
	}
        # Job was suspended by the user while pending;
	# this also sometimes occurs when checkpoint restarts fail repeatedly
        if ( $line =~ /Job was suspended by the user while pending\;/ )
        {
            confess(
		"A job was found in a suspended state (supposedly by the user but possibly automatically due to repeated failed attempts to restart from checkpoint).\n" .
		"The pipeline cannot proceed. " .
		"If you suspended this job, resume it and retry.\n" .
		"Otherwise, please fix manually by deleting the appropriate line in .jobs/*.w.jid\n" .
		"and removing the appropriate LSF output file(s) .jobs/*.w.<ID>.o\n"
		);
        }

        if ( $line =~ /The CPU time used is (\d+) seconds./ ) 
        { 
	    # CPU time can only go up, not down
            if ( !exists($$job{cpu_time}) or $$job{cpu_time} < $1 ) { $$job{cpu_time} = $1; }
        }
	
        if ( $line =~ /IDLE_FACTOR.*:\s*([0-9.]+)/ ) 
        { 
            # IDLE_FACTOR(cputime/runtime):   2.53
	    $$job{idle_factor} = $1; 
        }
    }
    
    # if wall time has not been reported yet, set it to 0
    if ( !exists($$job{wall_time}) ) { $$job{wall_time} = 0; }
    
    # if cpu time has not been reported yet, set it to wall time (as a guess)
    if ( !exists($$job{cpu_time}) ) { $$job{cpu_time} = $$job{wall_time}; }
    
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
        my $queue = $$job{queue};
	my $cpus = 1;
	if ( exists($$job{cpus}) ) { $cpus = $$job{cpus}; }
	my $runtime_limit_seconds = undef;
	if ( exists($$job{runtime_limit_seconds}) )
	{
	    $runtime_limit_seconds = $$job{runtime_limit_seconds};
	}
	elsif ( exists($queue_limits_seconds{$queue}) )
	{
	    $runtime_limit_seconds = $queue_limits_seconds{$queue};
	}
        if ( defined($runtime_limit_seconds) && ( ( ($$job{cpu_time} / $cpus) > $runtime_limit_seconds ) || ( $$job{wall_time} > $runtime_limit_seconds ) ) )
        {
	    warn(
		"\tc  job $$job{id} ($$job{lsf_id}) is currently running and has exceeded the runtime limit of $runtime_limit_seconds seconds (cpu time: $$job{cpu_time} wall time: $$job{wall_time})\n"
		);
	    # we have overrun either cpu time or wall time
            if ( exists($$job{chkpnt_dir}) )
            {
                # we are using checkpointing and are over the queue time limit
		
                # process chkpnt.log to check if a recent checkpoint has been started
		my $chkpnt_log_file = $$job{chkpnt_dir}."/chkpnt.log";
		my @chkpnt_log_lines;
		if ( -f $chkpnt_log_file ) 
		{
		    open(my $fh, '<', $chkpnt_log_file) or confess("$chkpnt_log_file: $!");
		    # read all lines into an array so we can process them in reverse order
		    foreach my $line (<$fh>) 
		    {
			chomp $line;
			push @chkpnt_log_lines, $line;
		    }
		    close $fh or confess("$chkpnt_log_file: $!");
		}
                # Example entries from chkpnt.log:
                # ########### begin to checkpoint ############
                # Fri Feb 21 15:40:13 2014 : Echkpnt : main() : the LSB_ECHKPNT_METHOD = blcr
                # Fri Feb 21 15:40:13 2014 : Echkpnt : main() : the LSB_ECHKPNT_METHOD_DIR =
                # Fri Feb 21 15:40:13 2014 : Echkpnt : main() : the echkpntProgPath is : /usr/local/lsf/9.1/linux2.6-glibc2.3-x86_64/etc/echkpnt.blcr
                # ########### Echkpnt end checkpoint ###########
                #
                # ########### begin to checkpoint ############
                # Fri Feb 21 16:23:12 2014 : Echkpnt : main() : the LSB_ECHKPNT_METHOD = blcr
                # Fri Feb 21 16:23:12 2014 : Echkpnt : main() : the LSB_ECHKPNT_METHOD_DIR =
                # Fri Feb 21 16:23:12 2014 : Echkpnt : main() : the echkpntProgPath is : /usr/local/lsf/9.1/linux2.6-glibc2.3-x86_64/etc/echkpnt.blcr
                # Fri Feb 21 16:23:12 2014 : Echkpnt : main() : the echkpnt.blcr fail,the exit value is 1
                # ########### Echkpnt end checkpoint ###########

		# wind back from last line looking for last "begin to checkpoint" and "end to checkpoint" lines
		my $last_begin_i = -1;
		my $last_end_i = -1;
		my $chkpnt_fail_count = 0;
		for ( my $i=@chkpnt_log_lines-1; $i>=0; $i-- )
		{
		    if ( $last_begin_i == -1 && $chkpnt_log_lines[$i] =~ m/^#+\s+begin.*checkpoint\s+#+\s*$/i ) 
		    {
			$last_begin_i = $i;
		    }
		    elsif ( $last_end_i == -1 && $chkpnt_log_lines[$i] =~ m/^#+\s+.*end.*checkpoint\s+#+\s*$/i ) 
		    {
			$last_end_i = $i;
		    }
		    if ( $chkpnt_log_lines[$i] =~ m/echkpnt[.]blcr fail/ )
		    {
			$chkpnt_fail_count++;
		    }
		}

		if ( $chkpnt_fail_count >= 3 )
		{
		    warn(
			"\tc  job $$job{id} ($$job{lsf_id}) has failed to checkpoint $chkpnt_fail_count times.\n" .
			"\tc  killing it so that it can be restarted.\n"
			);
		    kill_job($job);
		    return;
		}
		
		if ( $last_begin_i >= 0 )
		{
		    # have a started checkpoint, get start time
		    my $timestamp_line = $chkpnt_log_lines[$last_begin_i+1];
		    my $chkpnt_start_time;
		    # Fri Feb 21 16:23:12 2014 : Echkpnt : main() : the LSB_ECHKPNT_METHOD = blcr
		    if ( $timestamp_line =~ m/^\w+\s+(\w+)\s+(\d+) (\d+):(\d+):(\d+)\s+(\d+)\s*:/ ) 
		    {
			$chkpnt_start_time = DateTime->new(month=>$months{$1}, day=>$2, hour=>$3, minute=>$4, second=>$5, year=>$6)->epoch;
		    } 
		    else 
		    {
			confess("Have started checkpoint for job $$job{id} ($$job{lsf_id}), but could not find timestamp [$timestamp_line] in checkpoint log file ($chkpnt_log_file)\n");
		    }
		    
		    if ( $last_end_i < $last_begin_i ) # N.B. this relies on initial valie of last_end_i being -1
		    {
			# checkpoint in progress, check how long it has been in progress
			my $chkpnt_elapsed_seconds = $$job{cur_time} - $chkpnt_start_time;
			my $max_chkpnt_time_seconds = $max_chkpnt_time_s_per_mb * $$job{mem_mb};
			if ( $chkpnt_elapsed_seconds > $max_chkpnt_time_seconds )
			{
			    warn(
				"\tc  job $$job{id} ($$job{lsf_id}) has been checkpointing for a long time (${chkpnt_elapsed_seconds}s) and will likely never finish.\n" .
				"\tc  killing it so that it can be restarted.\n"
				);
			    kill_job($job);
			}
			else 
			{
			    # allow it to finish (do nothing for now)
			    warn(
				"\tc  job $$job{id} ($$job{lsf_id}) is currently checkpointing (it has been ${chkpnt_elapsed_seconds}s so far), waiting for it to finish...\n"
				);
			}
		    }
		    elsif ( $last_end_i > $last_begin_i ) 
		    {
			my $chkpnt_context_file = $$job{chkpnt_dir}."/jobstate.context";
			my $chkpnt_context_mtime = undef;
			if ( -e $chkpnt_context_file && -s $chkpnt_context_file && (my @st = stat($chkpnt_context_file)) )
			{
			    my ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size,$atime,$mtime,$ctime,$blksize,$blocks) = @st;
			    $chkpnt_context_mtime = $mtime;
			}
			else
			{
			    warn("\tc  could not stat checkpoint context file $chkpnt_context_file, checkpoint probably failed\n");
			}
			
			if ( defined($chkpnt_context_mtime) && ( $chkpnt_context_mtime >= ( $chkpnt_start_time - 30 ) ) ) # allow 30 second drift between timestamp and log time
			{
			    # the latest checkpoint was successful, kill job so it can be restarted
			    warn(
				"\tc  job $$job{id} ($$job{lsf_id}) is over queue time limit but has a recent checkpoint.\n" .
				"\tc  killing it so that it can be restarted.\n"
				);
			    kill_job($job);
			}
			else 
			{
			    # the latest checkpoint did not result in an updated jobstate.context
			    warn(
				"\tc  job $$job{id} ($$job{lsf_id}) appears to have failed its latest checkpoint, attempting checkpoint-and-kill now.\n"
				);
			    
			    # ask it to checkpoint-and-kill
			    my $bchkpnt_cmd = "bchkpnt -k '$$job{lsf_id}'";
			    my @out = `$bchkpnt_cmd`;
			    # Job <1878095[1]> is being checkpointed
			    if ( scalar @out!=1 || !($out[0]=~/^Job <([\d\[\]]+)> is being checkpointed/) )
			    {
				warn(
				    "Expected different output from bchkpnt.\n" .
				    "The bchkpnt command was:\n" .
				    "\t$bchkpnt_cmd\n" .
				    "The output was: " . join('',@out) . "\n"
				    );
			    }
			}
		    }
		}
		else 
		{
		    # no checkpoints yet -- this is unexpected (perhaps the job is running with more cpus than expected?)
		    warn( 
			"\tc  job $$job{id} ($$job{lsf_id}) has run past the queue limit but has not yet attempted to checkpoint.\n" .
			"\tc  This could be because you are using more CPUs than expected.\n" . 
			"\tc  This job may soon be killed by LSF.\n" .
			"\tc  You could try to save it by initiating a checkpoint manually using the command: bchkpnt '$$job{lsf_id}'\n" 
			);
		}
	    }
	    else # !exists($$job{chkpnt_dir})
	    {
		# use bswitch instead of checkpoint/restart, for some reason only after 1.3x the cpu time limit has gone
		if ( exists($queue_limits_seconds{$queue}) && ( $$job{cpu_time}*1.3 > $queue_limits_seconds{$queue} ) )
		{
		    my $bswitch;
		  QUEUE: for my $q (sort {$queue_limits_seconds{$a} <=> $queue_limits_seconds{$b}} keys %queue_limits_seconds)
		  {
		      if ( $$job{cpu_time}*1.3 > $queue_limits_seconds{$q} ) { next; }
		      warn("\ts changing queue of the job $$job{lsf_id} from $$job{queue} to $q\n");
		      my $cmd = "bswitch $q '$$job{lsf_id}'";
		      print STDERR "calling $cmd\n";
		      my $out = `$cmd`;
		      print STDERR "output: [$out]\n";
		      $$job{queue} = $q;
		      last QUEUE;
		  }
		}
	    } # not using checkpointing
        } # over time
    } # running
    return;
}

sub kill_job
{
    my ($job) = @_;

    my $bkill_cmd = "bkill -s KILL '$$job{lsf_id}'";
    my @out = `$bkill_cmd`;
    # Job <1878080> is being terminated
    # Job <1878081[1]> is being terminated
    if ( scalar @out!=1 || !($out[0]=~/^Job <([\d\[\]]+)> is being terminated/) )
    {
	warn(
	    "Expected different output from bkill.\n" .
	    "The bkill command was:\n" .
	    "\t$bkill_cmd\n" .
	    "The output was: " . join('',@out) . "\n"
	    );
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
            if ( $1 eq 'Exited' ) { $$out{status} = $Error; $$out{nfailures}++; }
            if ( $1 eq 'Done' ) { $$out{status} = $Done; $$out{nfailures} = 0; }
        }
	# TERM_OWNER: job killed by owner.
	if ( $line =~ /^TERM_OWNER:/ && $$out{nfailures} > 0 )
	{
	    $$out{nfailures}--;
	}
	# TERM_RUNLIMIT: job killed after reaching LSF run time limit.
	if ( $line =~ /^TERM_RUNLIMIT:/ && $$out{nfailures} > 0 )
	{
	    $$out{nfailures}--;
	}
	# TERM_CHKPNT
	if ( $line =~ /^TERM_CHKPNT/ && $$out{nfailures} > 0 )
	{
	    $$out{nfailures}--;
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
    my $mem_mb;
    if ( ! exists($$opts{memory}) )
    {
	$$opts{memory} = $default_memlimit_mb;
    }
    if ( $$opts{memory} ) 
    { 
        $mem_mb = int($$opts{memory});
        if ( $mem_mb )
        {
            my $units = get_lsf_limits_unit();
            my $lmem  = $units eq 'kB' ? $mem_mb*1000 : $mem_mb;
            $bsub_opts = sprintf " -M%d -R 'select[type==X86_64 && mem>%d] rusage[mem=%d]'", $lmem,$mem_mb,$mem_mb; 
        }
    }
    if ( !defined($$opts{runtime}) && defined($$opts{runtime_limit_seconds}) )
    {
	$$opts{runtime} = $$opts{runtime_limit_seconds};
    }
    if ( !defined($$opts{queue}) ) 
    {            
        if ( !$$opts{_run_with_blcr} && defined($$opts{runtime}) ) 
        { 
            if ( $$opts{runtime} <= $queue_limits_seconds{normal}/60 ) { $$opts{queue} = 'normal'; }
            elsif ( $$opts{runtime} <= $queue_limits_seconds{long}/60 ) { $$opts{queue} = 'long'; }
            else { $$opts{queue} = 'basement'; }
        }
        else 
        { 
	    # default to normal queue, and always use normal with checkpointing unless queue is explicitly set
            $$opts{queue} = 'normal';
        }
    }
    confess "queue not defined (this should not occur)\n" unless ( defined($$opts{queue}) );
    $bsub_opts .= " -q $$opts{queue}";

    if ( defined($$opts{cpus}) ) 
    {
        $bsub_opts .= " -n $$opts{cpus} -R 'span[hosts=1]'";
    }
    if ( $$opts{_run_with_blcr} ) 
    {
	my $runtime_limit_seconds = $queue_limits_seconds{$$opts{queue}};
	if ( defined($$opts{runtime_limit_seconds}) ) { $runtime_limit_seconds = $$opts{runtime_limit_seconds}; }
	if ( !defined($$opts{chkpnts_per_run}) ) 
	{
	    $$opts{chkpnts_per_run} = 1;
	}
	if ( !defined($$opts{chkpnt_period_minutes}) )
	{
	    my $chkpnts_per_run = int($$opts{chkpnts_per_run});
	    do 
	    {
		$$opts{chkpnt_period_minutes} = ((($runtime_limit_seconds - ($min_chkpnt_time_s_per_mb * $mem_mb)) / $chkpnts_per_run) / 60);
		# reduce checkpoints per run and repeat if calculated checkpoint period is below the minimum
		$chkpnts_per_run--;
	    } while( ($$opts{chkpnt_period_minutes} < $chkpnt_minimum_period_minutes) && ($chkpnts_per_run > 0) );
	    $$opts{chkpnts_per_run} = $chkpnts_per_run + 1;
	    # if checkpoint period is still below minimum, just set to minimum
	    if ( $$opts{chkpnt_period_minutes} < $chkpnt_minimum_period_minutes )
	    {
		warn(
		    "Checkpoint period for job ($$opts{id}) of $$opts{chkpnt_period_minutes}m is less than the minimum of ${chkpnt_minimum_period_minutes}m.\n" .
		    "Using the minimum.\n"
		    );
		$$opts{chkpnt_period_minutes} = $chkpnt_minimum_period_minutes;
	    }
	}
	my $cwd = `pwd`;
	chomp $cwd;
	$jids_file =~ /(.*\/\.jobs\/).*jid/;
	my $temp_dir = $1;
	my $chkpnt_dir = $cwd."/".$temp_dir."chkpnt";
        if ( !$$opts{chkpnt_period_minutes} ) 
	{ 
	    $$opts{chkpnt_period_minutes} = 350; 
	}
	else 
	{
	    $$opts{chkpnt_period_minutes} = int($$opts{chkpnt_period_minutes}); 
	}
	my $runtime_limit_minutes = int($runtime_limit_seconds/60);
        $bsub_opts .= " -k '$chkpnt_dir method=blcr $$opts{chkpnt_period_minutes}' -We $runtime_limit_minutes";
	if ( ! ( $cmd =~ m/^\s*cr[_]/ )  )
	{ 
	    $cmd = "cr_run $cmd"; 
	}
    }

    my $bsub_cmd  = qq[bsub -J '${job_name}[$bsub_ids]' -e $job_name.\%I.e -o $job_name.\%I.o $bsub_opts '$cmd'];

    # Submit to LSF
    print STDERR "$bsub_cmd\n";
    my @out = `$bsub_cmd 2>&1`;
    if ( scalar @out!=1 || !($out[0]=~/^Job <(\d+)> is submitted/) )
    {
	my $cwd = `pwd`;
	confess(
	    "Expected different output from bsub.\n" .
	    "The command was:\n" .
	    "\t$cmd\n" .
	    "The bsub_command was:\n" .
	    "\t$bsub_cmd\n" .
	    "The working directory was:\n" .
	    "\t$cwd\n" .
	    "The output was: " . join('',@out) . "\n"
	    );
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
    my ($jids_file, $job, $opts, $resources_changed) = @_;
    confess("RunnerLSF::restart_job called on job without id entry") unless exists($$job{id});

    my $id = $$job{id};
    my $restarted = 0;
### Commenting out brestart option as it seems to frequently produce un-checkpointable jobs after restart
    # if ( $resources_changed )
    # {
    # 	warn("\tc  resource requirements have changed for job ($id)\n");
    # }
    # else 
    # {
    # 	# resource requirements have not changed, try using brestart
    # 	warn("\tc  using brestart to restart job ($id) as resource requirements have not changed.\n");
    # 	my $brestart_opts = "";
	
    # 	if ( defined($$job{queue}) )
    # 	{
    # 	    $brestart_opts .= "-q $$job{queue} ";
    # 	}
	
    # 	if ( defined($$job{mem_mb}) )
    # 	{
    # 	    my $units = get_lsf_limits_unit();
    # 	    my $lmem = $units eq 'kB' ? ($$job{mem_mb}*1000) : $$job{mem_mb};
    # 	    $brestart_opts .= "-M $lmem ";
    # 	}
    # 	my @chkpnt_lsf_dir_comps = File::Spec->splitdir($$job{chkpnt_dir});
    # 	pop @chkpnt_lsf_dir_comps;
    # 	my $chkpnt_lsf_dir = File::Spec->catdir(@chkpnt_lsf_dir_comps);
    # 	my $brestart_cmd = qq[brestart $brestart_opts $chkpnt_lsf_dir '$$job{last_chkpnt_lsf_id}'];
	
    # 	if ( !defined('$$job{name}') )
    # 	{
    # 	    confess("Job name not defined for $id ($$job{last_chkpnt_lsf_id})");
    # 	}
	
    # 	# Submit to LSF
    # 	my @out = `$brestart_cmd 2>&1`;
    # 	if ( scalar @out!=1 || !($out[0]=~/^Job <(\d+)> is submitted/) )
    # 	{
    # 	    warn(
    # 		"Expected different output from brestart.\n" .
    # 		"The brestart_command was:\n" .
    # 		"\t$brestart_cmd\n" .
    # 		"The output was: " . join('',@out) . "\n"
    # 		);
    # 	    # Example: Checkpoint log is not found or is corrupted. Job not submitted.
    # 	}
    # 	else 
    # 	{
    # 	    # Write down info about the submitted command
    # 	    my $jid = $1;
    # 	    print STDERR "Job name set to $$job{name}\n";
    # 	    open(my $jids_fh, '>>', $jids_file) or confess("$jids_file: $!");
    # 	    print $jids_fh join("\t", $jid, $$job{name} , $brestart_cmd)."\n";
    # 	    close $jids_fh;
	    
    # 	    $restarted = 1;
    # 	}
    # }
    
    if ( ! $restarted )
    {
	# brestart failed or resource requirements have changed, submit using bsub and cr_restart
	my @id;
	push @id, $id;
	my $chkpnt_file = $$job{chkpnt_dir}."/jobstate.context";
	my $cmd = "cr_restart -f $chkpnt_file";
	my $job_name = $$job{name};
	warn("\tc  attempting to restart job $id using cr_restart.\n");
        run_array($jids_file, $job_name, $opts, $cmd, \@id);
    }
}

=head1 AUTHORS

petr.danecek@sanger
Allan Daly <ad7@sanger.ac.uk>
Joshua Randall <jr17@sanger.ac.uk>

=cut
1;

