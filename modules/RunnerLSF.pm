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

#my %queue_limits_seconds = ( basement=>1e9, long=>2880*60, normal=>720*60, small=>30*60 );
#my $chkpnt_minimum_period_minutes = 30;
my %queue_limits_seconds = ( basement=>1*60*60, long=>1*60*60, normal=>1*60*60, small=>30*60 );
my $chkpnt_minimum_period_minutes = 20;
my $default_memlimit_mb = 100;

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
	    if ( exists($$info{$id}{cpus}) )
	    {
		$jobs[$j]{cpus} = $$info{$id}{cpus};
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
      # roll up all continuation lines, signified by beginning with 21 spaces 
      # into a single (long) line beginning with the line before
      $output =~ s/\r//gm;
      $output =~ s/\n\s{21}//gm;
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

        # Tue Mar 19 13:58:23: Resource usage collected...
        if ( $line =~ /^\w+\s+(\w+)\s+(\d+) (\d+):(\d+):(\d+):\s+Resource usage collected/ ) 
        {
            if ( !exists($$job{started}) ) { confess("No wall time for job $$job{id}??", $line); }
            my $wall_time = DateTime->new(month=>$months{$1}, day=>$2, hour=>$3, minute=>$4, year=>$year)->epoch - $$job{started};
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
		    "Job ($$job{lsf_id}) appears to have died while checkpointing"
		    );
	    }
	    elsif ( $line =~ /exit code (\d+)/ ) 
	    {
		warn(
		    "Job ($$job{lsf_id}) had unexpected error status while checkpointing (exit code $1)\n"
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
        if ( exists($queue_limits_seconds{$queue}) && ( ( ($$job{cpu_time} / $cpus) > $queue_limits_seconds{$queue} ) || ( $$job{wall_time} > $queue_limits_seconds{$queue} ) ) )
        {
	    warn(
		"Job $$job{id} ($$job{lsf_id}) is currently running and has exceeded the $queue queue limit of $queue_limits_seconds{$queue} seconds (cpu time: $$job{cpu_time} wall time: $$job{wall_time})\n"
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
		}

		if ( $last_end_i < $last_begin_i ) # N.B. this relies on initial valie of last_end_i being -1
		{
		    # checkpoint in progress, allow it to finish (do nothing for now)
		    warn(
			"Running job $$job{id} ($$job{lsf_id}) is currently checkpointing, waiting for it to finish...\n"
			);
		}
		elsif ( $last_end_i > $last_begin_i ) 
		{
		    # have a completed checkpoint, see if it was successful
		    my $timestamp_line = $chkpnt_log_lines[$last_begin_i+1];
		    my $chkpnt_start_time;
		    # Fri Feb 21 16:23:12 2014 : Echkpnt : main() : the LSB_ECHKPNT_METHOD = blcr
		    if ( $timestamp_line =~ m/^\w+\s+(\w+)\s+(\d+) (\d+):(\d+):(\d+)\s+(\d+)\s*:/ ) 
		    {
			$chkpnt_start_time = DateTime->new(month=>$months{$1}, day=>$2, hour=>$3, minute=>$4, second=>$5, year=>$6)->epoch;
		    } 
		    else 
		    {
			confess("Have completed checkpoint for job $$job{id} ($$job{lsf_id}), but could not find timestamp [$timestamp_line]\n");
		    }
		    my $chkpnt_context_file = $$job{chkpnt_dir}."/jobstate.context";
		    my $chkpnt_context_mtime = undef;
		    if ( -e $chkpnt_context_file && -s $chkpnt_context_file && (my @st = stat($chkpnt_context_file)) )
		    {
			my ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size,$atime,$mtime,$ctime,$blksize,$blocks) = @st;
			$chkpnt_context_mtime = $mtime;
		    }
		    else
		    {
			warn("Could not stat checkpoint context file $chkpnt_context_file, checkpoint probably failed\n");
		    }
		    if ( defined($chkpnt_context_mtime) && ( $chkpnt_context_mtime >= ( $chkpnt_start_time - 30 ) ) ) # allow 30 second drift between timestamp and log time
		    {
			# the latest checkpoint was successful, kill job so it can be restarted
			warn(
			    "Running job $$job{lsf_id} is over queue time limit but has a recent checkpoint.\n" .
			    "Killing it so that it can be restarted.\n"
			    );
			my $cmd = "bkill -s KILL '$$job{lsf_id}'";
			print STDERR "Calling $cmd\n";
			my $out = `$cmd`;
			print STDERR "bkill output: $out\n";
		    }
		    else 
		    {
			# the latest checkpoint did not result in an updated jobstate.context
			warn(
			    "Job $$job{lsf_id} appears to have failed its latest checkpoint\n"
			    );
			
			# ask it to checkpoint-and-kill
			my $cmd = "bchkpnt -k '$$job{lsf_id}'";
			print STDERR "Calling $cmd\n";
			my $out = `$cmd`;
			print STDERR "bchkpnt output: $out\n";
		    }
		}
		else 
		{
		    # no checkpoints yet -- this is unexpected (perhaps the job is running with more cpus than expected?)
		    confess( 
			"\nA running job ($$job{lsf_id}) has run past the queue limit but has not yet attempted to checkpoint.\n" .
			"(according to $chkpnt_log_file).\n" .
			"You could try initiating a checkpoint manually using bchkpnt: bchkpnt '$$job{lsf_id}'\n\n" 
			);
		}
	    }
	    else # !exists($$job{chkpnt_dir})
	    {
		# use bswitch instead of checkpoint/restart, for some reason only after 1.3x the cpu time limit has gone
		if ( $$job{cpu_time}*1.3 > $queue_limits_seconds{$queue} )
		{
		    my $bswitch;
		    QUEUE: for my $q (sort {$queue_limits_seconds{$a} <=> $queue_limits_seconds{$b}} keys %queue_limits_seconds)
		    {
			if ( $$job{cpu_time}*1.3 > $queue_limits_seconds{$q} ) { next; }
			warn("Changing queue of the job $$job{lsf_id} from $$job{queue} to $q\n");
			my $cmd = "bswitch $q '$$job{lsf_id}'";
			print STDERR "calling $cmd\n";
			my $out = `$cmd`;
			print STDERR "output: [$out]\n";
			$$job{queue} = $q;
			last QUEUE;
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

            if ( !exists($out{memory_mb}) or $out{memory_mb}<$mem )
            {
                $out{memory_mb} = $mem;
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
    if ( ! exists($$opts{memory_mb}) )
    {
	$$opts{memory_mb} = $default_memlimit_mb;
    }
    if ( $$opts{memory_mb} ) 
    { 
        $mem_mb = int($$opts{memory_mb});
        if ( $mem_mb )
        {
            my $units = get_lsf_limits_unit();
            my $lmem  = $units eq 'kB' ? $mem_mb*1000 : $mem_mb;
            $bsub_opts = sprintf " -M%d -R 'select[type==X86_64 && mem>%d] rusage[mem=%d]'", $lmem,$mem_mb,$mem_mb; 
        }
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
	if ( !defined($$opts{chkpnts_per_run}) ) 
	{
	    $$opts{chkpnts_per_run} = 1;
	}
	if ( !defined($$opts{chkpnt_time_s_per_mb}) )
	{
            # estimated time to checkpoint: ~9 minutes per GB of RAM on lustre (scratch113), so 0.54 s/MB
	    $$opts{chkpnt_time_s_per_mb} = 0.54;
	}
	if ( !defined($$opts{chkpnt_period_minutes}) )
	{
	    my $chkpnts_per_run = int($$opts{chkpnts_per_run});
	    do 
	    {
		$$opts{chkpnt_period_minutes} = ($queue_limits_seconds{$$opts{queue}} - ($$opts{chkpnt_time_s_per_mb} * $mem_mb)) / $chkpnts_per_run;
		# reduce checkpoints per run and repeat if calculated checkpoint period is below the minimum
		$chkpnts_per_run--;
	    } while( ($$opts{chkpnt_period_minutes} < $chkpnt_minimum_period_minutes) && ($chkpnts_per_run > 0) );
	    $$opts{chkpnts_per_run} = $chkpnts_per_run + 1;
	    # if checkpoint period is still below minimum, just set to minimum
	    if ( $$opts{chkpnt_period_minutes} < $chkpnt_minimum_period_minutes )
	    {
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
        $bsub_opts .= " -k '$chkpnt_dir method=blcr $$opts{chkpnt_period_minutes}'";
	if ( ! $cmd =~ m/\s*cr_/ ) { $cmd = "cr_run $cmd"; }
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
    my ($jids_file, $job, $opts, $resources_changed) = @_;
    confess("RunnerLSF::restart_job called on job without id entry") unless exists($$job{id});

    my $id = $$job{id};
    if ( $resources_changed > 0 )
    {
	# resource requirements have changed, submit using bsub and cr_restart
	my @id;
	push @id, $id;
	my $chkpnt_file = $$job{chkpnt_dir}."/jobstate.context";
	my $cmd = "cr_restart -f $chkpnt_file";
	my $job_name = $$job{name};
	print STDERR "Resource requirements have changed, calling run_array to restart job ($id) using cr_restart.\n";
        run_array($jids_file, $job_name, $opts, $cmd, \@id);
    }
    else
    {
	# resource requirements have not changed, use brestart
	print STDERR "Using brestart to restart job ($id) as resource requirements have not changed.\n";
	my $brestart_opts = "";
	
	if ( defined($$job{queue}) )
	{
	    $brestart_opts .= "-q $$job{queue} ";
	}
	
	if ( defined($$job{mem_mb}) )
	{
	    my $units = get_lsf_limits_unit();
	    my $lmem = $units eq 'kB' ? ($$job{mem_mb}*1000) : $$job{mem_mb};
	    $brestart_opts .= "-M $lmem ";
	}
	
	my $brestart_cmd = qq[brestart $brestart_opts  '$$job{chkpnt_dir}' '$$job{last_chkpnt_lsf_id}'];
	
	if ( !defined('$$job{name}') )
	{
	    confess("Job name not defined for $id ($$job{last_chkpnt_lsf_id})");
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
	print STDERR "Job name set to $$job{name}\n";
	open(my $jids_fh, '>>', $jids_file) or confess("$jids_file: $!");
	print $jids_fh join("\t", $jid, $$job{name} , $brestart_cmd)."\n";
	close $jids_fh;
    }
}

=head1 AUTHORS

petr.danecek@sanger
Allan Daly <ad7@sanger.ac.uk>
Joshua Randall <jr17@sanger.ac.uk>

=cut
1;

