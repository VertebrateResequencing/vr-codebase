#
# Copyright (c) 2013, 2014 Genome Research Ltd. 
# 
# Authors: Allan Daly <ad7@sanger.ac.uk>
#          Joshua Randall <jr17@sanger.ac.uk>
#          Petr Danecek <pd3@sanger> 
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

# Checkpointing does not seem to be very robust, this is a list of unsolved
# errors one might encounter:
#   - 
#   - The job's error output file contains the error message:
#       Failed to open path/to/.jobs/name.w.chkp/[0-9]+/jobstate.context: No such file or directory 
#   - The checkpoint output file path/to/.jobs/name.w.chkp/[0-9]+/echkpnt.err contains the error:
#       Signal 18 (CONT) caught by ps (procps version 3.2.8).
#

package RunnerLSFCR;
use base qw(RunnerLSF);

use strict;
use warnings;
use Carp;
use Cwd;

sub new
{
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    # Estimated (minimum) time to checkpoint: ~2.5 minutes per GB of RAM on lustre (scratch113), so 0.15 s/MB 
    # But it is better to underestimate than overestimate - theoretical maximum is ~500MB/s, so 0.002s/MB
    $$self{default_limits}{chkpnt_MBps} = 0.002;

    return $self;
}

sub clean_jobs
{
    my ($self,$wfile,$ids,$all_done) = @_;
    for my $id (@$ids)
    {
        $self->SUPER::clean_jobs($wfile,$id);
        my @chks = $self->_read_checkpoints("$wfile.$id.chkp");
        for my $chk (@chks) { `rm -rf $$chk{chkdir}`; }
        unlink("$wfile.$id.chkp");
    }
    if ( $all_done )
    {
        `rm -rf $wfile.chkp`;
    }
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

    if ( $$job{status}!=$$self{Running} ) { return; }

    if ( !exists($$job{cpu_time}) ) { $$job{cpu_time} = 0; }
    my $wakeup_time = $$self{limits}{wakeup_interval} ? $$self{limits}{wakeup_interval} + 900 : 900;    # reserve at least 15 mins for the checkpoint
    my $chkpnt_time = $$job{mem_usage} * ($$self{limits}{chkpnt_MBps} ? $$self{limits}{chkpnt_MBps} : $$self{default_limits}{chkpnt_MBps});
    my $cur_queue   = $$job{queue};

    # Estimate how long it will take before we are called again. If we are near the time limit, create a checkpoint
    #if ( $$self{queue_limits}{$cur_queue}*60 > $$job{cpu_time} + $wakeup_time + $chkpnt_time ) { return; }

    if ( $$self{queue_limits}{$cur_queue}*60 > $$job{cpu_time} + $wakeup_time + $chkpnt_time ) { return; }

    my $chk_file = $jids_file;
    $chk_file =~ s/\.jid$//;
    $chk_file = "$chk_file.$$job{id}.chkp";
    my @chks = $self->_read_checkpoints($chk_file);

    if ( !@chks or $chks[-1]{chkdir} ne $$job{chkpnt_dir} ) 
    { 
        push @chks, { stamp=>time(), status=>'init', chkdir=>$$job{chkpnt_dir} };
    }
    if ( $chks[-1]{status} eq 'wait' )
    {
        # Waiting until bchkpnt finishes. Since the job is still running,
        # either it is taking too long to checkpoint or it failed.
        my $ret = $self->_check_bchkpnt_status($chks[-1]{chkdir});
        if ( $ret<0 ) { return; }   # chkpointing has not finished yet
        if ( $ret==0  )
        {
            # It ran OK and mark it as ready-to-use.
            $chks[-1]{status} = 'ready';
            $self->_write_checkpoints($chk_file,\@chks);
            return;
        }
        if ( $ret>=2 )
        {
            # Failed twice, fallback to bswitch
            $chks[-1]{status} = 'broken';
            $self->_write_checkpoints($chk_file,\@chks);
            $self->SUPER::_check_job($job);
            return;
        }
        $chks[-1]{status} = 'init';  # Failed, try again
    }
    if ( $chks[-1]{status} eq 'init' )
    {
        $chks[-1]{status} = 'wait';
        $self->_write_checkpoints($chk_file,\@chks);

        warn("Checkpointing job $$job{lsf_id}: $chks[-1]{chkdir}\n");
        sleep(1);   # prevent two many jobs checkpointed simultaneously

        # Errors:
        #   chkpnt.log
        #       Echkpnt : main() : the echkpntProgPath is : /usr/local/lsf/9.1/linux2.6-glibc2.3-x86_64/etc/echkpnt.blcr
        #       Echkpnt : main() : the echkpnt.blcr fail,the exit value is 255
        #   echkpnt.err
        #       Signal 18 (CONT) caught by ps (procps version 3.2.8).
        #       Please send bug reports to <feedback@lists.sf.net> or <albert@users.sf.net>
        #       Usage: cr_checkpoint [options] ID
        #
        my $cmd = "bchkpnt -k '$$job{lsf_id}'";
        my @out = `$cmd 2>&1`; 
        if ( scalar @out!=1 || !($out[0]=~/^Job <.+> is being checkpointed/) )
        {
            # Job <1526270[47]>: Job has not started yet
            if ( $out[0] =~ /Job has not started yet/ ) { return; }

            # Job <1878095[1]> is being checkpointed
            warn("Expected different output from bchkpnt. The bchkpnt command was:\n\t$cmd\nThe output was:\n". join('',@out));
        }
        if ( $? ) { $self->SUPER::_check_job($job); }    # fallback to bswitch
        return;
    }
}

# Returns 0 on last successful checkpoint, negative value when chkp not
# finished or the number of unsuccessful attempts.
sub _check_bchkpnt_status
{
    my ($class,$chkdir) = @_;
    my $file = "$chkdir/chkpnt.log";
    open(my $fh,'<',$file) or return -1;
    my $nbeg  = 0;
    my $nend  = 0;
    my $nfail = 0;
    my $fail  = 0;
    while (my $line=<$fh>)
    {
        # ########### begin to checkpoint ############
        if ( $line=~/^#+\s+begin.*checkpoint\s+#/ ) 
        { 
            $fail = 0;
            $nbeg++;
            next; 
        }
        # ########### Echkpnt end checkpoint ###########
        if ( $line=~/^#+\s+.+end.*checkpoint\s+#/ ) 
        { 
            if ( $fail ) { $nfail++; }
            $nend++; 
            next;
        }
        if ( $line=~/fail/ ) { $fail = 1; }
    }
    close($fh) or confess("close failed: $file");
    if ( !$nbeg or $nbeg!=$nend ) { return -1; }    # chkpnt not finished yet
    if ( !$fail ) { return 0; } # last chkpnt succeeded
    return $nfail;   # number of failures
}
sub _read_checkpoints
{
    my ($class,$file) = @_;
    my @out = ();
    open(my $fh,'<',$file) or return (@out);
    while (my $line=<$fh>)
    {
        chomp($line);
        my ($stamp,$status,$chkdir) = split(/\t/,$line);
        push @out, { stamp=>$stamp, status=>$status, chkdir=>$chkdir };
    }
    close($fh) or confess("close failed: $file");
    return (@out);
}
sub _write_checkpoints
{
    my ($class,$file,$chks) = @_;
    open(my $fh,'>',$file) or confess("$file: $!");
    for my $chk (@$chks)
    {
        print $fh "$$chk{stamp}\t$$chk{status}\t$$chk{chkdir}\n";
    }
    close($fh) or confess("close failed: $file");
}

sub _checkpoint_ready
{
    my ($self,$job_name,$id) = @_;

    my @chks = $self->_read_checkpoints("$job_name.$id.chkp");
    if ( @chks && $chks[-1]{status} eq 'wait' )
    {
        if ( !$self->_check_bchkpnt_status($chks[-1]{chkdir}) )
        {
            $chks[-1]{status} = 'ready';
            $self->_write_checkpoints("$job_name.$id.chkp", \@chks);
        }
    }
    if ( !@chks ) { return undef; }
    return $chks[-1]{status} eq 'ready' ? "$chks[-1]{chkdir}/jobstate.context" : undef;
}

sub run_jobs
{
    my ($self,$jids_file,$job_name,$cmd,$ids) = @_;
    if ( !scalar @$ids ) { confess("No IDs given??\n"); }

    my $cwd = getcwd();

    $cmd =~ s/{JOB_INDEX}/\$LSB_JOBINDEX/g;
    my $bsub_opts = $self->_create_bsub_opts_string();

    my (@run_ids,@restart_ids,%chks);
    for my $id (@$ids)
    {
        my $chk = $self->_checkpoint_ready($job_name,$id);
        if ( $chk ) { push @restart_ids,$id; $chks{$id} = $chk;  }
        else { push @run_ids,$id; }
    }

    my $cr_opts  = qq[-k '$cwd/$job_name.chkp method=blcr'];
    while ( @run_ids )
    {
        # NB: if the command stdout is connected to the same place where LSF
        #   writes the job info, the LSF info will be truncated on restart,
        #   redirect stdout to /dev/null
        my $bsub_ids = $self->_create_bsub_ids_string($job_name,\@run_ids);
        my $bsub_cmd = qq[bsub -rn -J '${job_name}[$bsub_ids]' -e $job_name.\%I.e -o $job_name.\%I.o $bsub_opts $cr_opts 'cr_run $cmd > /dev/null'];
        $self->_bsub_command($jids_file,$job_name,$bsub_cmd);
    }
    for my $id ( @restart_ids )
    {
        # Cannot submit multiple jobs at once via job array, cr_restart cannot be run from a script to stay checkpointable via bchkpnt.
        # The --run-on-fail-temp makes sure that on temporary errors (such as PIDs being in use) will not be reported as failures.
        my $bsub_cmd = qq[bsub -rn -J '${job_name}[$id]' -e $job_name.\%I.e -o $job_name.\%I.o $bsub_opts $cr_opts 'cr_restart --run-on-fail-temp=/bin/true $chks{$id} > /dev/null'];
        $self->_bsub_command($jids_file,$job_name,$bsub_cmd);
    }
}


=head1 AUTHORS

Allan Daly <ad7@sanger.ac.uk>
Joshua Randall <jr17@sanger.ac.uk>
petr.danecek@sanger

=cut

1;

