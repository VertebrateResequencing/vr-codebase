package VertRes::Pipeline;

use strict;
use warnings;
use Carp;   # this is for confess in lock_file and unlock_file. eventually should go somewhere else

use base qw(VertRes::Base);
use LSF;
use Fcntl qw(:DEFAULT :flock);

our $Yes     = 0;
our $No      = 1;
our $Error   = 2;
our $Running = 3;

our $default_options =
{
    # General options, common to all modules
    'check_status'   => 0,
    'clean'          => 0,
    'exit_on_errors' => 1,
    'global_lock'    => 'Pipeline.lock',
    'lane_path'      => '',
    'list_tasks'     => 0,
    'log'            => '_log',
    'prefix'         => '_',
    'rerun'          => 0,
    'stepwise'       => 1,
    'task'           => '',

    'Error'          => $Error,
    'Yes'            => $Yes,
    'No'             => $No,
    'Running'        => $Running,
};


=head2 new

        Options    :
                    check_status    .. if set to 1, nothing will be run, only the status of the lane will be reported
                    clean           .. delete all files provided by actions
                    exit_on_errors  .. if set to 0, an unfinished task will be run again even if the LSF job previously failed
                    global_lock     .. lock file for the lane, to prevent multiple run_lane subroutines running
                    lane_path       .. optional, can be overrided by run_lane
                    list_tasks      .. don't do anything, only list available tasks
                    log             .. this is used when logfile is not given (may be both relative and absolute path)
                    prefix          .. to distinguish between multiple variants of the same module, but running with different options
                    rerun           .. if set to 1, always run again. Implies that the task option is set.
                    stepwise        .. execute only one task at time
                    task            .. if set, execute only the given tasks (separated by commas, the order is ignored)
        Example    : See TrackQC.pm for an example.

=cut

sub new 
{
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(%$default_options, @args);
    return $self;
}



=head2 run_lane

        Arg [1]    : the lane to be processed (absolute path). Can be ommited when supplied with new().
        Example    : See TrackQC.pm for an example.
        Returntype : $Yes if the lane has finished, $Running when the lane is running or queued to LSF,
                        $No when unfinished and not running, $Error if an error occured and not running.

=cut

sub run_lane
{
    my ($self,$lane_path) = @_;

    $lane_path = $$self{'lane_path'} unless $lane_path;
    if ( !$lane_path ) { $self->throw("Missing parameter: the lane to be run.\n") }
    if ( !-d $lane_path ) { $self->throw("The directory does not exist:\n\t$lane_path\n") }
    if ( !($lane_path =~ m{^/}) ) { $self->throw("The lane path must be an absolute path.\n") }

    if ( $$self{'rerun'} && !$$self{'task'} ) { $self->throw("The option \"rerun\" cannot be used without \"task\".\n") }
    if ( $$self{'task'} )
    {
        my @tasks = split(/\s*,\s*/,$$self{'task'});
        $$self{'task'} = {};
        for my $task (@tasks) { $$self{'task'}{$task}=1; }
    }

    $$self{'logfile'} = ($$self{'log'} =~ m{^/}) ? $$self{'log'} : "$lane_path/$$self{'log'}";
    $self->log_file($$self{'logfile'});

    my $read_only = ($$self{'check_status'} || $$self{'list_tasks'}) ? 1 : 0;

    # The cleaning should be probably called directly, not via run_lane...? 
    #   If some of the read-only options is set, choose to do the read only things.
    if ( $$self{clean} && !$read_only )
    {
        $self->clean();
        return $Yes;
    }

    # This lock file prevents multiple running instances of Pipeline for one lane (e.g. one
    #   is run from cron and one manually). Non-blocking lock is requested only for read only
    #   operations.
    #
    my $global_lock_file = "$lane_path/$$self{'global_lock'}";
    my $locked_fh = lock_file($global_lock_file,$read_only ? 0 : 1);
    if ( !$locked_fh && !$read_only ) { $self->throw("The directory is being currently worked on: $global_lock_file\n"); }


    # If the program exits abruptly by $self->throw or die call, the lock file must be cleaned. 
    #   Otherwise, the program would always exit on this lane.
    if ( $locked_fh ) { $self->register_for_unlinking($global_lock_file); }
    
    # The main fun begins.
    my $all_done = $Yes;

    $self->debug("$lane_path\n");
    for my $action (@{$self->{'actions'}})
    {
        # If the option 'list_tasks' is set, only list available actions.
        if ( $$self{'list_tasks'} ) { $self->debug("$$action{'name'}\n"); next }

        # If the option 'task' is set, only the given task will be run.
        if ( $$self{'task'} && !$$self{'task'}{$$action{'name'}} ) { next }

        # First check what is the status of the task - is the task finished, that is,
        #   are all target files in place and are they in the expected state?
        #
        # Note, action_lock is not a real lock, it just stores the LSF job IDs.
        #
        my $action_lock = "$lane_path/$$self{'prefix'}$$action{'name'}.jids";
        my $status = $self->is_finished($lane_path,$action);
        if ( $status == $Yes ) 
        { 
            $self->debug("The task \"$$action{'name'}\" completed.\n");
            unlink($action_lock);
            next;
        }
        if ( $status & $Error )
        {
            $self->warn("The task \"$$action{'name'}\" ended with an error:\n\t$lane_path\n");
            $all_done = $Error;
            last;
        }


        # If we are here, the task needs to be run. Check if it is currently running.
        #
        $status = $self->running_status($action_lock);
        if ( $status & $LSF::Running ) 
        { 
            $self->debug("The task \"$$action{'name'}\" is running.\n");
            $all_done = $Running; 
            if ( ! $$self{'stepwise'} ) { next; }
            last;
        }
        if ( $status & $LSF::Error )
        {
            # If exit_on_errors is set to 0, the action will be rerun. If it finishes
            #   successfully, the lock file will be deleted above after the is_finished()
            #   call. Only in case it fails repeatedly, we'll be executing this again. 
            #   The script calling pipeline must handle repeated failures.
            #
            $self->warn("The task \"$$action{'name'}\" ended with an error - some of the LSF jobs returned wrong status:\n\t$action_lock\n");
            if ( $$self{'exit_on_errors'} )
            {
                $self->warn("When problem fixed, remove the lockfile and rerun...\n");
                $all_done = $Error;
                last;
            }
            if ( $$self{check_status} )
            {
                $all_done = $Error;
                last;
            }
        }
        $all_done = $No;


        # If we are here, the task is not running and should be run. Check if it CAN be
        #   run, i.e., if all prerequisites exist.
        #
        $status = $self->can_be_run($lane_path,$action);
        if ( $status == $No ) { last }

        # OK, run the action!
        #
        $self->debug("The task $$action{'name'} will be run...\n");
        if ( $$self{'check_status'} ) { last }

        eval { $status = &{$$action{'action'}}($self,$lane_path,$action_lock); };
        if ($@) 
        { 
            if ( $@ eq $VertRes::Base::SIGNAL_CAUGHT_EVENT ) 
            { 
                # The pipeline received SIGTERM or SIGINT.
                die $@; 
            }

            $self->warn("The action \"$$action{'name'}\" failed.\n"); 
            $self->warn($@); 
            $all_done = $Error; 
            last; 
        }

        # Some actions are done immediately and return $Yes. For those, run next action.
        last unless (!$$self{'stepwise'} || ($status && $status==$Yes) );
    }

    if ( $locked_fh )
    {
        unlock_file($global_lock_file,$locked_fh);
        $self->unregister_for_unlinking($global_lock_file);
    }

    return $all_done;
}


sub what_files_are_missing
{
    my ($self, $path, $files) = @_;

    my @missing = ();
    for my $file (@$files)
    {
        my $file_path = index($file, '/') == 0 ? $file : "$path/$file";
        if ( !(-e $file_path) )
        {
            push @missing, $file_path;
        }
    }
    return \@missing;
}

=head2 clean

        Description : The default routine cleans all lock files. The pipelines should implement the cleaning.
        Returntype  : None

=cut

sub clean
{
    my ($self) = @_;

    if ( !$$self{'lane_path'} ) { $self->throw("Missing parameter: the lane to be cleaned.\n"); }

    for my $action (@{$self->{'actions'}})
    {
        my $action_lock = "$$self{'lane_path'}/$$self{'prefix'}$$action{'name'}.jids";
        if ( ! -e $action_lock ) { next }
        $self->debug("unlink $action_lock\n");
        unlink($action_lock) or $self->throw("Could not unlink $action_lock: $!");
    }

    return;
}


=head2 is_finished

        Description : The default routine checks if all files provided by the action are present.
        Arg [1]     : the lane to be processed (absolute path)
        Arg [2]     : the 'action' hash for the task in question
        Returntype  : $YES if the lane has finished, $No when unfinished.

=cut

sub is_finished
{
    my ($self,$lane_path,$action) = @_;

    if ( $$self{'rerun'} ) { return $No; }

    # First, check if all files are in place, if so, we are done. If more
    #   sofisticated checking needs to be done (like checking the consistency of
    #   the files), rewrite this routine in child.
    #
    # When the empty list is returned by 'provides', what_files_are_missing
    #   will return an empty list. This is interpreted that nothing is required
    #   and the task is finished.
    #
    # When 0 or '' is returned by 'provides', it is assumed that the task is not
    #   finished and must be run. It is the caller responsibility to ensure that
    #   the task will not be run again and again.
    #
    my $provides = &{$$action{'provides'}}($self,$lane_path);
    if ( !$provides ) { return $No; }

    my $missing = $self->what_files_are_missing($lane_path,$provides);

    if ( !scalar @$missing ) { return $Yes; }

    $self->debug("Some of the files provided by \"$$action{'name'}\" not in place:\n\t".join("\n\t",@$missing)."\n");

    return $No;
}


=head2 can_be_run

        Description : The default routine checks if all files required by the action are present.
        Arg [1]     : the lane to be processed (absolute path)
        Arg [2]     : the 'action' hash for the task in question
        Returntype  : $YES if the lane has all prerequisites, $No if something is missing.

=cut

sub can_be_run
{
    my ($self,$lane_path,$action) = @_;

    # First, check if all files are in place, if so, we are done. If more
    #   sofisticated checking needs to be done (like checking the consistency of
    #   the files), rewrite this routine in child.
    #
    my $requires = &{$$action{'requires'}}($self,$lane_path);
    my $missing  = $self->what_files_are_missing($lane_path,$requires);

    if ( !scalar @$missing ) { return $Yes; }

    $self->debug("Some of the files required by \"$$action{'name'}\" not in place:\n\t".join("\n\t",@$missing));

    return $No;
}


=head2 running_status

        Description : The default routine checks the lock file for LSF job ID's and checks their status.
        Arg [1]     : the lock file of the task
        Returntype  : LSF::is_job_running return codes.

=cut

sub running_status
{
    my ($self,$lock_file) = @_;

    my $status;
    eval { $status = LSF::is_job_running($lock_file); };
    if ( $@ )
    {
        if ( $@ eq $VertRes::Base::SIGNAL_CAUGHT_EVENT ) 
        { 
            # The pipeline received SIGTERM or SIGINT.
            die $@; 
        }
        return $LSF::Error;
    }
    return $status;
}


sub lock_file
{
    my ($lock_file,$block) = @_;

    my $operation = $block ? LOCK_EX : LOCK_EX|LOCK_NB;

    sysopen(my $fh, $lock_file, O_RDWR|O_CREAT) or confess "$lock_file: $!";
    my $locked = flock($fh, $operation);

    if ( !$locked ) 
    {
        close($fh);
        return 0;
    }

    return $fh;
}


# So far, no checking is being done.
sub unlock_file
{
    my ($lock_file,$lock_fh) = @_;
    if ( -e $lock_file ) { unlink $lock_file; }
    flock($lock_fh,LOCK_UN);
    close($lock_fh) or confess "close $lock_fh: $!";
}

1;

