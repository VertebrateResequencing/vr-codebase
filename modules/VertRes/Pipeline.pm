package VertRes::Pipeline;

use strict;
use warnings;
use Carp;   # this is for confess in lock_file and unlock_file. eventually should go somewhere else

use base qw(VertRes::Base);
use LSF;
use Fcntl qw(:DEFAULT :flock);

our $Yes   = 0;
our $No    = 1;
our $Error = 2;

our $default_options =
{
    # General options, common to all modules
    'check_status'   => 0,
    'clean'          => 0,
    'exit_on_errors' => 1,
    'global_lock'    => 'Pipeline.lock',
    'lane'           => '',
    'list_tasks'     => 0,
    'log'            => '_log',
    'prefix'         => '_',
    'rerun'          => 0,
    'stepwise'       => 1,
    'task'           => '',

    'Error'          => $Error,
    'Yes'            => $Yes,
    'No'             => $No,
};


=head2 new

        Options    :
                    check_status    .. if set to 1, nothing will be run, only the status of the lane will be reported
                    clean           .. delete all files provided by actions
                    exit_on_errors  .. if set to 0, an unfinished task will be run again even if the LSF job previously failed
                    global_lock     .. lock file for the lane, to prevent multiple run_lane subroutines running
                    lane            .. optional, can be overrided by run_lane
                    list_tasks      .. don't do anything, only list available tasks
                    log             .. this is used when logfile is not given (may be both relative and absolute path)
                    prefix          .. to distinguish between multiple variants of the same module, but running with different options
                    rerun           .. if set to 1, always run again. Implies that the task option is set.
                    stepwise        .. execute only one task at time
                    task            .. if set, execute only the given tasks (separated by commas, the order is ignored)
        Example    : See TrackDummy.pm for "Hello World" example and TrackQC.pm for real-life example.

=cut

sub new 
{
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(%$default_options, @args);
    return $self;
}


=head2 run_lanes

        Arg [1]     : file with a list of lanes to be processed
        Description : A convenience wrapper for the run_lane routine - not used.
        Example    : See test-qc and run-qc-pipeline for an example.
        Returntype : $YES if all of the lanes have finished, $No when some is unfinished and $Error if an error occured.

=cut

sub run_lanes
{
    my ($self,$lane_list_file) = @_;

    my $done = $Yes;

    open(my $fh,'<',$lane_list_file) or $self->throw("$lane_list_file: $!");
    while (my $lane=<$fh>)
    {
        if ( $lane=~/^\s*$/ ) { next }

        chomp($lane);
        my $status = $self->run_lane($lane);
        if ( $status==$Error ) { $done=$Error }
        elsif ( $status==$No && $done!=$Error ) { $done=$No; }
    }
    close($fh) or $self->throw("$lane_list_file: $!");

    return $done;
}


=head2 run_lane

        Arg [1]    : the lane to be processed (absolute path). Can be ommited when supplied with new().
        Example    : See TrackDummy.pm for "Hello World" example and TrackQC.pm for real-life example.
        Returntype : $YES if the lane has finished, $No when unfinished, $Error if an error occured.

=cut

sub run_lane
{
    my ($self,$lane_path) = @_;

    $lane_path = $$self{'lane'} unless $lane_path;
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

    
    # The cleaning should be called directly, not via run_lane.
    if ( $$self{clean} )
    {
        $self->clean();
        return $Yes;
    }

    # This lock file prevents multiple running instances of Pipeline for one lane (e.g. one
    #   is run from cron and one manually). 
    my $global_lock_file = "$lane_path/$$self{'global_lock'}";
    my $locked_fh = lock_file($global_lock_file);
    if ( !$locked_fh ) {  $self->throw("The directory is being currently worked on: $global_lock_file\n"); }


    # If the program exits abruptly by $self->throw or die call, the lock file must be cleaned. 
    #   Otherwise, the program would always exit on this lane.
    $self->register_for_unlinking($global_lock_file);
    
    
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
        my $action_lock = "$lane_path/$$self{'prefix'}$$action{'name'}.lock";
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
            $all_done = $No; 
            next unless $$self{'stepwise'};
            last;
        }
        if ( $status & $LSF::Error )
        {
            $self->warn("The task \"$$action{'name'}\" ended with an error - some of the LSF jobs returned wrong status:\n\t$action_lock\n");
            if ( $$self{'exit_on_errors'} )
            {
                $self->warn("When problem fixed, remove the lockfile and rerun...\n");
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
            $self->warn("The action \"$$action{'name'}\" failed.\n"); 
            $self->warn($@); 
            $all_done = $Error; 
            last; 
        }

        # Some actions are done immediately. For those, run next action.
        last unless (!$$self{'stepwise'} || $status==$Yes);
    }

    unlock_file($global_lock_file,$locked_fh);
    $self->unregister_for_unlinking($global_lock_file);

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

    if ( !$$self{'lane'} ) { $self->throw("Missing parameter: the lane to be cleaned.\n"); }

    for my $action (@{$self->{'actions'}})
    {
        my $action_lock = "$$self{'lane'}/$$self{'prefix'}$$action{'name'}.lock";
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
    my $provides = &{$$action{'provides'}}($self,$lane_path);
    my $missing  = $self->what_files_are_missing($lane_path,$provides);

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
    if ( $@ ) { return $LSF::Error }
    return $status;
}


sub lock_file
{
    my ($lock_file) = @_;

    sysopen(my $fh, $lock_file, O_RDWR|O_CREAT) or confess "$lock_file: $!";
    my $locked = flock($fh, LOCK_EX|LOCK_NB);

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
    close($lock_fh) or confess "close $lock_fh: $!";
    if ( -e $lock_file ) { unlink $lock_file; }
}

1;

