=head1 NAME

VertRes::Pipeline 

=head1 SYNOPSIS

use VertRes::Pipeline;

=head1 DESCRIPTION



=head1 AUTHOR

Petr Danecek pd3@sanger.ac.uk

=cut

package VertRes::Pipeline;

use strict;
use warnings;
use Carp;   # this is for confess in lock_file and unlock_file. eventually should go somewhere else

use base qw(VertRes::Base);
use VertRes::LSF;
use File::Spec;
use File::Copy;
use File::Basename;
use Fcntl qw(:DEFAULT :flock);
use VertRes::Utils::FileSystem;
use Bio::VertRes::Permissions::ModifyPermissions;
use Bio::VertRes::Permissions::Groups;

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
    bless($self,$class);
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
        for my $task (@tasks) 
        { 
            if ( !$self->_has_action($task) ) { $self->throw(qq[The task "$task" not recognised.\n]); }
            $$self{'task'}{$task}=1; 
        }
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
        if ( $status & $VertRes::LSF::Running ) 
        { 
            $self->debug("The task \"$$action{'name'}\" is running.\n");
            $all_done = $Running; 
            if ( ! $$self{'stepwise'} ) { next; }
            last;
        }
        if ( $status & $VertRes::LSF::Error )
        {
            # If exit_on_errors is set to 0, the action will be rerun. If it finishes
            #   successfully, the lock file will be deleted above after the is_finished()
            #   call. Only in case it fails repeatedly, we'll be executing this again. 
            #   The script calling pipeline must handle repeated failures.
            #
            if ( $$self{check_status} )
            {
                $all_done = $Error;
                last;
            }
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
            if ( $@ eq $VertRes::Base::SIGNAL_CAUGHT_EVENT ) 
            { 
                # The pipeline received SIGTERM or SIGINT.
                die $@; 
            }

            $self->warn("The action \"$$action{'name'}\" failed:\n$@"); 
            $all_done = $Error; 
            last; 
        }

        # Some actions are done immediately and return $Yes. For those, run next action.
        last unless (!$$self{'stepwise'} || (defined $status && $status==$Yes) );
        $self->debug("The task \"$$action{'name'}\" completed.\n");
        unlink($action_lock);
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
                        If task is set, all files listed by call of 'provides' for this task will be deleted.
                        This will effectively force to rerun the task once next time.
        Returntype  : None

=cut

sub clean
{
    my ($self) = @_;

    if ( !$$self{'lane_path'} ) { $self->throw("Missing parameter: the lane to be cleaned.\n"); }

    for my $action (@{$self->{'actions'}})
    {
        my $action_lock = "$$self{'lane_path'}/$$self{'prefix'}$$action{'name'}.jids";
        if ( -e $action_lock ) 
        { 
            $self->debug("unlink $action_lock\n");
            unlink($action_lock) or $self->throw("Could not unlink $action_lock: $!");
        }

        # Clean the files for each listed task.
        if ( $$self{task} && $$self{'task'}{$$action{name}} && (my $provides=&{$$action{'provides'}}($self,$$self{lane_path})) )
        {
            for my $file (@$provides)
            {
                my $file_path = index($file, '/') == 0 ? $file : "$$self{lane_path}/$file";
                if ( !(-e $file_path) ) { next; }
                $self->debug("unlink $file_path\n");
                unlink($file_path) or $self->throw("Could not unlink $file_path: $!");
            }
        }
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

    $self->debug("Some of the files required by \"$$action{'name'}\" not in place:\n\t".join("\n\t",@$missing)."\n");

    return $No;
}


=head2 running_status

        Description : The default routine checks the lock file for LSF job ID's and checks their status.
        Arg [1]     : the lock file of the task
        Returntype  : VertRes::LSF::is_job_running return codes.

=cut

sub running_status
{
    my ($self,$lock_file) = @_;

    my $status;
    eval { $status = VertRes::LSF::is_job_running($lock_file); };
    if ( $@ )
    {
        if ( $@ eq $VertRes::Base::SIGNAL_CAUGHT_EVENT ) 
        { 
            # The pipeline received SIGTERM or SIGINT.
            die $@; 
        }
        return $VertRes::LSF::Error;
    }
    return $status;
}


sub lock_file
{
    my ($lock_file,$block) = @_;

    # A change was requested: never block.
    #   my $operation = $block ? LOCK_EX : LOCK_EX|LOCK_NB;
    my $operation = LOCK_EX|LOCK_NB;

    # check the disk isn't full (below sysopen doesn't fail when disk is full;
    # it just hangs instead)
    my (undef, $lock_dir) = fileparse($lock_file);
    my $bytes_left = VertRes::Utils::FileSystem->new->disk_available($lock_dir);
    unless ($bytes_left > 1000) {
        warn "got $bytes_left bytes free on $lock_dir\n";
        return 0;
    }

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

# For a given bsub job name (arg 3 to VertRes::LSF::run), if the bsub o/e files exist
# they will be moved to .previous files. If a .previous file exists, its
# content will be moved to an .archive file. Returns the abs path to the
# .previous error file so you can parse it for errors before proceeding.
#
# With extra option set to true, will move all data from any existing current
# or previous bsub file to the .archive file: for use when an action completed
# successfully and you want to tidy up
sub archive_bsub_files {
    my ($self, $lane_path, $job_name, $final) = @_;
    
    my $prev_error;
    foreach my $suffix ('e', 'o') {
        my $bsub_file = File::Spec->catfile($lane_path, $job_name.'.'.$suffix);
        my $previous_bsub_file = $bsub_file.'.previous';
        my $archive_bsub_file = File::Spec->catfile($lane_path, '.'.$job_name.'.'.$suffix.'.archive');
        
        if (-s $bsub_file) {
            # don't archive files greater than 1GB
            if (-e $previous_bsub_file && -s $previous_bsub_file < 1000000000) {
                $self->_move_file_content($previous_bsub_file, $archive_bsub_file);
            }
            else {
                unlink($previous_bsub_file);
            }
            
            move($bsub_file, $previous_bsub_file) || $self->throw("Couldn't move $bsub_file to $previous_bsub_file");
        }
        
        if ($final) {
            # don't archive files greater than 1GB
            if (-e $previous_bsub_file && -s $previous_bsub_file < 1000000000) {
                $self->_move_file_content($previous_bsub_file, $archive_bsub_file);
            }
            else {
                unlink($previous_bsub_file);
            }
        }
        
        if ($suffix eq 'e' && -s $previous_bsub_file) {
            $prev_error = $previous_bsub_file;
        }
    }
    
    return $prev_error;
}

sub _move_file_content {
    my ($self, $source, $dest) = @_;
    
    if (-s $source) {
        open(my $ifh, $source) || $self->throw("Could not open '$source'");
        open(my $ofh, '>>', $dest) || $self->throw("Could not append to '$dest'");
        while (<$ifh>) {
            print $ofh $_;
        }
        close($ofh);
        close($ifh);
        
        unlink($source);
    }
}


sub _has_action
{
    my ($self,$task) = @_;
    if ( $task eq 'clean' ) { return 1; }
    for my $action (@{$self->{'actions'}})
    {
        if ( $$action{name} eq $task ) { return 1; }
    }
    return 0;
}


sub get_unix_group
{
	my ( $self ) = @_;
	return undef if(!(defined($$self{unix_group}) && defined($$self{octal_permissions}) ));
	
	my $unix_group = $$self{unix_group};
	my $vrtrack = VRTrack::VRTrack->new( $$self{db} ) or $self->throw("Could not connect to the database\n");
	my $vrlane = VRTrack::Lane->new_by_name( $vrtrack, $$self{lane} ) or $self->throw("No such lane in the DB: [$$self{lane}]\n");
	my %lane_objs = $vrtrack->lane_hierarchy_objects($vrlane);
	if(defined($lane_objs{project}) )
    {
        if(defined($lane_objs{project}->data_access_group) )
        {
            my $obj = Bio::VertRes::Permissions::Groups->new();
            my $group = $obj->is_member_of_group($lane_objs{project}->data_access_group);
            $unix_group = $group if(defined($group));
        }
    }
	return $unix_group;
}

sub update_file_permissions
{
	my ( $self, $lane_path ) = @_;
	return unless(defined($$self{octal_permissions}));
	
	my $unix_group = $self->get_unix_group;
	if(defined($unix_group) )
	{
        my $change_permissions_obj = Bio::VertRes::Permissions::ModifyPermissions->new(
            input_directories => [$lane_path],
            group             => $unix_group,
		        threads           => 0,
            octal_permissions => $$self{octal_permissions});
        $change_permissions_obj->update_permissions;	
    }
}

sub umask_str
{
	my ( $self ) = @_;
	my $umask_string = "";
	if($$self{umask})
	{
		$umask_string = ' umask '.$$self{umask}."; ";
	}
	return $umask_string;
}


1;

