=head1 NAME

Runner.pm   - A simple module for quick development of scripts and pipelines which can be run in both serial and parallel mode.

=head1 SYNOPSIS

    # The code "test-runner" below shows a simple pipeline which creates
    # three files in your home directory (named "Hello.1", "Hello.2", and "Hello.3").
    # When all files are created, the message "All done!" will be printed. The
    # pipeline can be run in
    #   - crontab mode (exits at checkpoints when some of the files are not finished)
    #       test-runner
    #   - daemon mode (waits at checkpoints with 1 minute sleep intervals)
    #       test-runner +loop 60
    #   - serial mode (jobs are not submitted to LSF but are run locally)
    #       test-runner +local
    #
    # The test-runner code:

    #!/usr/bin/env perl
    use strict;
    use warnings;
    
    # Create new runner object and run it
    my $runner = myRunner->new();
    $runner->run();
    
    exit;
    
    #------------------------
    
    package myRunner;
    use base qw(Runner);
    use strict;
    use warnings;
    
    # The user must define at least this method
    sub main
    {
        my ($self) = @_;
        for my $file qw(1 2 3)
        {
            # When run in parallel mode (default), the jobs will be submitted
            #   to farm by the spawn call. The arguments are: 
            #   - method    .. subroutine to be called (defined by the user)
            #   - done_file .. the file to be created by the method
            #   - params    .. arbitrary number of arguments which will be passsed to the method

            my $method    = "touch";
            my $done_file = "$ENV{HOME}/Hello.$file";
            my @params    = ();
            $self->spawn($method,$done_file,@params);
        }
        # Checkpoint, wait until all the above files are finished
        $self->wait;

        print STDERR "All done!\n";
        $self->all_done;
    }

    sub touch
    {
        my ($self,$file) = @_;
        print STDERR "touch $file\n";
        # Pretend that it takes some time to complete this task
        sleep(10);
        `touch $file`;
    }

    sub new
    {
        my ($class,@args) = @_;
        my $self = $class->SUPER::new(@args);
        return $self;
    }

=head1 METHODS

=cut

package Runner;
use strict;
use warnings;
use Carp;
use Storable qw(nstore retrieve dclone);
use File::Temp;
use Data::Dumper;

sub new
{
    my ($class,@args) = @_;
    my $self = @args ? {@args} : {};
    bless $self, ref($class) || $class;
    $$self{_status_codes}{DONE} = 1;
    $$self{_farm} = 'LSF';
    $$self{_farm_options} = { runtime=>600 };
    $$self{_running_jobs} = {};
    $$self{_nretries} = 1;
    $$self{usage} = 
        "Runner.pm arguments:\n" .
        "   +help               Summary of commands\n" .
        "   +config <file>      Configuration file\n" .
        "   +local              Do not submit jobs to LSF, but run serially\n" .
        "   +loop <int>         Run in daemon mode with <int> sleep intervals\n" .
        "   +maxjobs <int>      Maximum number of simultaneously running jobs\n" .
        "   +retries <int>      Maximum number of retries. When negative, the runner eventually skips the task rather than exiting completely. [$$self{_nretries}]\n" .
        "   +run <file>         Run the freezed object created by spawn\n" .
        "   +sampleconf         Print a working configuration example\n" .
        "   +show <file>        Print the content of the freezed object created by spawn\n" .
        "   +verbose            Print debugging messages\n" .
        "\n";
    return $self;
}

=head2 run

    About : The main runner method which parses runner's command line parameters and calls the main() method defined by the user.
    Args  : The system command line options are prefixed by "+" to distinguish from user-module options and must come first, before the user-module options 
                +help
                    Summary of commands
                +config <file>
                    Optional configuration file for overriding defaults
                +local
                    Do not submit jobs to LSF, but run serially
                +loop <int>
                    Run in daemon mode with <int> sleep intervals
                +maxjobs <int>
                    Maximum number of simultaneously running jobs
                +retries <int>
                    Maximum number of retries
                +run <file>
                    Run the freezed object created by spawn
                +sampleconf
                    Print a working config file example
                +show <file>
                     Print the content of the freezed object created by spawn
                +verbose   
                    Print debugging messages
                
=cut

sub run
{
    my ($self) = @_;

    # Parse runner system parameters
    while (defined(my $arg=shift(@ARGV)))
    {
        if ( $arg eq '+help' ) { $self->throw(); }
        if ( $arg eq '+config' ) { $self->_read_config(shift(@ARGV)); next; }
        if ( $arg eq '+sampleconf' ) { $self->_sample_config(); next; }
        if ( $arg eq '+loop' ) { $$self{_loop}=shift(@ARGV); next; }
        if ( $arg eq '+maxjobs' ) { $$self{_maxjobs}=shift(@ARGV); next; }
        if ( $arg eq '+retries' ) { $$self{_nretries}=shift(@ARGV); next; }
        if ( $arg eq '+verbose' ) { $$self{_verbose}=1; next; }
        if ( $arg eq '+local' ) { $$self{_run_locally}=1; next; }
        if ( $arg eq '+show' ) 
        { 
            $arg = shift(@ARGV);
            my $obj = retrieve($arg); 
            print Dumper($obj);
            exit;
        }
        if ( $arg eq '+run' ) 
        { 
            $arg = shift(@ARGV);
            $self->_revive($arg);
            exit;
        }
        unshift(@ARGV,$arg);
        last;
    }

    # Run the user's module once or multiple times
    while (1)
    {
        my $pid = fork();
        if ( !$pid ) { $self->main(); }
        else
        {
            # If killed, killed the child process too.
            $SIG{TERM} = $SIG{INT} = sub { kill 9,$pid; die "Signal caught, killing the child $pid.\n"; };
        }
        wait();
        if ( !$$self{_loop} or $? ) { return; }
        $self->debugln("sleeping...");
        sleep($$self{_loop});
    }
}


# Read the user config file, the values set there will override variables set in the user's clone of Runner.
# Note that developers could be easily forced to document *all* options if desired.
sub _read_config
{
    my ($self,$config) = @_;

    if ( !exists($$self{_sampleconf}) or !length($$self{_sampleconf}) )
    {
        $self->throw("No config file parameters are accepted by this script.");
    }

    my %x = do "$config";
    while (my ($key,$value) = each %x)
    {
        if ( !ref($value) ) 
        { 
            $$self{$key} = $value;
            next;
        }
        $$self{$key} = dclone($value);
    }
    
    $$self{_mtime} = $self->mtime($config);
}

=head2 mtime

    About : Stat the mtime of a file
    Usage : $self->mtime('/path/to/file');
    Args  : path to file

=cut

sub mtime
{
    my ($self, $file) = @_;
    my $mtime = (stat($file))[9];
    $mtime || $self->throw("Could not stat mtime for file $file\n");
    return $mtime;
}

sub _sample_config
{
    my ($self) = @_;
    if ( !exists($$self{_sampleconf}) )
    {
        $self->throw("Sample config not available, the '_sampleconf' key not set. This should be fixed!\n");
    }
    if ( !length($$self{_sampleconf}) )
    {
        print "# No config file parameters accepted by this script.\n";
        $self->all_done;
    }
    print $$self{_sampleconf};
    $self->all_done;
}

=head2 set_limits

    About : Set time and memory requirements for computing farm
    Usage : $self->set_limits(memory=>1_000, runtime=>24*60);
    Args  : <memory>
                Expected memory requirements [MB] or undef to unset
            <runtime>
                Expected running time [minutes] or undef to unset
                
=cut

sub set_limits
{
    my ($self,%args) = @_;
    $$self{_farm_options} = { %{$$self{_farm_options}}, %args };
}


=head2 get_limits

    About : get limits set for computing farm
    Usage : $self->get_limits(memory);
    Args  : See set_limits
                
=cut

sub get_limits
{
    my ($self,$arg) = @_;
    return exists($$self{_farm_options}{$arg}) ? $$self{_farm_options}{$arg} : undef;
}


=head2 spawn

    About : Schedule a job for execution. When +maxjobs limit is exceeded, the runner will exit via self->wait call.
    Usage : $self->spawn("method",$done_file,@params);
    Args  : <func_name>
                The method to be run
            <file>
                The file to be created by the method. If exists, the task is completed and
                spawn returns immediately.
            <array>
                Arbitrary number of parameters to be passed to the method. Note that these
                are stored by Storable and thus the same limitations apply. Passing for example,
                complex objects with CODE refs will not work.
    Returns : 0
                Job was submitted to farm
              1
                Job is finished
                
=cut

sub spawn
{
    my ($self,$call,@args) = @_;

    if ( !$self->can($call) ) { $self->throw("No such method: [$call]\n"); }

    # Register the file for the next checkpoint
    my $done_file = $args[0];
    push @{$$self{_checkpoints}}, $done_file;

    # If the file is there, no need to run anything
    if ( $self->is_finished($done_file) ) { return 1; }

    # If the file needs to be skipped, then skip it
    my $basename = $self->_get_temp_prefix($done_file);
    if ( -e $basename . '.s' )
    {
        # This is currently the only way to clean the skip files: run with +retries set to positive value
        if ( $$self{_nretries}<0 ) { return 1; }
        $self->debugln("Cleaning skip file: $basename.s");
        unlink($basename . '.s');
    }

    # Store all necessary information to run the task in a temp file
    $$self{_store}{call} = $call;
    $$self{_store}{args} = \@args;
    $$self{_store}{done_file} = $done_file;
    my $tmp_file = $basename . '.r';
    nstore($self,$tmp_file);

    # With '+local', the jobs will be run serially
    if ( $$self{_run_locally} ) 
    {
        my $cmd = qq[$0 +run $tmp_file];
        $self->debugln("$call:\t$cmd");
        system($cmd);
        return 1;
    }

    # Otherwise submit to farm 
    else
    {
        # Check if the number of running jobs should be kept low. In case there are too many jobs already, let through
        #   only jobs which previously failed, i.e. are registered as running.
        if ( exists($$self{_maxjobs}) && scalar keys %{$$self{_running_jobs}} >= $$self{_maxjobs} && !exists($$self{_running_jobs}{$done_file}) )
        {
            $self->wait;
            return 1;
        }
        $$self{_running_jobs}{$done_file} = 1;

        $self->_spawn_to_farm($tmp_file);
        return 0;
    }
}

sub _spawn_to_farm
{
    my ($self,$freeze_file) = @_;

    my $farm = $$self{_farm};
    eval "require $farm ";
    if ( $@ ) { $self->throw("require $farm\n$@"); }

    my $prefix = $self->_get_temp_prefix($$self{_store}{done_file});
    my $farm_jobs_ids = $prefix.'.jid';
    my $status = $farm->can('is_job_running')->($farm_jobs_ids);

    # If the file is already running, return.
    if ( $status & eval "\$${farm}::Running" ) { return; }

    # If the file was already run and failed, check if it failed repeatedly
    if ( $status & eval "\$${farm}::Error" ) 
    { 
        my ($nfailures) = `wc -l $farm_jobs_ids`;
        $nfailures =~ s/\s+.+$//;
        if ( $nfailures > abs($$self{_nretries}) )
        {   
            my $msg = "The job failed repeatedly: $prefix.[oers], $$self{_store}{call}(" .join(',',@{$$self{_store}{args}}). ")"
                ."\n(Remove $prefix.jid to clean the status, increase +retries or run with negative value of +retries to skip this task.)\n";

            if ( $$self{_nretries}>=0 )
            {
                $self->throw($msg);
            }
            $self->warn($msg,"This was last attempt, excluding from normal flow. (Remove $prefix.s to clean the status.)\n");
            system("touch $prefix.s");
            return;
        }
        else
        {
            $self->warn("Running again, the previous attempt failed: $prefix.[oer], $$self{_store}{call}(" .join(',',@{$$self{_store}{args}}). ")\n");
        }
    }

    # Run the job
    my $cmd = qq[$0 +run $freeze_file];
    $self->debugln("$$self{_store}{call}:\t$cmd");
    my $ok;
    eval {
        $farm->can('run')->($farm_jobs_ids,'.',$prefix,$$self{_farm_options},$cmd);
        $ok = 1;
    };
    if ( !$ok )
    {
        if ( $$self{_nretries}<0 )
        {
            $self->warn($@,"This was last attempt, excluding from normal flow. (Remove $prefix.s to clean the status.)\n");
            system("touch $prefix.s");
            return;
        }
        $self->throw($@);
    }
    return;
}

=head2 wait

    About : Checkpoint, wait for all tasks to finish. 
    Usage : $self->spawn("method",$done_file1,@params1); 
            $self->spawn("method",$done_file2,@params2);
            $self->wait();
    Args  : <none>
                Without arguments, waits for files registered by previous spawn calls.
            <@files>
                Extra files to wait for, in addition to those registered by spawn.
=cut

sub wait
{
    my ($self,@files) = @_;

    for my $file (@{$$self{_checkpoints}},@files)
    {
        if ( ! $self->is_finished($file) ) 
        { 
            my $prefix = $self->_get_temp_prefix($file);
            $self->debugln("The job not finished: $prefix.*");
            exit; 
        }
    }
    $$self{_checkpoints} = [];
}


=head2 all_done

    About : Exit with "all done" status
    Usage : $self->all_done();
    Args  : None
                
=cut

sub all_done
{
    my ($self) = @_;
    $self->debugln("All done!");
    exit $$self{_status_codes}{DONE};
}

=head2 clean

    About : Clean all system files (in .jobs directories)
    Usage : $self->clean($dir);
    Args  : <@dirs>
                Directories to recursively clean from all .jobs subdirs
=cut

sub clean
{
    my ($self,@dirs) = @_;
    for my $dir (@dirs)
    {
        my $cmd = "find $dir -name .jobs | xargs rm -rf";
        $self->debugln($cmd);
        system($cmd);
    }
}

# Do not advertise this method, the user module should not need it
#
#   =head2 is_finished
#   
#       About : Check if the file is finished.
#       Usage : $self->is_finished('some/file');
#       Args  : <file list>
#                   The name of the file to check the existence of
#                   
#   =cut

sub is_finished
{
    my ($self,@files) = @_;
    my $all_finished = 1;
    for my $file (@files)
    {
        my $is_finished = -e $file;
        if ( !$is_finished ) { $all_finished=0; }
        elsif ( exists($$self{_running_jobs}{$file}) ) { delete($$self{_running_jobs}{$file}); }
    }
    return $all_finished;
}

# Run the freezed object created by spawn
sub _revive
{
    my ($self,$freeze_file) = @_;

    my $back = retrieve($freeze_file);
    if ( $$self{clean} ) { unlink($freeze_file); }

    while (my ($key,$value) = each %$back) { $$self{$key} = $value; }

    my $code = $self->can($$self{_store}{call});
    &$code($self,@{$$self{_store}{args}});

    # If we are here, the code finished successfully - remove the jids file to wipe out history
    my $prefix = $self->_get_temp_prefix($$self{_store}{done_file});
    my $farm_jobs_ids = $prefix.'.jid';
    unlink($farm_jobs_ids);
}

# Create a temporary prefix using the given template and create a temporary directory
sub _get_temp_prefix
{
    my ($self,$fname) = @_;
    if ( !($fname=~m{(^.+/)?([^/]+)$}) ) { $self->throw("FIXME: could not parse [$fname]\n"); }
    my $dir  = (defined $1 ? $1 : '')  . '.jobs';
    my $file = $2;
    if ( ! -d $dir ) 
    { 
        `mkdir -p $dir`; 
        if ( $? ) { $self->throw("Cannot create directory [$dir]: $!"); }
    }
    return "$dir/$file";
}

=head2 throw

    About : Throws an error.
    Args  : <array>
                The message to be printed. If no message is given, the usage will be printed instead.

=cut

sub throw
{
    my ($self,@msg) = @_;
    if ( scalar @msg ) { confess @msg; }
    die $$self{usage};
}

=head2 warn

    About : Print a warning message.
    Args  : <array>
                The message.

=cut

sub warn
{
    my ($self,@msg) = @_;
    print STDERR @msg;
}

=head2 debugln

    About : When the runner is run with the "+verbose" command line option, the debugging message will printed with newline appended.
    Args  : <array>
                The message.

=cut

sub debugln
{
    my ($self,@msg) = @_;
    if ( !$$self{_verbose} ) { return; }
    print STDERR @msg , "\n";
}

=head1 AUTHORS

petr.danecek@sanger

=cut

1;
