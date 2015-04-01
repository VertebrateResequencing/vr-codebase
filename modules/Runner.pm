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
        for my $file (qw(1 2 3))
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

        # Clean temporary files
        $self->clean;

        print STDERR "Mission accomplished!\n";
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
use Cwd;

sub new
{
    my ($class,@args) = @_;
    my $self = @args ? {@args} : {};
    bless $self, ref($class) || $class;
    $$self{_status_codes}{DONE}  = 111;
    $$self{_status_codes}{WAIT}  = 0;
    $$self{_status_codes}{ERROR} = 255;
    $$self{_farm} = 'LSF';
    $$self{_farm_options} = {};
    $$self{_running_jobs} = {};
    $$self{_nretries} = 1;
    $$self{_verbose} = 1;
    $$self{usage} = 
        "Runner.pm arguments:\n" .
        "   +help                   Summary of commands\n" .
        "   +config <file>          Configuration file\n" .
        "   +debug <file1> <file2>  Run the freezed object <file1> overriding with keys from <file2>\n" .
        "   +js <platform>          Job scheduler (lowercase allowed): LSF (bswitch), LSFCR (BLCR) [LSF]\n" .
        "   +kill                   Kill all running jobs\n" .
        "   +local                  Do not submit jobs to LSF, but run serially\n" .
        "   +lock <file>            Exit if another instance is already running\n" .
        "   +loop <int>             Run in daemon mode with <int> seconds sleep intervals\n" .
        "   +mail <address>         Email when the runner finishes\n" .
        "   +maxjobs <int>          Maximum number of simultaneously running jobs\n" .
        "   +nocache                When checking for finished files, do not rely on cached database and check again\n" .
        "   +retries <int>          Maximum number of retries. When negative, the runner eventually skips the task rather than exiting completely. [$$self{_nretries}]\n" .
        "   +run <file> <id>        Run the freezed object created by spawn\n" .
        "   +sampleconf             Print a working configuration example\n" .
        "   +show <file>            Print the content of the freezed object created by spawn\n" .
        "   +silent                 Decrease verbosity of the Runner module\n" .
        "\n";
    return $self;
}

=head2 run

    About : The main runner method which parses runner's command line parameters and calls the main() method defined by the user.
    Args  : The system command line options are prefixed by "+" to distinguish from user-module options and must come first, before the user-module options 
                +config <file>
                    Optional configuration file for overriding defaults
				+debug <file1> <file2>
					Run the freezed object <file1> overriding with keys from <file2>
                +help
                    Summary of commands
                +js <platform>
                    Job scheduler: LSF (bswitch to deal with overrun), LSFCR (BLCR)
                +kill
                    Kill all running jobs
                +local
                    Do not submit jobs to LSF, but run serially
                +lock <file>
                    Exit if another instance is already running
                +loop <int>
                    Run in daemon mode with <int> seconds sleep intervals.
                    Negative values can be used to request only one iteration
                    at a time and still indicate the interval, so that the job
                    scheduler can estimate if a job needs switching to a longer queue.
                +mail <address>
                    Email to send when the runner is done
                +maxjobs <int>
                    Maximum number of simultaneously running jobs
                +nocache
                    When checking for finished files, do not rely on cached data and check again
                +retries <int>
                    Maximum number of retries
                +run <file>
                    Run the freezed object created by spawn
                +sampleconf
                    Print a working config file example
                +show <file>
                    Print the content of the freezed object created by spawn
                +silent
                    Decrease verbosity of the Runner module
                
=cut

sub run
{
    my ($self) = @_;
    my @args = @ARGV;

    $$self{_about} = "Working directory: " . getcwd() . "\nCommand line: $0 " . join(' ',@args) . "\n";

    # Parse runner system parameters. Allow mixing + and - parameters
    my @argv = ();
    while (defined(my $arg=shift(@ARGV)))
    {
        if ( substr($arg,0,1) ne '+' ) { push @argv, $arg; next; }
        if ( $arg eq '+help' ) { $self->throw(); }
        if ( $arg eq '+config' ) { $self->_read_config(shift(@ARGV)); next; }
        if ( $arg eq '+sampleconf' ) { $self->_sample_config(); next; }
        if ( $arg eq '+loop' ) { $$self{_loop}=shift(@ARGV); next; }
        if ( $arg eq '+lock' ) { $$self{_lock}=shift(@ARGV); next; }
        if ( $arg eq '+kill' ) { $$self{_kill_jobs}=1; next; }
        if ( $arg eq '+maxjobs' ) { $$self{_maxjobs}=shift(@ARGV); next; }
        if ( $arg eq '+mail' ) { $$self{_mail}=shift(@ARGV); next; }
        if ( $arg eq '+nocache' ) { $$self{_nocache}=1; next; }
        if ( $arg eq '+retries' ) { $$self{_nretries}=shift(@ARGV); next; }
        if ( $arg eq '+verbose' ) { $$self{_verbose}=1; next; }
        if ( $arg eq '+silent' ) { $$self{_verbose}=0; next; }
        if ( $arg eq '+js' ) 
        { 
            $$self{_farm}=shift(@ARGV); 
            if ( $$self{_farm} eq 'lsf' ) { $$self{_farm} = 'LSF'; }
            elsif ( $$self{_farm} eq 'lsf-cr' ) { $$self{_farm} = 'LSFCR'; }
            elsif ( $$self{_farm} eq 'lsfcr' ) { $$self{_farm} = 'LSFCR'; }
            elsif ( $$self{_farm} eq 'LSF-CR' ) { $$self{_farm} = 'LSFCR'; }
            next; 
        }
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
            my $run_file = shift(@ARGV);
            my $job_id   = shift(@ARGV);
            if ( $job_id eq '--' )
            {
                $self->_revive($run_file);
            }
            else
            {
                $self->_revive_array($run_file, $job_id);
            }
            exit;
        }
        if ( $arg eq '+debug' ) 
        { 
            my $file1 = shift(@ARGV);
			my $file2 = shift(@ARGV);
            $self->_revive($file1,$file2);
            exit;
        }
    }
    @ARGV = @argv;

    $self->create_lock();

    # Run the user's module once or multiple times
    while (1)
    {
        # The parent loops and forks, the childs run main() in user's module
        my $pid = fork();
        if ( !$pid ) { $self->main(); return; }
        else
        {
            # If killed, kill the child process too.
            $SIG{TERM} = $SIG{INT} = sub { kill 9,$pid; die "Signal caught, killing the child $pid.\n"; };
        }
        wait();
        my $status = $?>>8;
        if ( $status ) 
        {
            $self->remove_lock();
            if ( $status==$$self{_status_codes}{DONE} ) { exit $status; }
            # Exit with the correct status immediately if the user module fails. Note that +retries applies only to spawned jobs.
            die "\n"; 
        }
        if ( !$$self{_loop} or $$self{_loop}<0 )
        { 
            $self->remove_lock();
            return; 
        }
        $self->debugln($$self{_about}, "sleeping for $$self{_loop} seconds...");
        sleep($$self{_loop});
    }
}

sub create_lock
{
    my ($self) = @_;
    if ( !exists($$self{_lock}) ) { return; }

    my $lock = $$self{_lock};
    if ( -e $lock )
    {
        # Find out the PID of the running pipeline
        open(my $fh,'<',$lock) or $self->throw("$lock; $!");
        while (my $pid=<$fh>)
        {
            chomp($pid);
            if ( !($pid=~/^\d+$/) ) { $self->throw("Could not parse lock file: $lock, $pid\n"); }

            # Is the process still running?
            my $running = kill 0, $pid;
            if ( $running )
            {
                $self->warn("\nAnother process is running ($pid), exiting.\n\n");
                exit $$self{_status_codes}{WAIT};
            }

            $self->warn("\nIgnoring an old lock file, PID $pid is not running.\n\n");
            last;
        }
        close($fh);
    }

    open(my $fh,'>',$lock) or usage(qq[$lock: $!]);
    print $fh $$ . "\n";
    close($fh);
}
sub remove_lock
{
    my ($self) = @_;
    if ( !exists($$self{_lock}) ) { return; }

    my $lock = $$self{_lock};
    if ( -e $lock ) { unlink($lock); }
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
	if ( !-e $config ) { $self->throw("The file does not exist: $config\n"); }

    open(my $fh,'<',$config) or $self->throw("$config: $!");
    my @config_lines = <$fh>;
    close($fh) or $self->throw("close failed: $config");
    my $config_str = join('',@config_lines);
    my $x = eval "{ $config_str }";
    if ( $@ ) { $self->throw("eval $config: $@\n"); }
    while (my ($key,$value) = each %$x)
    {
        if ( !ref($value) ) 
        { 
            $$self{$key} = $value;
            next;
        }
        $$self{$key} = dclone($value);
    }
    $$self{_config} = $config;
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
    Args  : <cpus>
                Number of processors to use
            <memory>
                Expected memory requirements [MB] or undef to unset
            <runtime>
                Expected running time [minutes] or undef to unset
            <queues>
                Hash with farm queue names (keys) and maximum runtime limits in
                seconds (values)
            <wakeup_interval>
                Make the job scheduler aware of the pipeline's polling interval
                [seconds] so that it can estimate if a job exceeds the runtime
                limit before next wake-up call. The command line parameter
                +loop overrides this value.

            Plus any job-scheduler specific options
                
=cut

sub set_limits
{
    my ($self,%args) = @_;
    $$self{_farm_options} = { %{$$self{_farm_options}}, %args };
    if ( exists($$self{_revived_file}) )
    {
        my $limits_fname = "$$self{_revived_file}.$$self{_revived_job}.limits";
        open(my $fh,'>',$limits_fname) or $self->throw("$limits_fname: $!");
        print $fh Dumper($$self{_farm_options});
        close($fh);
    }
}

# When revived, the user module can indirectly request increase of system limits in the next run, see set_limits
sub _update_limits
{
    my ($self,$wfile,$id) = @_;
    my $limits_fname = "$wfile.$id.limits";
    if ( !-e $limits_fname ) { return; }
    my $limits = do "$limits_fname";
    if ( $@ ) { $self->throw("do $limits_fname: $@\n"); }
    $self->inc_limits(%$limits);
}

=head2 inc_limits

    About : increase limits if lower than requested
    Usage : $self->inc_limits(memory=>10_000);
    Args  : See set_limits
                
=cut

sub inc_limits
{
    my ($self,%args) = @_;
    for my $key (keys %args)
    {
        if ( !exists($$self{_farm_options}{$key}) or $$self{_farm_options}{$key}<$args{$key} ) 
        { 
            $self->debugln("increasing limit, $key set to $args{$key}");
            $$self{_farm_options}{$key} = $args{$key};
        }
    }
}

=head2 get_limits

    About : get limits set for computing farm
    Usage : $self->get_limits('memory');
    Args  : See set_limits
                
=cut

sub get_limits
{
    my ($self,$arg) = @_;
    if ( ! defined $arg ) { return %{$$self{_farm_options}}; }
    return exists($$self{_farm_options}{$arg}) ? $$self{_farm_options}{$arg} : undef;
}

=head2 freeze

    About : freeze the runner object
    Usage : $self->freeze();
    Args  : <file>
                Targe checkpoint file name (optional)

=cut

sub freeze
{
    my ($self,$arg) = @_;
    my $rfile = $self->_get_temp_prefix($arg) . '.r';
    nstore($self,"$rfile.part");
    if ( -e $rfile )
    {
        my $md5_ori = `md5sum $rfile`;
        my $md5_new = `md5sum $rfile.part`;
        $md5_ori =~ s/\s+.*$//;
        $md5_new =~ s/\s+.*$//;
        if ( $md5_ori eq $md5_new )
        {
            unlink("$rfile.part");
            return $rfile;
        }
    }
    rename("$rfile.part",$rfile) or $self->throw("rename $rfile.part $rfile: $!");
    return $rfile;
}

=head2 spawn

    About : Schedule a job for execution. When +maxjobs limit is exceeded, the runner will exit via self->wait call.
                When the job fails too many times (controllable by the +retries option), the pipeline exists. With 
                negative value of +retries, a skip file '.s' is created and the job is reported as finished.
                The skip files are cleaned automatically when +retries is set to a positive value. A non cleanable 
                variant is force skip '.fs' file which is never cleaned by the pipeline and is created/removed 
                manually by the user. When 'fs' is cleaned by the user, the pipeline must be run with +nocache
                in order to notice the change. 
                Note that 'spawn' only schedules the tasks and the jobs are submitted to the farm by the 'wait' call.
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
                
=cut

sub spawn
{
    my ($self,$call,@args) = @_;

    if ( !$self->can($call) ) { $self->throw("No such method: [$call]\n"); }

    # Register the job
    my $job = {
        done_file => $args[0],
        call      => $call,
        args      => \@args,
    };
    push @{$$self{_checkpoints}}, $job;
}


# Consult cached database of finished files. Only if file is not present in the database,
#   a stat call is performed.
sub _is_marked_as_finished
{
    my ($self,$job) = @_;

    my $call      = $$job{call};
    my $done_file = $$job{done_file};
    my $wfile     = $$job{wait_file};
    my $is_dirty  = 0;
    my $is_done   = 0;

    # First time here, read the plain text "database" of finished files
    if ( !exists($$self{_jobs_db}{$done_file}) && !$$self{_nocache} )
    {
        if ( -e $wfile )
        {
            open(my $fh,'<',$wfile) or $self->throw("$wfile: $!");
            while (my $line=<$fh>)
            {
                chomp($line);
                if ( !($line=~/^([01sf])\t(\d+)\t(.*)$/) ) { $self->throw("Could not parse $wfile: $line\n"); }
                my $done = $1;
                my $id   = $2;  
                my $file = $3;
                $$self{_jobs_db}{$file}{finished} = $done;
                $$self{_jobs_db}{$file}{wfile}    = $wfile;
                $$self{_jobs_db}{$file}{call}     = $call;
                $$self{_jobs_db}{$file}{id}       = $id;
                $$self{_jobs_db}{$file}{dfile}    = $file;

                if ( !exists($$self{_max_ids}{$call}) or $$self{_max_ids}{$call}<$id ) { $$self{_max_ids}{$call} = $id; }
            }
            close($fh);
        }
        else { $is_dirty = 1; }
    }

    # No cache exists, init the job, set its ID and control file locations
    if ( !exists($$self{_jobs_db}{$done_file}) )
    {
        if ( !exists($$self{_jobs_db}{$done_file}{id}) ) { $$self{_jobs_db}{$done_file}{id} = ++$$self{_max_ids}{$call}; }
        $$self{_jobs_db}{$done_file}{call}  = $call;
        $$self{_jobs_db}{$done_file}{wfile} = $wfile;
        $$self{_jobs_db}{$done_file}{dfile} = $done_file;
        $$self{_jobs_db}{$done_file}{finished} = 0;
    }
    my $sfile  = "$wfile.$$self{_jobs_db}{$done_file}{id}.s";
    my $fsfile = "$wfile.$$self{_jobs_db}{$done_file}{id}.fs";

    # If skip file exists and +retries >0, the skip file will be deleted
    if ( $$self{_nretries}>0 && $$self{_jobs_db}{$done_file}{finished} eq 's' ) 
    { 
        $$self{_jobs_db}{$done_file}{finished} = 0; 
        $is_dirty = 1;
    }
    if ( $$self{_jobs_db}{$done_file}{finished} ) 
    { 
        if ( $$self{_jobs_db}{$done_file}{finished} eq 'f' ) { $self->debugln("Skipping the job upon request, a force skip file exists: $fsfile"); }
        if ( $$self{_jobs_db}{$done_file}{finished} eq 's' ) { $self->debugln("Skipping the job, a skip file exists: $sfile"); }
        return 2;
    }

    # Stat on non-existing files is cheap
    if ( $self->is_finished($done_file) ) 
    { 
        $is_dirty = 1;
        $is_done  = 1; 
    }
    # If the file needs to be skipped, then skip it by reporting that it's done even if it is not
    elsif ( -e $fsfile )
    {
        # The only way to clean a force skip file is to remove it manually
        $self->debugln("Skipping the job upon request, a force skip file exists: $fsfile");
        $is_done  = 'f';
        $is_dirty = 1;
    }
    elsif ( -e $sfile )
    {
        # This is currently the only way to clean the skip files: run with +retries set to positive value
        if ( $$self{_nretries}<0 ) 
        { 
            $self->debugln("Skipping the job, a skip file exists: $sfile");
            $is_done = 's'; 
        }
        else
        {
            $self->debugln("Cleaning skip file: $sfile");
            unlink($sfile);
        }
        $is_dirty = 1;
    }
    $$self{_jobs_db}{$done_file}{finished} = $is_done;
    $$self{_jobs_db_dirty} += $is_dirty;
    return $$self{_jobs_db}{$done_file}{finished};
}


sub _mark_as_finished
{
    my ($self) = @_;
    if ( !$$self{_jobs_db_dirty} ) { return; }

    my %wfiles = ();
    for my $job (values %{$$self{_jobs_db}})
    {
        push @{$wfiles{$$job{wfile}}}, $job;
    }
    for my $wfile (keys %wfiles)
    {
        $self->_mkdir($wfile);
        open(my $fh,'>',"$wfile.part") or $self->throw("$wfile.part: $!");
        my @clean_ids = ();
        my $all_done  = 1;
        for my $job (sort {$$a{id}<=>$$b{id}} @{$wfiles{$wfile}})
        {
            print $fh "$$job{finished}\t$$job{id}\t$$job{dfile}\n"; 
            if ( $$job{finished} ) { push @clean_ids, $$job{id}; }
            else { $all_done = 0; }
        }
        close($fh);
        if ( $$self{_js} ) 
        { 
            # Clean all jobs associated with this wfile
            $$self{_js}->clean_jobs($wfile,\@clean_ids,$all_done); 
        }
        rename("$wfile.part",$wfile) or $self->throw("rename $wfile.part $wfile: $!");
    }
}


sub _get_unfinished_jobs
{
    my ($self) = @_;

    my %calls = ();
    for my $job (@{$$self{_checkpoints}})
    {
        push @{$calls{$$job{call}}}, $job;
    }

    # Determine the base directory (common prefix) which holds the list of completed jobs
    my %wfiles = ();
    my $nprn_done = 0;
    my $nprn_pend = 0;
    for my $call (keys %calls)
    {
        my @list = @{$calls{$call}};
        my $min_len = length($list[0]{done_file});

        for (my $i=1; $i<@list; $i++)
        {
            my $len = length($list[$i]{done_file});
            if ( $len < $min_len ) { $min_len = $len; }

            my $j;
            for ($j=0; $j<$min_len; $j++)
            {
                if ( substr($list[0]{done_file},$j,1) ne substr($list[$i]{done_file},$j,1) ) { last; }
            }
            $min_len = $j;
        }
        my $dir = $self->_get_temp_dir(substr($list[0]{done_file},0,$min_len));
        for (my $i=0; $i<@{$calls{$call}}; $i++)
        {
            my $job = $calls{$call}[$i];
            $$job{wait_file} = "$dir/$call.w";
            my $ret;
            if ( ($ret=$self->_is_marked_as_finished($job)) )
            {
                if ( $nprn_done < 2 ) { $self->debugln("\to  $$job{done_file} .. " . ($ret ne '1' ? 'cached' : 'done')); $nprn_done++; }
                elsif ( $nprn_done < 3 ) { $self->debugln("\to  ...etc..."); $nprn_done++; }
                splice(@{$calls{$call}}, $i, 1);
                $i--;
            }
            else
            {
                my $id = $$self{_jobs_db}{$$job{done_file}}{id};
                if ( exists($wfiles{$$job{wait_file}}{$id}) ) 
                { 
                    $self->throw("The target file name is not unique: $$job{done_file}\n",Dumper($wfiles{$$job{wait_file}}{$id}{args},$$job{args})); 
                }
                if ( $nprn_pend < 2 ) { $self->debugln("\tx  $$job{done_file} .. unfinished"); $nprn_pend++; }
                elsif ( $nprn_pend < 3 ) { $self->debugln("\tx  ...etc..."); $nprn_pend++; }
                $wfiles{$$job{wait_file}}{$id} = $job;
            }
        }
        if ( !@{$calls{$call}} ) { delete($calls{$call}); }
    }
    $self->_mark_as_finished();
    return \%wfiles;
}


=head2 wait

    About : Checkpoint, submit the jobs scheduled by 'spawn' to the farm and wait for all tasks to finish. 
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

    # Initialize job scheduler (LSF, LSFCR, ...). This needs to be done here in
    # order for $js->clean_job() to work when $self->_get_unfinished_jobs() is called
    if ( !$$self{_run_locally} && !$$self{_js} )
    {
        my $farm = 'Runner' . $$self{_farm};
        eval 
        {
            require "$farm.pm";
            $$self{_js} = $farm->new();
        };
        if ( $@ ) { $self->throw("require $farm\n$@"); }
        if ( $$self{_maxjobs} ) { $$self{_js}->set_max_jobs($$self{_maxjobs}); }
    }

    # First check the files passed to wait() explicitly
    for my $file (@files)
    {
        if ( ! $self->is_finished($file) ) 
        { 
            $self->debugln("The file not finished: $file");
            exit; 
        }
    }
    if ( !exists($$self{_checkpoints}) or !scalar @{$$self{_checkpoints}} ) { return; }

    my (@caller) = caller(0);
    $self->debugln("Checking the status of ", scalar @{$$self{_checkpoints}}," job(s) .. $caller[1]:$caller[2]");
    my $jobs = $self->_get_unfinished_jobs();

    $$self{_checkpoints} = [];
    if ( !scalar keys %$jobs ) 
    { 
        $self->debugln("\t-> done");
        return; 
    }
    my $n = 0;
    for my $wfile (keys %$jobs) { $n += scalar keys %{$$jobs{$wfile}}; }
    $self->debugln("\t-> $n unfinished");

    # With '+local', the jobs will be run serially
    if ( $$self{_run_locally} ) 
    {
        for my $wfile (keys %$jobs)
        {
            $$self{_store} = $$jobs{$wfile};
            my $rfile = $self->freeze($wfile);
            for my $id (sort {$a<=>$b} keys %{$$jobs{$wfile}})
            {
                $self->_mkdir($$jobs{$wfile}{$id}{done_file});
                my $cmd = qq[$0 +run $rfile $id];
                $self->debugln("$$jobs{$wfile}{$id}{call}:\t$cmd");
                system($cmd);
                if ( $? ) { $self->throw("The command exited with a non-zero status $?\n"); }
            }
        }
        return;
    }

    # Spawn to farm
    my $js = $$self{_js};
    if ( $$self{_loop} )
    { 
        $$self{_farm_options}{wakeup_interval} = abs($$self{_loop}); 
        $js->set_limits(%{$$self{_farm_options}});
    }

    my $is_running = 0;

    # Each wfile corresponds to a single task group, each typically having multiple parallel jobs
    for my $wfile (keys %$jobs)     
    {
        my $prefix = $self->_get_temp_prefix($wfile);
        my $jobs_id_file = $prefix . '.jid';

        my @ids = sort { $a<=>$b } keys %{$$jobs{$wfile}};
        $self->debugln("\t.. ", scalar keys %$jobs > 1 ? scalar @ids."x\t$wfile" : "$wfile");

        my $tasks = $js->get_jobs($jobs_id_file,\@ids);
        my $is_wfile_running = 0;
        my $warned = 0;
        for (my $i=0; $i<@$tasks; $i++)
        {
            my $must_run = 1;
            my $done_file = $$jobs{$wfile}{$ids[$i]}{done_file};
            my $task = $$tasks[$i];

            if ( !defined $task )
            {
                $self->throw("\nCould not determine status of $i-th job: [$done_file] [$wfile] $ids[$i]\n");
            }

            # If the job is already running, don't check for errors - these could be from previous run anyway.
            if ( $js->job_running($task) ) 
            { 
                $must_run = 0;
                $is_running++;
                $is_wfile_running = 1;
                if ( $$self{_kill_jobs} ) { $js->kill_job($task); next; }
            }
        
            elsif ( $$self{_kill_jobs} ) { next; }

            # With very big arrays, it takes long time for is_marked_as_finished to complete
            #   and therefore we have to take farm's Done status seriously. Check again if
            #   the file has not appeared in the meantime, stat on non-existent files is fast anyway.
            elsif (  $js->job_done($task) )
            {
                if ( $$self{_nocache} && !$self->is_finished($done_file) ) { $must_run = 1; }
                else { $must_run = 0; }
            }

            # If the job has been already ran and failed, check if it failed repeatedly
            elsif ( $js->job_failed($task) ) 
            { 
                my $nfailures = $js->job_nfailures($task);
                if ( $nfailures > abs($$self{_nretries}) )
                {
                    if ( $$self{_nretries} < 0 )
                    {
                        my $sfile = "$wfile.$ids[$i].s";
                        $self->warn("\nThe job failed repeatedly (${nfailures}x) and +retries is negative, skipping: $wfile.$ids[$i].[eos]\n\n");
                        system("touch $sfile");
                        if ( $? ) { $self->throw("The command exited with a non-zero status $?: touch $sfile\n"); }
                        $must_run = 0;
                    }
                    else
                    {
                        my $msg = 
                            "The job failed repeatedly, ${nfailures}x: $wfile.$ids[$i].[eo]\n" .
                            "(Remove $jobs_id_file to clean the status, increase +retries or run with negative value of +retries to skip this task.)\n";

                        $self->_send_email('failed', "The runner failed repeatedly\n", $$self{_about}, "\n", $msg);
                        $self->throw($msg);
                    }
                }
                elsif ( !$warned )
                {
                    $self->warn("\nRunning again, the job failed or needs restart: $wfile.$ids[$i].[eo]\n\n");
                    $warned = 1;
                }

                # Increase memory limits if necessary: by a set minimum or by a percentage, which ever is greater
                my %limits = $js->past_limits($ids[$i],$wfile);
                if ( exists($limits{MEMLIMIT}) )
                { 
                    my $mem = $limits{memory}*1.3 > $limits{memory}+1_000 ? $limits{memory}*1.3 : $limits{memory}+1_000;
                    $self->inc_limits(memory=>$mem); 
                }
            }

            if ( $must_run )
            { 
                if ( $$self{_maxjobs} && $$self{_maxjobs}<$is_running ) { last; }
                $self->_update_limits($wfile, $ids[$i]);
                next; 
            }

            splice(@ids, $i, 1);
            splice(@$tasks, $i, 1);
            $i--;
        }
        if ( $$self{_kill_jobs} ) { next; }
        if ( !@ids ) 
        { 
            if ( !$is_wfile_running ) { unlink($jobs_id_file); }
            next; 
        }
        if ( $$self{_maxjobs} )
        {
            if ( $$self{_maxjobs} <= $is_running ) { last; }
            if ( $$self{_maxjobs} < $is_running + @ids ) { splice(@ids, $$self{_maxjobs} - $is_running - @ids); }
        }
        $is_running += scalar @ids;

        $$self{_store} = $$jobs{$wfile};
        my $rfile = $self->freeze($wfile);
        my $cmd = qq[$0 +run $rfile {JOB_INDEX}];
        $self->debugln("$wfile:\t$cmd");

        for my $id (@ids)
        {
            $self->_mkdir($$jobs{$wfile}{$id}{done_file});
        }

        my $ok;
        eval 
        {
            $js->set_limits(%{$$self{_farm_options}});
            $js->run_jobs($jobs_id_file,$prefix,$cmd,\@ids);
            $ok = 1;
        };
        if ( !$ok )
        {
            $self->throw($@);
        }
    }
    if ( $$self{_kill_jobs} ) { $self->all_done; }
    if ( $is_running ) { exit; }
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
    $self->_send_email('done', "The runner has finished, all done!\n", $$self{_about});
    exit $$self{_status_codes}{DONE};
}

sub _send_email
{
    my ($self,$status, @msg) = @_;
    if ( !exists($$self{_mail}) ) { return; }
    open(my $mh,"| mail -s 'Runner report: $status' $$self{_mail}");
    print $mh join('',@msg) . "\n";
    close($mh);
}

=head2 clean

    About : Clean all system files (in .jobs directories)
    Usage : $self->clean($dir);
    Args  : <@dirs>
                Directories to recursively clean from all .jobs subdirs leaving a single tarball instead

=cut

sub clean
{
    my ($self,@dirs) = @_;

    if ( !@dirs ) { return; }

    # Create the tarball, existing file will not be overwritten
    my $tarball = $dirs[0] . '/cleaned-job-outputs.tgz';
    if ( !-e $tarball )
    {
        my $dirs = join(' ',@dirs);
        my $cmd = "find $dirs -name .jobs | tar -T - -czf $tarball";
        $self->debugln($cmd);
        system($cmd);
    }
    for my $dir (@dirs)
    {
        my $cmd = "find $dir -name .jobs | xargs rm -rf";
        $self->debugln($cmd);
        system($cmd);
    }
}

=head2 is_finished

    About : Check if the file is finished.
    Usage : $self->is_finished('some/file');
    Args  : <file list>
                The name of the file to check the existence of
                
=cut

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

sub _is_storable
{
    my ($self,$file) = @_;
    for (my $i=0; $i<3; $i++)
    {
        if ( $i>0 ) { sleep 5; }
        my @out = `file -b $file`;
        if ( $? ) { next; }
        if ( $out[0] =~ /^perl Storable/ ) { return 1; }
        return 0;
    }
    $self->throw("Cannot [file -b $file]: $!");
}

sub _revive_array
{
    my ($self,$freeze_file,$job_id) = @_;
    
    my $back;
    for (my $i=0; $i<3; $i++)
    {
        if ( $i>0 ) { sleep 5; }
        eval { $back = retrieve($freeze_file); };
        if ( !$@ ) { last; }
    }
    if ( $@ ) { $self->throw("retrieve() threw an error: $freeze_file\n$@\n"); }
    if ( $$self{clean} ) { unlink($freeze_file); }

    $self = $back;

    if ( !exists($$self{_store}{$job_id}) ) { $self->throw("No such job $job_id in $freeze_file\n"); }
    my $job = $$self{_store}{$job_id};
    $$self{_revived_file} = $freeze_file;
    $$self{_revived_job}  = $job_id;

    my $code = $self->can($$job{call});
    &$code($self,@{$$job{args}});
}


# Run the freezed object created by spawn. The $config_file argument can supply custom
#	settings for overriding and debugging settings in the freezed file.
sub _revive
{
    my ($self,$freeze_file,$config_file) = @_;

    my $back;
    eval { $back = retrieve($freeze_file); };
    if ( $@ ) { $self->throw("retrieve() threw an error: $freeze_file\n$@\n"); }
    if ( $$self{clean} ) { unlink($freeze_file); }

    $self = $back;

	if ( defined $config_file )
	{
		my %x = do "$config_file";
        if ( $@ ) { $self->throw("do $config_file: $@\n"); }
		while (my ($key,$value) = each %x) { $$self{$key} = $value; }
	}
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
    if ( !($fname=~m{([^/]+)$}) ) { $self->throw("FIXME: could not parse [$fname]\n"); }
    my $file = $1;
    my $dir  = $self->_get_temp_dir($fname);
    if ( ! -d $dir ) 
    { 
        `mkdir -p $dir`; 
        if ( $? ) { $self->throw("Cannot create directory [$dir]: $!"); }
    }
    return "$dir/$file";
}

sub _get_temp_dir
{
    my ($self,$fname) = @_;
    my @items = split(m{/}, $fname);
    my $len = length($fname);
    if ( substr($fname,$len-1,1) ne '/' ) { splice(@items, -1); }
    if ( $items[-1] eq '.jobs' ) { splice(@items,-1); }
    return join('/',@items,'.jobs');
}

sub _mkdir
{
    my ($self,$fname) = @_;
    $fname =~ s{[^/]+$}{};
    if ( !-e $fname ) { `mkdir -p $fname`; }
    return $fname;
}


=head2 cmd

    About : Executes a command via bash in the -o pipefail mode. 
    Args  : <string>
                The command to be executed
            <hash>
                Optional arguments: 
                - verbose           .. print command to STDERR before executing [0]
                - require_status    .. throw if exit status is different [0]

=cut

sub cmd
{
    my ($self,$cmd,%args) = @_;

    if ( $$self{verbose} or $args{verbose} ) { print STDERR $cmd,"\n"; }

    # Why not to use backticks? Perl calls /bin/sh, which is often bash. To get the correct
    #   status of failing pipes, it must be called with the pipefail option.

    my $kid_io;
    my $pid = open($kid_io, "-|");
    if ( !defined $pid ) { $self->throw("Cannot fork: $!"); }

    my @out;
    if ($pid) 
    {
        # parent
        @out = <$kid_io>;
        close($kid_io);
    } 
    else 
    {      
        # child
        exec('/bin/bash', '-o','pipefail','-c', $cmd) or $self->throw("Failed to run the command [/bin/sh -o pipefail -c $cmd]: $!");
    }

    my $exit_status = $?;
    my $status = exists($args{require_status}) ? $args{require_status} : 0;
    if ( $status ne $exit_status ) 
    {
        my $msg;
        if ( $? & 0xff )
        {
            $msg = "The command died with signal ".($? & 0xff);
        }
        else
        {
            $msg = "The command exited with status ".($? >> 8)." (expected $status)";
        }
        $msg .= ":\n\t$cmd\n\n";
        if ( @out ) {  $msg .= join('',@out,"\n\n"); }
        $self->throw($msg); 
    }
    return @out;
}

=head2 throw

    About : Throws an error.
    Args  : <array>
                The message to be printed. If no message is given, the usage will be printed instead.

=cut

sub throw
{
    my ($self,@msg) = @_;
    $! = $$self{_status_codes}{ERROR};
    if ( scalar @msg ) { confess "\n[". scalar gmtime() ."]\n", @msg; }
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

sub debug
{
    my ($self,@msg) = @_;
    if ( !$$self{_verbose} ) { return; }
    print STDERR @msg;
}

=head1 AUTHORS

petr.danecek@sanger

=cut

1;
