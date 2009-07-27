=head1 NAME

VertRes::Wrapper::WrapperI - interface for wrapper modules

=head1 SYNOPSIS

package VertRes::Wrapper::MyWrapper;
use base qw(VertRes::Wrapper::WrapperI);

sub new {
    my ($class, @args) = @_;

    my $self = $class->SUPER::new(@args,
                                  exe      => 'myexe',
                                  params   => [qw(param1 param2)],
                                  switches => [qw(switch1 switch2)]);

    return $self;
}

1;

package main;
use VertRes::Wrapper::MyWrapper;

my $wrapper = VertRes::Wrapper::MyWrapper->new(param2 => 17, switch1 => 1);
$wrapper->run('input_file', 'output_file');

# that last call to run() ends up calling:
# myexe -param2 17 -switch1 input_file output_file

=head1 DESCRIPTION

Make a wrapper (a module that runs an external executable) that can get/set all
the programs args easily.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Wrapper::WrapperI;

use strict;
use warnings;
use Cwd qw(abs_path);
use IPC::Open2;

use base qw(VertRes::Base);

our %allowed_run_methods = (bsub => 1, open => 1, open_to => 1, open_2 => 1, system => 1);

=head2 new

 Title   : new
 Usage   : my $self = $class->SUPER::new(@args, params => [qw(param1)]);
 Function: Instantiates your object, including setting all your params and
           switches.
 Returns : $self hash-ref blessed into your class
 Args    : quiet      => boolean
           exe        => string
           params     => [params]
           switches   => [switches]
           run_method => bsub|system|open

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    # set our own args
    $self->_set_from_args(\@args,
                          methods => [qw(quiet exe params switches run_method)]);
    
    # now set the params and switches of the exe our child class is wrapping,
    # creating the get/setter methods if necessary
    $self->_set_params_and_switches_from_args(@args);
    
    $self->_set_run_status(-2);
    
    return $self;
}

sub _set_params_and_switches_from_args {
    my $self = shift;
    
    my @methods;
    foreach my $hash_ref ($self->params(), $self->switches()) {
        $hash_ref || next;
        push(@methods, keys %{$hash_ref});
    }
    if (@methods) {
        my @args = @_;
        # first reset to undefs and 0s
        my $hash_ref = $self->params();
        foreach my $param (keys %{$hash_ref || {}}) {
            unshift(@args, ($param => undef));
        }
        $hash_ref = $self->switches();
        foreach my $switch (keys %{$hash_ref || {}}) {
            unshift(@args, ($switch => 0));
        }
        
        # now set the user's args
        $self->_set_from_args(\@args,
                              methods => \@methods,
                              create => 1,
                              case_sensitive => 1);
    }
}

=head2 exe

 Title   : exe
 Usage   : my $exe = $obj->exe();
 Function: Get/set the executable name.
 Returns : executable name string 
 Args    : none to get, string to set

=cut

sub exe {
    my $self = shift;
    if (@_) {
        $self->{_exe} = shift;
    }
    return $self->{_exe};
}

=head2 version

 Title   : version
 Usage   : my $version = $obj->version()
 Function: Returns the program version (if available).
 Returns : string representing version of the program 
 Args    : n/a

=cut

sub version {
    my $self = shift;
    return;
}

=head2 quiet

 Title   : quiet
 Usage   : $obj->quiet(1);
           if ($factory->quiet()) { ... }
 Function: Get/set the quiet state. Can be used by wrappers to control if
           program output is printed to the console or not.
 Returns : boolean
 Args    : none to get, boolean to set

=cut

sub quiet {
    my $self = shift;
    if (@_) { $self->{quiet} = shift }
    return $self->{quiet} || 0;
}

=head2 run

 Title   : run
 Usage   : $wrapper->run();
 Function: Support generic running of exe() with executable parameters set via
           _get_params_string(). Prior to doing anything calls
           $self->_pre_run(@_), giving wrapper module authors the opportunity to
           do something special first, like call _set_params_string() with
           custom args, minimising the need to override this method. _pre_run
           should return args suitable for passing to the program without any
           parameters. It then chooses how to run the exe based on
           run_method(). By default it will run via bsub.
 Returns : the value of $self->_post_run(@output_of_run), which depends on the
           run_method()
 Args    : passed on to the exe() call after the _get_params_string() params

=cut

sub run {
    my $self = shift;
    my @extra_args = $self->_pre_run(@_);
    my $params = $self->_get_params_string();
    
    my $exe = $self->exe() || $self->throw("exe is unknown, can't run!");
    
    my $run_method = $self->run_method();
    $run_method = '_'.$run_method.'_run';
    $self->_set_run_status(0);
    my @result = $self->_post_run($self->$run_method($exe, $params, @extra_args));
    
    return wantarray ? @result : $result[0];
}

=head2 run_method

 Title   : via_bsub
 Usage   : $wrapper->run_method('bsub');
 Function: Wrapper module authors can choose how to run the executable if they
           don't want to override run(). End-users can also override.
 Returns : string
 Args    : none to get, string to set:
           system (run in a system call - default)
           bsub (run via bsub using the options set in bsub_options())
           open (run in an open() call that returns a filehandle you must
                 read through - useful for piping straight through to a parser)
           open_to (run in an open() call where one of the file args is a
                    filehandle that will be piped into the executable)
           open_2 (run in an IPC::Open2 call where input is a supplied filehanle
                   like for open_to and you are returned a filehandle like open)

=cut

sub run_method {
    my ($self, $method) = @_;
    
    if ($method && exists $allowed_run_methods{$method}) {
        $self->{_run_method} = $method;
    }
    
    return $self->{_run_method} || 'system';
}

sub _pre_run {
    my $self = shift;
    return @_;
}

sub _post_run {
    my $self = shift;
    $self->_set_run_status(1);
    $self->check_output_files;
    return @_;
}

sub _bsub_run {
    my ($self, $exe, $params, @extra_args) = @_;
    my $bo = $self->bsub_options();
    
    my $optional = '';
    $optional .= ' -R '.$bo->{R} if $bo->{R};
    $optional .= ' -M '.$bo->{M} if $bo->{M};
    
    my $command = "bsub -J $bo->{J} -e $bo->{e} -o $bo->{o} -q $bo->{q}$optional '$exe$params @extra_args'";
    $self->debug("will run command '$command'");
    
    system($command) && $self->throw("$exe call ($command) crashed: $? | $!");
    
    return;
}

sub _open_run {
    my ($self, $exe, $params, @extra_args) = @_;
    
    my $redirect = $self->quiet ? ' 2> /dev/null ' : '';
    my $command = $exe.$params." @extra_args".$redirect;
    $self->debug("will run command '$command'");
    
    open(my $pipe, "$command |") || $self->throw("$exe call ($command) failed to start: $? | $!");
    
    return $pipe;
}

sub _open_to_run {
    my ($self, $exe, $params, @extra_args) = @_;
    
    my $fh;
    foreach (@extra_args) {
        if (ref($_) && ref($_) eq 'GLOB') {
            $fh = $_;
            $_ = '-';
        }
    }
    
    my $redirect = $self->quiet ? ' 2> /dev/null ' : '';
    my $command = $exe.$params." @extra_args".$redirect;
    $self->debug("will run command '$command'");
    
    open(my $pipe, "| $command") || $self->throw("$exe call ($command) failed to start: $? | $!");
    
    while (<$fh>) {
        print $pipe $_;
    }
    close($fh);
    close($pipe);
    
    return;
}

sub _open_2_run {
    my ($self, $exe, $params, @extra_args) = @_;
    
    my $in_fh;
    foreach (@extra_args) {
        if (ref($_) && ref($_) eq 'GLOB') {
            $in_fh = $_;
            $_ = '-';
        }
    }
    
    my $redirect = $self->quiet ? ' 2> /dev/null ' : '';
    my $command = $exe.$params." @extra_args".$redirect;
    $self->debug("will run command '$command'");
    
    my ($out_fh);
    my $pid = open2($out_fh, $in_fh, $command);
    push(@{$self->{_open_2_pids}}, $pid);
    print "got pid $pid\n";
    
    return $out_fh;
}

sub _system_run {
    my ($self, $exe, $params, @extra_args) = @_;
    
    my $redirect = $self->quiet ? ' 2> /dev/null ' : '';
    my $command = $exe.$params." @extra_args".$redirect;
    $self->debug("will run command '$command'");
    
    system($command) && $self->throw("$exe call ($command) crashed: $? | $!");
    
    return;
}

=head2 run_status

 Title   : run_status
 Usage   : my $status = $wrapper->run_status();
 Function: Discover what happened with your run.
 Returns : int, with meanings as follows:
           -2 = not yet attempted to run
           -1 = ran and failed, but repeating the run might result in success
           0  = ran and failed with a permanent error
           1  = ran seemingly to success (the best you can get when running via
                                          bsub or open)
           2  = ran definitely to success
 Args    : n/a

=cut

sub run_status {
    my $self = shift;
    return $self->{_run_status};
}

sub _set_run_status {
    my ($self, $status) = @_;
    $self->{_run_status} = $status;
}

=head2 params

 Title   : params
 Usage   : $wrapper->params([qw(param1 param2)]);
 Function: Get/set all the paramaters that are understood by the executable
           being wrapped.
 Returns : hashref of params (where keys are method names for getting/setting
           the param values and values are how those params should be provided
           to the executable - usually the same as the key)
 Args    : array ref of method names to call, or hash ref where keys are method
           names and values are how those names should be output in the params
           string

=cut

sub params {
    my ($self, $ref) = @_;
    
    if ($ref) {
        $self->{_params} = ref($ref) eq 'HASH' ? $ref : {map { $_ => $_ } @{$ref}};
        delete $self->{_params_string};
    }
    
    return $self->{_params};
}

=head2 switches

 Title   : switches
 Usage   : $wrapper->switches([qw(switch1 switch2)]);
 Function: Get/set all the switches that are understood by the executable
           being wrapped.
 Returns : hashref of switches (where keys are method names for getting/setting
           the param values and values are how those params should be provided
           to the executable - usually the same as the key)
 Args    : array ref of method names to call, or hash ref where keys are method
           names and values are how those names should be output in the params
           string

=cut

sub switches {
    my ($self, $ref) = @_;
    
    if ($ref) {
        $self->{_switches} = ref($ref) eq 'HASH' ? $ref : {map { $_ => $_ } @{$ref}};
        delete $self->{_params_string};
    }
    
    return $self->{_switches};
}

=head2  _set_params_string()

 Title   : _set_params_string
 Usage   : $self->_set_params_string(dash => 1);
 Function: For internal use by wrapper modules to build parameter strings
           suitable for sending to the program being wrapped. For each method
           name defined in params() and switches(), calls the method and adds
           the method name (as modified by optional things) along with its value
           (unless a switch) to the parameter string.
 Example : $self->params([qw(window evalue_cutoff)]);
           $self->switches([qw(simple large all)]);
           $self->_set_params_string(-double_dash => 1,
                                     -underscore_to_dash => 1);
           my $params = $self->_get_params_string();
           If window() and simple() had not been previously called, but
           evalue_cutoff(0.5), large(1) and all(0) had been called, $params
           would be ' --evalue-cutoff 0.5 --large'
 Returns : n/a
 Args    : join => string      # define how parameters and their values are
                                 joined, default ' '. (eg. could be '=' for
                                 param=value)
           lc => boolean       # lc() method names prior to output in string
           dash => boolean     # prefix all method names with a single dash
           double_dash => bool # prefix all method names with a double dash
           mixed_dash => bool  # prefix single-character method names with a
                               # single dash, and multi-character method names
                               # with a double-dash
           underscore_to_dash => boolean # convert all underscores in method
                                           names to dashes

=cut

sub _set_params_string {
    my ($self, %args) = @_;
    
    my ($join, $lc, $d, $dd, $md, $utd) =
       ($args{join}, $args{lc}, $args{dash}, $args{double_dash},
        $args{misxed_dash}, $args{underscore_to_dash});
    $self->throw("-dash, -double_dash and -mixed_dash are mutually exclusive") if (defined($d) + defined($dd) + defined($md) > 1);
    $join ||= ' ';
    
    my $params = $self->params();
    my $switches = $self->switches();
    
    my $param_string = '';
    for my $hash_ref ($params, $switches) {
        while (my ($method, $method_out) = each %{$hash_ref}) {
            my $value = $self->$method();
            next unless (defined $value);
            next if (exists $switches->{$method} && ! $value);
            
            $method_out = lc($method_out) if $lc;
            my $method_length = length($method_out) if $md;
            $method_out = '-'.$method_out if ($d || ($md && ($method_length == 1)));
            $method_out = '--'.$method_out if ($dd || ($md && ($method_length > 1)));
            $method_out =~ s/_/-/g if $utd;
            
            # quote values that contain spaces
            if (exists $params->{$method} && $value =~ /^[^'"\s]+\s+[^'"\s]+$/) {
                $value = '"'.$value.'"';
            }
            
            $param_string .= ' '.$method_out.(exists $switches->{$method} ? '' : $join.$value);
        }
    }
    
    $self->{_params_string} = $param_string;
}

=head2 _get_params_string

 Title   : _get_params_string
 Usage   : my $params = $wrapper->_get_params_string();
 Function: Get the parameter string for supplying to the executable. If
           _set_params_string() had not previously been called, calls
           _set_params_string(dash => 1) as the default first.
 Returns : params string
 Args    : n/a

=cut

sub _get_params_string {
    my $self = shift;
    unless (defined $self->{_params_string}) {
        $self->_set_params_string(dash => 1);
    }
    return $self->{_params_string};
}

=head2 bsub_options

 Title   : bsub_options
 Usage   : my $bus_options_hash_ref = $wrapper->bsub_options();
 Function: Get/set bsub options.
 Returns : hash ref, keys as below
 Args    : none to get, key => value pairs to set:
           J => jobname (default is something unique and sensible)
           o => output file (bsub output file, default jobname.o in current dir)
           e => error file (bsub error output file, default jobname.e ")
           q => queue (default normal)
           M => memory size in bytes (default unset)
           R => memory selection statement (default unset)

=cut

sub bsub_options {
    my $self = shift;
    
    unless (defined $self->{_bsub_opts}) {
        my $class = ref($self);
        $class =~ s/VertRes::Wrapper:://;
        my $user = $ENV{USER} || 'vr';
        my $J = "vr_wrap_$class.$user.$$";
        $self->{_bsub_opts} = {J => $J,
                               o => $J.'.o',
                               e => $J.'.e',
                               q => 'normal'};
    }
    
    if (@_) {
        my %user_args = @_;
        foreach my $arg (qw(J o e q M R)) {
            $self->{_bsub_opts}->{$arg} = $user_args{$arg} || next;
        }
    }
    
    return $self->{_bsub_opts};
}

=head2 register_output_file_to_check

 Title   : register_output_file_to_check
 Usage   : $wrapper->register_output_file_to_check('out_file');
 Function: Set a filepath that you expect to be written to by the executable
           you are wrapping; these files are, by default, checked for size
           with check_output_files() in _post_run().
 Returns : n/a
 Args    : filepath

=cut

sub register_output_file_to_check {
    my ($self, $out_file) = @_;
    if ($out_file) {
        push(@{$self->{_out_files_to_check}}, $out_file);
    }
}

=head2 check_output_files

 Title   : check_output_files
 Usage   : $wrapper->check_output_files();
 Function: Checks to see if the files registered by
           register_output_file_to_check() exist and have some size: a basic
           kind of check to see if your executable has produced some output.
           If any output file is empty or non-existant, the run_status will
           become -1.
 Returns : n/a (all files to check will also be unregistered afterwards)
 Args    : n/a

=cut

sub check_output_files {
    my $self = shift;
    if (defined $self->{_out_files_to_check}) {
        foreach my $file (@{$self->{_out_files_to_check}}) {
            unless (-s $file) {
                $self->_set_run_status(-1);
            }
        }
    }
    delete $self->{_out_files_to_check};
}

# is file2 older than file1?
sub _is_older {
    my ($self, $file1, $file2) = @_;
    -e $file1 || return 0;
    -e $file2 || return 0;
    
    # get the stats of the actual files, not potentially their symlinks
    my @stat = stat(abs_path($file1));
    my $first_mtime = $stat[9];
    @stat = stat(abs_path($file2));
    my $second_mtime = $stat[9];
    
    return $second_mtime < $first_mtime;
}

sub DESTROY {
    my $self = shift;
    
    foreach my $pid (@{$self->{_open_2_pids} || []}) {
        waitpid $pid, 0;
    }
}

1;
