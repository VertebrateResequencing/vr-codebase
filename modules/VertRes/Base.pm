=head1 NAME

VertRes::Base - base module that all VertRes modules inherit from

=head1 SYNOPSIS

package VertRes::MyNewModule;

use base qw(VertRes::Base);

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    $self->_set_from_args(\@args, methods => [qw(mymethod1 mymethod2)]);

    # MyNewModule-specific stuff goes here

    return $self;
}

sub example_sub {
    my $self = shift;

    # throw an exception
    $self->throw("This is an exception") if $bad_thing_happened;

    # catch and handle a possible exception
    eval {
        my $other = VertRes::othermodule->new()
        $other->dangerous_method();
    };
    if ($@) {
        # got an exception, do something about it...
    }
    else {
        # didn't get an exception, do something else...
    }

    # warn about something
    $self->warn("This is a warning") if $something_not_quite_right;

    # print a debugging message for those who want to read it
    $self->debug("This is a debugging message");

    # change the verbosity of another object
    my $other = VertRes::othermodule->new(verbose => 1);
    $other->verbose(2);

    # turn on logging of the above kinds of messages
    $self->write_logs(1);
    $self->log_file('/my/special/logfile/location');

    # only log a message to the log file; don't print to STDERR
    $self->log('my log message');

    # if you do something that will require special consideration when
    # destroying ourselves, write a method to handle that and register it
    $self->register_for_cleanup('my_cleanup_method');
}

=head1 DESCRIPTION

Provides fundamental needs for all modules: error and debug handling, logging,
object creation and destruction.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Base;

use strict;
use warnings;
use Carp qw(carp cluck confess);
use Time::Format;
use Cwd qw(getcwd);
use File::Spec;
use IO::Capture::Stderr;

our $GLOBAL_VERBOSITY;
our $VERSION = '0.1';
$VERSION = eval $VERSION;

=head2 new

 Title   : new
 Usage   : my $self = $class->SUPER::new(@args);
 Function: Generic new method, also sets up VertRes essentials. Auto-sets args
           and calls methods as per _set_from_args().
 Returns : $self hash-ref blessed into your class
 Args    : verbose => int
           log_file => path (default ~/.vertres.log)
           write_logs => boolean

=cut

sub new {
    my $class = shift;
    my $self = {};
    bless $self, ref($class) || $class;
    
    # setup default log settings
    $self->{_log_file} = File::Spec->catfile($ENV{HOME}, '.vertres.log');
    $self->{_write_logs} = 0;
    
    $self->_set_from_args(\@_, methods => [qw(verbose log_file write_logs)]);
    
    return $self;
}

=head2 _set_from_args

 Usage     : $self->_set_from_args(\%args, methods => \@methods);
 Purpose   : Takes a hash of user-supplied args whose keys match method names,
             and calls the method supplying it the corresponding value. If a
             key doesn't match a method name, and neither force and create are
             in effect, the key => value pair is set in the $self hash.
 Example   : $self->_set_from_args(\%args, methods => [qw(sequence id desc)]);
             Where %args = (sequence    => $s,
	                    description => $d,
	                    ID          => $i);

             the above _set_from_args() calls the following methods:
             $self->sequence($s);
             $self->id($i);
             ( $self->description($i) is not called because 'description' wasn't
               one of the given methods; instead $self->{description} = $d is
               set )
 Argument  : \%args | \@args : a hash ref or associative array ref of arguments
                               where keys are any-case strings corresponding to
                               method names, and values are the values the
                               method should be supplied. If keys contain 
                               internal hyphens (eg. to separate multi-word 
                               args) they are converted to underscores, since
                               method names cannot contain dashes.
             methods => []   : (optional) only call methods with names in this
                               array ref. Can instead supply a hash ref where
                               keys are method names (of real existing methods
                               unless -create is in effect) and values are array
                               refs of synonyms to allow access to the method
                               using synonyms. If there is only one synonym it
                               can be supplied as a string instead of a single-
                               element array ref.
             force => bool   : (optional, default 0) call methods that don't
                               seem to exist, ie. let AUTOLOAD handle them.
             create => bool  : (optional, default 0) when a method doesn't
                               exist, create it as a simple getter/setter
                               (combined with -methods it would create all the
                               supplied methods that didn't exist, even if not
                               mentioned in the supplied %args).
             code => '' | {} : (optional) when creating methods use the supplied
                               code (a string which will be evaulated as a sub).
                               The default code is a simple get/setter.
                               Alternatively you can supply a hash ref where
                               the keys are method names and the values are
                               code strings. The variable '$method' will be
                               available at evaluation time, so can be used in
                               your code strings. Beware that the strict pragma
                               will be in effect.
             case_sensitive => bool : require case sensitivity on the part of
                                      user (ie. a() and A() are two different
                                      methods and the user must be careful
                                      which they use). Default 0 (insenstive).

 Comments  : The \%args argument will usually be the args received during new()
             from the user. The user is allowed to get the case wrong, include
             0 or more than one hyphens as a prefix, and to include hyphens as
             multi-word arg separators: '--an-arg' => 1, -an_arg => 1 and
             An_Arg => 1 are all equivalent, calling an_arg(1). However, in
             documentation users should only be told to use the standard form
             an_arg to avoid confusion. A possible exception to this is a
             wrapper module where '--an-arg' is what the user is used to
             supplying to the program being wrapped.

             Another issue with wrapper modules is that there may be an
             argument that has meaning both to VertRes and to the program, eg.
             verbose. The recommended way of dealing with this is to leave
             verbose to set the VertRes verbosity whilst requesting users use
             an invented program_verbose (or similar) to set the program
             verbosity. This can be resolved back with
             VertRes::Wrapper::WrapperI's _setparams() method and code along
             the lines of:
             my %methods = map { $_ => $_ } @LIST_OF_ALL_ALLOWED_PROGRAM_ARGS;
             delete $methods{'verbose'};
             $methods{'program_verbose'} = 'verbose';
             my $param_string = $self->_setparams(-methods => \%methods);
             system("$exe $param_string");

=cut

sub _set_from_args {
    my ($self, $args, %own_args) = @_;
    $self->throw("a hash/array ref of arguments must be supplied") unless ref($args);
    
    my ($methods, $force, $create, $code, $case);
    if (keys %own_args) {
        $methods = $own_args{methods};
        $force = $own_args{force};
        $create = $own_args{create};
        $code = $own_args{code};
        $case = $own_args{case_sensitive};
    }
    my $default_code = 'my $self = shift;
                        if (@_) { $self->{\'_\'.$method} = shift }
                        return $self->{\'_\'.$method};';
    
    my %method_names = ();
    my %syns = ();
    if ($methods) {
        my @names;
        if (ref($methods) eq 'HASH') {
            @names = keys %{$methods};
            %syns = %{$methods};
        }
        else {
            @names = @{$methods};
            %syns = map { $_ => $_ } @names;
        }
        %method_names = map { $case ? $_ : lc($_) => $_ } @names;
    }
    
    # deal with hyphens
    my %orig_args = ref($args) eq 'HASH' ? %{$args} : @{$args};
    my %args;
    while (my ($method, $value) = each %orig_args) {
        $method =~ s/^-+//;
        $method =~ s/-/_/g;
        $args{$method} = $value;
    }
    
    # create non-existing methods on request
    if ($create) {
        unless ($methods) {
            %syns = map { $_ => $case ? $_ : lc($_) } keys %args;
        }
        
        foreach my $method (keys %syns) {
            $self->can($method) && next;
            
            my $string = $code || $default_code;
            if (ref($code) && ref($code) eq 'HASH') {
                $string = $code->{$method} || $default_code;
            }
            
            my $sub = eval "sub { $string }";
            $self->throw("Compilation error for $method : $@") if $@;
            
            no strict 'refs';
            *{ref($self).'::'.$method} = $sub;
        }
    }
    
    # create synonyms of existing methods
    while (my ($method, $syn_ref) = each %syns) {
        my $method_ref = $self->can($method) || next;
        
        foreach my $syn (@{ ref($syn_ref) ? $syn_ref : [$syn_ref] }) {
            next if $syn eq $method;
            $method_names{$case ? $syn : lc($syn)} = $syn;
            next if $self->can($syn);
            no strict 'refs';
            *{ref($self).'::'.$syn} = $method_ref;
        }
    }
    
    # set values for methods
    my %set_as_method = %{$self->{_methods_already_set} || {}};
    while (my ($method, $value) = each %args) {
        $method = $method_names{$case ? $method : lc($method)} || ($methods ? next : $method);
        $self->can($method) || next unless $force;
        $self->$method($value);
        $set_as_method{$method} = 1;
    }
    $self->{_methods_already_set} = \%set_as_method;
    
    # set values for vars
    while (my ($var, $value) = each %args) {
        next if $set_as_method{$var};
        $self->{$var} = $value;
    }
}

=head2 verbose

 Title   : verbose
 Usage   : $self->verbose(1); # per-object
           VertRes::Base::verbose(1); # global
 Function: Sets verbose level for how ->warn() behaves
           -1 = no warning
            0 = standard, small warning
	    0.5 = tiny warning with no line numbers
            1 = warning with stack trace
            2 = warning becomes throw
 Returns : The current verbosity setting (integer between -1 to 2)
 Args    : -1,0,1 or 2

=cut

sub verbose {
    my ($self, $got_input);
    if (@_) {
        # (we can pick up undef as an arg when called as a class method
        #  this way)
        $self = shift if ref $_[0];
        $got_input = 1 if (! $self || @_);
    }
    
    if ($got_input) {
        my $value = shift;
        if (ref $self) {
            $self->{_verbose} = $value;
        }
        else {
            # set globally if not an instantiation
            $GLOBAL_VERBOSITY = $value;
        }
    }
    
    if (defined $GLOBAL_VERBOSITY) {
        return $GLOBAL_VERBOSITY;
    }
    
    return $self->{_verbose} || 0;
}

=head2 warn

 Title   : warn
 Usage   : $obj->warn("warning message");
 Function: Warns about something; strength of warning determined by ->verbose().
           If logging is turned on will also output warning to log file.
 Returns : n/a
 Args    : A string giving a warning message

=cut

sub warn {
    my ($self, $message) = @_;
    
    my $verbose = $self->verbose();
    return if $verbose <= -1;
    
    if ($verbose >= 2) {
        $self->throw($message);
    }
    elsif ($verbose == 1) {
        cluck($message);
    }
    elsif ($verbose == 0.5) {
	CORE::warn($message."\n");
    }
    else {
        carp($message);
    }
    
    $self->log($message);
}

=head2 debug

 Title   : debug
 Usage   : $obj->debug("This is debugging output");
 Function: Warns a debugging message when verbose is > 0.
           If logging is turned on will also output message to log file.
 Returns : n/a
 Args    : Message string to debug about

=cut

sub debug {
    my ($self, $message) = @_;
    
    if ($self->verbose > 0) {
        $self->log($message);
        $message .= "\n" unless $message =~ /\n$/;
        CORE::warn $message;
    }
}

=head2 throw

 Title   : throw
 Usage   : $obj->throw("throwing exception message");
 Function: Throws an exception, which, if not caught with an eval or
           a try block will provide a nice stack trace to STDERR
           with the message. (Uses Carp's confess().) Automatically includes
           the date and current working directory.
           If logging is turned on will also output throw message to log file.
 Returns : n/a
 Args    : A string giving a descriptive error message

=cut

sub throw {
    my ($self, $message) = @_;
    
    $message ||= '[no message]';
    $self->log($message, 2);
    
    my $cwd = getcwd();
    $message = "FATAL ERROR on $time{'yyyy/mm/dd hh:mm:ss'} in $cwd\n-------------------\n$message\n";
    confess($message);
}

=head2 log

 Title   : log
 Usage   : $obj->log('message');
 Function: Appends message to file set in log_file() if write_logs() is on.
 Returns : n/a
 Args    : log message, optional int to decide the strength of the message, as
           per verbose(), using current value of verbose() as default.

=cut

sub log {
    my ($self, $message, $verbose) = @_;
    return unless $self->write_logs();
    $verbose ||= $self->verbose();
    return unless $verbose >= -1;
    
    my $prefix = '';
    
    $self->{_log_count}++;
    if ($self->{_log_count} == 1) {
        $prefix = "----------\n$time{'yyyy/mm/dd'}\n----------\n";
    }
    
    $prefix .= "$time{'hh:mm:ss'} | ";
    
    # for some reason Carp's shortmess and longmess methods aren't returning
    # the same thing as what carp()/croak() produce! hacky work-around...
    # (we die instead of throw since throw calls log()...)
    my $capture = IO::Capture::Stderr->new();
    $capture->start();
    $verbose >= 1 ? cluck($message) : carp($message);
    $capture->stop();
    
    my $file = $self->log_file;
    my $opened = open(my $fh, ">>", $file);
    if ($opened) {
        print $fh $prefix, $capture->read, "\n";
        close($fh);
    }
    else {
        carp("Unable to write to log file '$file', disabling logging!");
        $self->write_logs(0);
    }
}

=head2 write_logs

 Title   : write_logs
 Usage   : $obj->write_logs(1);
 Function: Turn the writing of logs on or off.
 Returns : boolean (current value for writing logs, true is on)
 Args    : boolean (optional, to set: true means logs should be written)

=cut

sub write_logs {
    my $self = shift;
    
    if (@_) {
        $self->{_write_logs} = shift;
    }
    
    return $self->{_write_logs};
}

=head2 log_file

 Title   : log_file
 Usage   : $obj->log_file($filename);
 Function: Get/set the log filename.
 Returns : path of (potential) log file
 Args    : path of log file to set. Defaults to 

=cut

sub log_file {
    my $self = shift;
    
    if (@_) {
        $self->{_log_file} = shift;
    }
    
    return $self->{_log_file};
}

=head2 register_for_cleanup

 Title   : register_for_cleanup
 Usage   : $self->register_for_cleanup($method_name);
 Function: Store a method that will be called upon object destruction.
 Returns : n/a
 Args    : method name string

=cut

sub register_for_cleanup {
    my ($self, $method) = @_;
    if ($method && $self->can($method)) {
        $self->{'_cleanup_methods'}->{$method} = 1;
    }
}

=head2 unregister_for_cleanup

 Title   : unregister_for_cleanup
 Usage   : $self->unregister_for_cleanup($method_name);
 Function: Unstore a method previously registered so that it won't be called
           upon object destruction.
 Returns : n/a
 Args    : method name string

=cut

sub unregister_for_cleanup {
    my ($self, $method) = @_;
    delete $self->{'_cleanup_methods'}->{$method};
}

=head2 register_for_unlinking

 Title   : register_for_unlinking
 Usage   : $self->register_for_unlinking($filename);
 Function: Store a filepath that will be unlinked upon object destruction.
 Returns : n/a
 Args    : filepath(s)

=cut

sub register_for_unlinking {
    my ($self, @files) = @_;
    foreach my $file (@files) {
        $self->{'_unlink_files'}->{$file} = 1;
    }
}

=head2 unregister_for_unlinking

 Title   : unregister_for_unlinking
 Usage   : $self->unregister_for_unlinking($filename);
 Function: Unstore a file previously registered so that it won't be unlinked
           upon object destruction.
 Returns : n/a
 Args    : filepath(s)

=cut

sub unregister_for_unlinking {
    my ($self, @files) = @_;
    foreach my $file (@files) {
        delete $self->{'_unlink_files'}->{$file};
    }
}

=head2 DESTROY

 Title   : DESTROY
 Usage   : n/a
 Function: Our DESTROY method calls all our cleanup methods and unlinks
           registered files.
 Returns : n/a
 Args    : n/a

=cut

sub DESTROY {
    my $self = shift;
    
    my @cleanup_methods = keys %{$self->{'_cleanup_methods'} || {}};
    foreach my $method (@cleanup_methods) {
        $self->$method();
    }
    
    my @files_to_unlink = keys %{$self->{'_unlink_files'} || {}};
    foreach my $file (@files_to_unlink) {
        unlink($file);
    }
}

=head2 Signal handling

    By default, DESTROY is called for each object when __DIE__, SIGINT, or SIGTERM 
    is received. In case that the default DESTROY handlers are not flexible enough,
    do NOT override $SIG{INT} or $SIG{TERM} handlers. Either create a child of Base
    and modify the DESTROY handler, or add static set_die_handler and unset_die_handler
    calls to Base. 
    
    The signals should be probably trapped only when requested, not by default. Maybe
    the subroutine register_for_unlinking should do this, as it this was introduced
    only to really unlink certain files?

=cut

our $SIGNAL_CAUGHT_EVENT = "Signal caught.\n";
$SIG{TERM} = sub { die $SIGNAL_CAUGHT_EVENT; };   # Pipeline.pm and pipeline script rely on this string.
$SIG{INT}  = sub { die $SIGNAL_CAUGHT_EVENT; };

1;
