=head1 NAME

VertRes::Parser::LSF - parse LSF output files

=head1 SYNOPSIS

use VertRes::Parser::LSF;

my $parser = VertRes::Parser::LSF->new(file => 'job.o');

# How many LSF reports are there in the file?
my $n = $parser->nrecords();

# Loop through the records and get status and memory
for (my $i=0; $i<$n; $i++) {
    my $status = $parser->get('status',$i);
    my $memory = $parser->get('memory',$i);
    my $queue  = $parser->get('queue',$i);
}


# shortcuts to get the most recent stats in the LSF file 
my $status = $parser->status;     # the status of the last executed job
my $wall_time = $parser->time;    # wall time (s)
my $cpu_time = $parser->cpu_time; # cpu time (s)
my $idle = $parser->idle_factor;  # the idle factor
my $memory = $parser->memory;     # the memory (MB)
my $cmd = $parser->cmd;           # the command line that was executed

# get the array reference that will hold the most recently requested result
my $result_holder = $pars->result_holder();

# loop through the LSF output, getting results
while ($parser->next_result()) {
    my $time = $result_holder->[3];
    # etc.
}

=head1 DESCRIPTION

A parser for LSF output files. 

=head1 AUTHOR

petr.danecek@sanger & bix@sendu.me.uk

=cut

package VertRes::Parser::LSF;

use strict;
use warnings;
use base qw(VertRes::Parser::ParserI);
use DateTime;

our %months = qw(Jan 1
                 Feb 2
                 Mar 3
                 Apr 4
                 May 5
                 Jun 6
                 Jul 7
                 Aug 8
                 Sep 9
                 Oct 10
                 Nov 11
                 Dec 12);

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Parser::LSF->new(file => 'filename');
 Function: Build a new VertRes::Parser::LSF object.
 Returns : VertRes::Parser::LSF object
 Args    : file => filename -or- fh => filehandle

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    return $self;
}

=head2 result_holder

 Title   : result_holder
 Usage   : my $result_holder = $obj->result_holder()
 Function: Get the data structure that will hold the last result requested by
           next_result()
 Returns : array ref, where the elements are:
           [0]  cmd
           [1]  status
           [2]  memory
           [3]  time
           [4]  cpu_time
           [5]  idle_factor
           [6]  queue
 Args    : n/a

=cut

=head2 next_result

 Title   : next_result
 Usage   : while ($obj->next_result()) { # look in result_holder }
 Function: Parse the next line from the sequence.index.
 Returns : boolean (false at end of output; check the result_holder for the
           actual result information)
 Args    : n/a

=cut

sub next_result {
    my $self = shift;
    
    # just return if no file set
    my $fh = $self->fh() || return;
    my $fh_id = $self->_fh_id;
    
    # we're only interested in the small LSF-generated report, which may be
    # prefaced by an unlimited amount of output from the program that LSF ran,
    # so we go line-by-line to find our little report
    my ($found_report_start, $found_report_end, $next_is_cmd);
    my ($started, $finished, $cmd, $mem, $status,$queue);
    my $cpu = 0;
    while (<$fh>) {
        if (/^Sender: LSF System/) {
            $found_report_start = 1;
        }
        elsif (/^The output \(if any\) is above this job summary/) {
            $found_report_end = 1;
            last;
        }
        
        if ($found_report_start) {
            if (/^Started at \S+ (.+)$/) { $started = $1; }
            elsif (/^Job was executed.+in queue \<([^>]+)\>/) { $queue = $1; }
            elsif (/^Results reported at \S+ (.+)$/) { $finished = $1; }
            elsif (/^# LSBATCH: User input/) { $next_is_cmd = 1; }
            elsif ($next_is_cmd) {
                if (/^-{60}$/) {
                    $next_is_cmd = 0;
                }
                else {
                    $cmd .= $_;
                }
            }
            elsif (/^Successfully completed/) { $status = 'OK'; }
            elsif (/^Cannot open your job file/) { $status = 'unknown'; }
            elsif (! $status && /^Exited with exit code/) { $status = 'exited'; }
            elsif (/^TERM_\S+ job killed by/) { $status = 'killed'; }
            elsif (/^TERM_([^:]+):/) { $status = $1; }
            elsif (/^\s+CPU time\s+:\s+(\S+)/) { $cpu = $1; }
            elsif (/^\s+Max Memory\s+:\s+(\S+)\s+(\S+)/) { 
                $mem = $1;
                if ($2 eq 'KB') { $mem /= 1024; }
                elsif ($2 eq 'GB') { $mem *= 1024; }
            }
        }
    }
    
    # if we didn't see a whole LSF report, assume eof
    unless ($found_report_start && $found_report_end) {
        $self->{"saw_last_record_$fh_id"} = 1;
        return;
    }
    unless ($status) {
        warn "a status was not parsed out of a result in ".$self->file."\n";
    }
    
    chomp($cmd) if $cmd;
    
    # calculate wall time and idle factor
    my $date_regex = qr/(\w+)\s+(\d+) (\d+):(\d+):(\d+)/;
    my ($smo, $sd, $sh, $sm, $ss) = $started =~ /$date_regex/;
    my ($emo, $ed, $eh, $em, $es) = $finished =~ /$date_regex/;
    my $dt = DateTime->new(year => 2010, month => $months{$smo}, day => $sd, hour => $sh, minute => $sm, second => $ss);
    my $st = $dt->epoch;
    $dt = DateTime->new(year => 2010, month => $months{$emo}, day => $ed, hour => $eh, minute => $em, second => $es);
    my $et = $dt->epoch;
    my $wall = $et - $st;
    my $idle = sprintf("%0.2f", ($cpu < 1 ? 1 : $cpu) / ($wall < 1 ? 1 : $wall));
    
    # fill in both the result_holder and maintain in memory all results in a
    # seperate data structure so get() will work the old way
    $self->{_result_holder}->[0] = $cmd;
    $self->{_result_holder}->[1] = $status;
    $self->{_result_holder}->[2] = $mem;
    $self->{_result_holder}->[3] = $wall;
    $self->{_result_holder}->[4] = $cpu;
    $self->{_result_holder}->[5] = $idle;
    $self->{_result_holder}->[6] = $queue;
    unless ($self->{"saw_last_record_$fh_id"}) {
        push(@{$self->{results}}, {cmd => $cmd,
                                   status => $status,
                                   memory => $mem,
                                   time => $wall,
                                   cpu_time => $cpu,
                                   idle_factor => $idle,
                                   queue => $queue});
    }
    
    return 1;
}

sub get {
    my ($self,$key,$idx) = @_;
    if ( !defined($idx) ) { $idx=-1; }
    if ( !exists($$self{results}) or !scalar @{$$self{results}} ) {
        my $fh = $self->fh() || return;
        $self->_save_position;
        my $fh_id = $self->_fh_id;
        unless ($self->{"saw_last_record_$fh_id"}) {
            while ($self->next_result) {
                next;
            }
        }
        $self->_restore_position;
        unless (defined $self->{results}) {
            $self->warn("No records in the LSF output file? [$$self{file}]");
            return;
        }
    }
    if ( !defined $key ) { return undef; }  # This can really happen, see nrecords().
    my $record = $$self{results}[$idx];
    if ( !defined $record ) { $self->throw("The index out of bounds: $idx [$$self{file}]"); }
    if ( !exists($$record{$key}) ) { $self->throw("The key '$key' not present for [$idx] in [$$self{file}]"); }
    return $$record{$key};
}

sub nrecords {
    my ($self) = @_;
    if ( !exists($$self{results}) ) { $self->get(); }
    return scalar @{$self->{results} || []};
}

=head2 status

 Title   : status
 Usage   : my $status = $obj->status();
 Function: Get the status of the last job reported in the LSF file.
 Returns : string (OK|exited|killed|MEMLIMIT|RUNLIMIT)
 Args    : n/a

=cut

sub status {
    my $self = shift;
    return $self->get('status', -1);
}

=head2 time

 Title   : time
 Usage   : my $time = $obj->time();
 Function: Get the wall-time of the last job reported in the LSF file.
 Returns : int (s)
 Args    : n/a

=cut

sub time {
    my $self = shift;
    return $self->get('time', -1);
}

=head2 cpu_time

 Title   : cpu_time
 Usage   : my $time = $obj->cpu_time();
 Function: Get the cpu-time of the last job reported in the LSF file.
 Returns : real number (s)
 Args    : n/a

=cut

sub cpu_time {
    my $self = shift;
    return $self->get('cpu_time', -1);
}

=head2 idle_factor

 Title   : idle_factor
 Usage   : my $idle_factor = $obj->idle_factor();
 Function: Compare cpu time to wall time to see what proportion was spent
           waiting on disc.
 Returns : real number 0-1 inclusive (1 means no time spent waiting on disc, 0
           means the cpu did nothing and we spent all time waiting on disc)
 Args    : n/a

=cut

sub idle_factor {
    my $self = shift;
    return $self->get('idle_factor', -1);
}

=head2 memory

 Title   : memory
 Usage   : my $memory = $obj->memory();
 Function: Get the max memory used of the last job reported in the LSF file.
 Returns : int (s)
 Args    : n/a

=cut

sub memory {
    my $self = shift;
    return $self->get('memory', -1);
}

=head2 cmd

 Title   : cmd
 Usage   : my $cmd = $obj->cmd();
 Function: Get the command-line of the last job reported in the LSF file.
 Returns : string
 Args    : n/a

=cut

sub cmd {
    my $self = shift;
    return $self->get('cmd', -1);
}

=head2 queue

 Title   : queue
 Usage   : my $queue = $obj->queue();
 Function: Get the command-line of the last job reported in the LSF file.
 Returns : string
 Args    : n/a

=cut

sub queue {
    my $self = shift;
    return $self->get('queue', -1);
}

1;
