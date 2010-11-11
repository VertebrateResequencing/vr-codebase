=head1 NAME

VertRes::Parser::LSF - parse LSF output files

=head1 SYNOPSIS

use VertRes::Parser::LSF;

my $parser = VertRes::Parser::LSF->new(file => 'job.o');
$parser->get('status');     # the status of the last executed job
$parser->get('time');       # running time
$parser->get('memory');     # the memory in MB

# Query all records: 0 is the first run, $n-1 the last.
my $n = $parser->nrecords();
for (my $i=0; $i<$n; $i++) 
{
    print $parser->get('status',$i) . "\n";
}

=head1 DESCRIPTION

A parser for LSF output files. 

=head1 AUTHOR

petr.danecek@sanger

=cut

package VertRes::Parser::LSF;

use strict;
use warnings;
use base qw(VertRes::Base);

sub new 
{
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    if ( !$$self{file} ) { $self->throw("Missing the 'file' parameter."); }
    open(my $fh,'<',$$self{file}) or $self->throw("$$self{file}: $!");
    $self->_parse($fh);
    close($fh);
    
    return $self;
}

# Currently collects/looks at the following:
#
#   Start of the record:
#       Your job looked like:
#   Status:
#       Successfully completed
#       Exited with exit code 9.
#       TERM_OWNER: job killed by owner.
#       TERM_RUNLIMIT: job killed after reaching LSF run time limit.
#       TERM_MEMLIMIT: job killed after reaching LSF memory usage limit.
#   Collects:
#       CPU time   :    170.19 sec.
#       Max Memory :       251 MB
#
sub _parse
{
    my ($self,$fh) = @_;

    my @results;
    my $record = {};
    while (my $line=<$fh>)
    {
        # New records are recognised by "Your job looked like:". 
        #   If anything above is needed, please modify.
        if ( $line=~/^Your job looked like:/ && scalar keys %$record ) 
        { 
            push @results, $record;
            $record = {};
            next;
        }
        
        if ( $line=~/^Successfully completed./ ) { $$record{status} = 'OK'; next; }
        elsif ( !exists($$record{status}) && $line=~/^Exited with exit code/ ) { $$record{status} = 'exited'; next; }
        elsif ( $line=~/^TERM_\S+ job killed by/ ) { $$record{status} = 'killed'; next; }
        elsif ( $line=~/^TERM_([^:]+):/ ) { $$record{status} = $1; next; }
        elsif ( $line=~/^\s+CPU time\s+:\s+(\S+)/ ) { $$record{time}=$1; next; }
        elsif ( $line=~/^\s+Max Memory\s+:\s+(\S+)\s+(\S+)/ ) 
        { 
            $$record{memory}=$1;
            if ( $2 eq 'KB' ) { $$record{memory} /= 1024; }
            elsif ( $2 eq 'GB' ) { $$record{memory} *= 1024; }
            next; 
        }
    }
    if ( scalar keys %$record ) { push @results,$record; }
    $$self{results} = \@results;
}

sub get
{
    my ($self,$key,$idx) = @_;
    if ( !defined($idx) ) { $idx=-1; }
    if ( !exists($$self{results}) or !scalar @{$$self{results}} ) { $self->throw("No records in the LSF output file? [$$self{file}]"); }
    my $record = $$self{results}[$idx];
    if ( !defined $record ) { $self->throw("The index out of bounds: $idx [$$self{file}]"); }
    if ( !exists($$record{$key}) ) { $self->throw("The key '$key' not present for [$idx] in [$$self{file}]"); }
    return $$record{$key};
}

sub nrecords
{
    my ($self) = @_;
    if ( !exists($$self{results}) ) { return 0; }
    return scalar @{$$self{results}};
}

1;
