=head1 NAME

VertRes::Parser::ParserI - interface for parsing files

=head1 SYNOPSIS

use base qw(VertRes::Parser::ParserI);

=head1 DESCRIPTION

Make a parser that takes a file or filehandle and will later go through the
filehandle, storing desired data in result_holder().

This system is built for speed and low memory. To avoid object creation and
method call overhead you call for the result data structure once, and the
results within are updated as you loop through results with the next_result()
method.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Parser::ParserI;

use strict;
use warnings;

use base qw(VertRes::IO);

=head2 new

 Title   : new
 Usage   : my $self = $class->SUPER::new(@args);
 Function: Instantiates your object, opening user's file/ filehandle.
 Returns : $self hash-ref blessed into your class
 Args    : file => filename -or- fh => filehandle

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    # we'll store the latest results in this data structure
    $self->{_result_holder} = [];
    
    return $self;
}

=head2 result_holder

 Title   : result_holder
 Usage   : my $result_holder = $obj->result_holder()
 Function: Get the data structure that will hold the last result requested by
           next_result()
 Returns : array ref
 Args    : n/a

=cut

sub result_holder {
    my $self = shift;
    return $self->{_result_holder};
}

=head2 next_result

 Title   : next_result
 Usage   : while ($obj->next_result()) { # look in result_holder }
 Function: Parse the next result.
 Returns : boolean (false at end of output; check the result_holder for the
           actual result information)
 Args    : n/a

=cut

sub next_result {
    my $self = shift;
    $self->throw("This is just an interface; this method should have been overridden");
}

sub _save_position {
    my $self = shift;
    
    my $fh = $self->fh() || return;
    
    my $tell = tell($fh);
    if ($tell == -1) {
        $self->warn("this parsing method doesn't work on piped input");
        return;
    }
    my @current_results = @{$self->{_result_holder}};
    
    $self->{_tell} = $tell;
    $self->{_current_results} = \@current_results;
    
    return 1;
}

sub _restore_position {
    my $self = shift;
    
    my $fh = $self->fh() || return;
    my @current_results = @{$self->{_current_results}};
    
    $self->seek($fh, $self->{_tell}, 0);
    for my $i (0..$#current_results) {
        $self->{_result_holder}->[$i] = $current_results[$i];
    }
    
    return 1;
}

1;
