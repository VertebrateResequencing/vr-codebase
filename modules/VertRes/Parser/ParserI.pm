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

=head2 _save_position

 Title   : _save_position
 Usage   : $self->_save_position()
 Function: Internal method for parser authors. Saves the current filehandle
           position; for use before seeking.
 Returns : boolean (true if this is a seekable filehandle)
 Args    : n/a

=cut

sub _save_position {
    my $self = shift;
    
    my $fh = $self->fh() || return;
    
    my $tell = tell($fh);
    if ($tell == -1) {
        $self->warn("this parsing method doesn't work on piped input");
        return;
    }
    my $current_results = ref($self->{_result_holder}) eq 'ARRAY' ? [@{$self->{_result_holder}}] : {%{$self->{_result_holder}}};
    
    $self->{_tell} = $tell;
    $self->{_current_results} = $current_results;
    
    return 1;
}

=head2 _get_header

 Title   : _get_header
 Usage   : $self->_get_header()
 Function: Internal method for parser authors. This does nothing. If your file
           format has a header, you should implement it by overriding this.
           In your implementation, _set_header_parsed() should be called after
           successfully parsing your header.
           It should return true if the header was parsed, false if not.
 Returns : boolean
 Args    : n/a

=cut

sub _get_header {
    my $self = shift;
    unless ($self->_header_parsed) {
        $self->_set_header_parsed(0);
    }
    return 0;
}

=head2 _header_parsed

 Title   : _header_parsed
 Usage   : if ($self->_header_parsed()) { #... }
 Function: Internal method for parser authors. Ask if the header of the current
           file has been parsed.
 Returns : boolean
 Args    : n/a

=cut

sub _header_parsed {
    my $self = shift;
    
    my $fh = $self->fh() || return;
    my $fh_id = $self->_fh_id;
    
    return defined $self->{'_got_header'.$fh_id};
}

=head2 _set_header_parsed

 Title   : _set_header_parsed
 Usage   : $self->_set_header_parsed()
 Function: Internal method for parser authors. Ask if the header of the current
           file has been parsed.
 Returns : boolean
 Args    : n/a (optionally, supply a filehandle position as from tell() that
           denotes the start of the first record after the header; by default
           it assumes the current filehandle position matches this criteria)

=cut

sub _set_header_parsed {
    my ($self, $tell) = @_;
    
    my $fh = $self->fh() || return;
    my $fh_id = $self->_fh_id;
    
    unless (defined $tell) {
        $tell = tell($fh);
    }
    
    $self->{'_got_header'.$fh_id} = $tell;
}

=head2 _seek_first_result

 Title   : _seek_first_result
 Usage   : $self->_seek_first_result()
 Function: Internal method for parser authors. Seeks back to before the first
           result (ie. after a header if present) so that next_result() will
           behave as if it was called for the first time on a file.
           You should call _save_position() first, then call this, then do your
           work, then call _restore_position().
 Returns : n/a
 Args    : n/a

=cut

sub _seek_first_result {
    my $self = shift;
    
    my $fh = $self->fh() || return;
    my $fh_id = $self->_fh_id;
    
    $self->_get_header();
    my $tell;
    if (defined $self->{'_got_header'.$fh_id}) {
        $tell = $self->{'_got_header'.$fh_id};
        if ($tell == -1) {
            $self->warn("this parsing method doesn't work on piped input");
            return;
        }
    }
    else {
        $tell = 0;
    }
    
    $self->seek($tell, 0);
}

=head2 _restore_position

 Title   : _restore_position
 Usage   : $self->_restore_position()
 Function: Internal method for parser authors. Restores the current filehandle
           position to the location when _save_position() was last called; for
           use after you've seeked somewhere and done some work.
           If _save_position() was called at the true start of the file, this
           actually calls _seek_first_result().
 Returns : n/a
 Args    : n/a

=cut

sub _restore_position {
    my $self = shift;
    
    $self->fh() || return;
    
    if ($self->{_tell} == 0) {
        # we might have saved position before parsing the header, but now have
        # the flag set that we've parsed the header; don't seek back before
        # the header!
        $self->_seek_first_result();
        return;
    }
    
    $self->seek($self->{_tell}, 0);
    if (ref($self->{_result_holder}) eq 'ARRAY') {
        my @current_results = @{$self->{_current_results}};
        for my $i (0..$#current_results) {
            $self->{_result_holder}->[$i] = $current_results[$i];
        }
    }
    else {
        while (my ($key, $val) = each %{$self->{_current_results}}) {
            $self->{_result_holder}->{$key} = $val;
        }
    }
    
    return 1;
}

1;
