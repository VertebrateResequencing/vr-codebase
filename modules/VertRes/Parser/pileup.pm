=head1 NAME

VertRes::Parser::pileup - parse pileup format alignments

=head1 SYNOPSIS

use VertRes::Parser::pileup;

# create object, supplying pileup output file or filehandle
my $pars = VertRes::Parser::pileup->new(file => 'pileup.txt');

# get the array reference that will hold the most recently requested result
my $result_holder = $pars->result_holder();

# loop through the output, getting results
while ($pars->next_result()) {
    ...
}

=head1 DESCRIPTION

A straightforward parser for pileup alignments.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Parser::pileup;

use strict;
use warnings;

use base qw(VertRes::Parser::ParserI);


=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Parser::pileup->new(file => 'filename');
 Function: Build a new VertRes::Parser::pileup object.
 Returns : VertRes::Parser::pileup object
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
           [0] sequence name
           [1] base position
           [2] reference base
           [3] number of reads covering the site
           [4] read bases
           [5] base qualities
           (for more info on [4] and [5], see
            http://samtools.sourceforge.net/pileup.shtml)
 Args    : n/a

=cut

=head2 next_result

 Title   : next_result
 Usage   : while ($obj->next_result()) { # look in result_holder }
 Function: Parse the next line from the pileup file.
 Returns : boolean (false at end of output; check the result_holder for the
           actual result information)
 Args    : n/a

=cut

sub next_result {
    my $self = shift;
    
    # get the next line
    my $fh = $self->fh() || return;
    my $line = <$fh> || return;
    
    my @data = split(qr/\s+/, $line);
    @data || return;
    
    for my $i (0..$#data) {
        $self->{_result_holder}->[$i] = $data[$i];
    }
    
    return 1;
}

1;
