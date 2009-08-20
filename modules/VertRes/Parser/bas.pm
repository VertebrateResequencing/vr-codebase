=head1 NAME

VertRes::Parser::bas - parse bas files

=head1 SYNOPSIS

use VertRes::Parser::bas;

# create object, supplying bas file or filehandle
my $pars = VertRes::Parser::bas->new(file => 'raw.bam.bas');

# get the array reference that will hold the most recently requested result
my $result_holder = $pars->result_holder();

# loop through the output, getting results
while ($pars->next_result()) {
    my $total_bases = $result_holder->[6];
    # etc.
}

=head1 DESCRIPTION

A straightforward parser for bas files, as produced by
VertRes::Utils::Sam->bas().

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Parser::bas;

use strict;
use warnings;

use base qw(VertRes::Parser::ParserI);

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Parser::bas->new(file => 'filename');
 Function: Build a new VertRes::Parser::bas object.
 Returns : VertRes::Parser::bas object
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
           [0] bam_filename
           [1] project
           [2] sample
           [3] platform
           [4] library
           [5] readgroup
           [6] #_total_bases
           [7] #_mapped_bases
           [8] #_total_reads
           [9] #_mapped_reads
           [10] #_mapped_reads_paired_in_sequencing
           [11] #_mapped_reads_properly_paired
           [12] %_of_mismatched_bases
           [13] average_quality_of_mapped_bases
           [14] mean_insert_size
           [15] insert_size_sd
           [16] median_insert_size
           [17] insert_size_median_absolute_deviation
 Args    : n/a

=cut

=head2 next_result

 Title   : next_result
 Usage   : while ($obj->next_result()) { # look in result_holder }
 Function: Parse the next line from the bas file.
 Returns : boolean (false at end of output; check the result_holder for the
           actual result information)
 Args    : n/a

=cut

sub next_result {
    my $self = shift;
    
    # get the next line
    my $fh = $self->fh() || return;
    my $line = <$fh> || return;
    chomp($line);
    
    my @data = split(qr/\t/, $line);
    @data || return;
    
    # header?
    if ($data[17] eq 'insert_size_median_absolute_deviation') {
        return $self->next_result;
    }
    
    for my $i (0..$#data) {
        $self->{_result_holder}->[$i] = $data[$i];
    }
    
    return 1;
}

1;
