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
    my $total_bases = $result_holder->[7];
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
           [1] bam md5 checksum
           [2] study
           [3] sample
           [4] platform
           [5] library
           [6] readgroup
           [7] #_total_bases
           [8] #_mapped_bases
           [9] #_total_reads
           [10] #_mapped_reads
           [11] #_mapped_reads_paired_in_sequencing
           [12] #_mapped_reads_properly_paired
           [13] %_of_mismatched_bases
           [14] average_quality_of_mapped_bases
           [15] mean_insert_size
           [16] insert_size_sd
           [17] median_insert_size
           [18] insert_size_median_absolute_deviation
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
    
    # old bas files didn't have md5 checksum in column 2, and oldest didn't
    # have bam filename in column 1
    if (@data == 19) {
        for my $i (0..$#data) {
            $self->{_result_holder}->[$i] = $data[$i];
        }
    }
    elsif (@data == 18) {
        my $j = 0;
        for my $i (0..$#data) {
            if ($i == 1) {
                $j++;
            }
            $self->{_result_holder}->[$j] = $data[$i];
            $j++;
        }
    }
    elsif (@data == 17) {
        my $j = 2;
        for my $i (0..$#data) {
            $self->{_result_holder}->[$j] = $data[$i];
            $j++;
        }
    }
    else {
        $self->throw("Unexpected number of columns (".scalar(@data)."); is this really a bas file?\n$line");
    }
    
    # header?
    if ($self->{_result_holder}->[18] eq 'insert_size_median_absolute_deviation') {
        # initialise everything to unknown/0 so that if user looks at result
        # holder without checking that next_result returned true, they don't
        # get the header values
        for my $i (0..6) {
            $self->{_result_holder}->[$i] = 'unknown';
        }
        for my $i (7..18) {
            $self->{_result_holder}->[$i] = 0;
        }
        $self->_set_header_parsed() unless $self->_header_parsed();
        return $self->next_result;
    }
    
    return 1;
}

=head2 total_reads

 Title   : total_reads
 Usage   : my $total_reads = $obj->total_reads();
 Function: Get the total reads of all readgroups reported in the bas file.
 Returns : int
 Args    : n/a

=cut

sub total_reads {
    my $self = shift;
    return $self->_total(9);
}

sub _total {
    my ($self, $index) = @_;
    
    my $fh = $self->fh() || return;
    my $fh_id = $self->_fh_id;
    $self->_save_position || return;
    
    $self->_seek_first_result();
    
    my $total = 0;
    my $rh = $self->result_holder;
    while ($self->next_result) {
        $total += $rh->[$index];
    }
    
    $self->_restore_position;
    
    return $total;
}

=head2 mapped_reads

 Title   : mapped_reads
 Usage   : my $mapped_reads = $obj->mapped_reads();
 Function: Get the total mapped reads of all readgroups reported in the bas
           file.
 Returns : int
 Args    : n/a

=cut

sub mapped_reads {
    my $self = shift;
    return $self->_total(10);
}

=head2 percent_mapped

 Title   : percent_mapped
 Usage   : my $percent_mapped = $obj->percent_mapped();
 Function: Get the percent of mapped reads across all readgroups reported in
           the bas file
 Returns : float
 Args    : n/a

=cut

sub percent_mapped {
    my $self = shift;
    my $total = $self->total_reads;
    my $mapped = $self->mapped_reads;
    return (100 / $total) * $mapped;
}

1;
