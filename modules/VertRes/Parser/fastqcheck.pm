=head1 NAME

VertRes::Parser::fastqcheck - parse fastqcheck files

=head1 SYNOPSIS

use VertRes::Parser::fastqcheck;

# create object, supplying fastqcheck output file or filehandle
my $pars = VertRes::Parser::fastqcheck->new(file => 'read.fastq.gz.fastqcheck');

# get basic summary info:
my $num_of_sequences = $pars->num_sequences();
my $total_length_of_all_sequences = $pars->total_length();
my $average_length_of_a_sequence = $pars->avg_length();
my $length_of_longest_sequence = $pars->max_length();
my ($total_sd, $per_base_sd) = $pars->standard_deviations();
my ($bases,$avg_quals) = $pars->avg_base_quals();

# get the array reference that will hold the most recently requested result
my $result_holder = $pars->result_holder();

# loop through the output, getting results
while ($pars->next_result()) {
    my $position = $result_holder->[0];
    my $a_percent = $result_holder->[1];
    # etc.
}

=head1 DESCRIPTION

A straightforward parser for fastqcheck files.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Parser::fastqcheck;

use strict;
use warnings;

use base qw(VertRes::Parser::ParserI);


=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Parser::fastqcheck->new(file => 'filename');
 Function: Build a new VertRes::Parser::fastqcheck object.
 Returns : VertRes::Parser::fastqcheck object
 Args    : file => filename -or- fh => filehandle

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    return $self;
}

=head2 num_sequences

 Title   : num_sequences
 Usage   : my $num_of_sequences = $obj->num_sequences();
 Function: Get the number of sequences that the fastqcheck file is summarising.
 Returns : int
 Args    : n/a (new() or file() must allready have been supplied with a
                filename or filehandle)

=cut

sub num_sequences {
    my $self = shift;
    return $self->_getter('num_sequences');
}

=head2 total_length

 Title   : total_length
 Usage   : my $total_length_of_all_sequences = $obj->total_length();
 Function: Get the total length of sequences that the fastqcheck file is
           summarising.
 Returns : int
 Args    : n/a (new() or file() must allready have been supplied with a
                filename or filehandle)

=cut

sub total_length {
    my $self = shift;
    return $self->_getter('total_length');
}

=head2 avg_length

 Title   : avg_length
 Usage   : my $average_length_of_a_sequence = $obj->avg_length();
 Function: Get the average length of a sequence.
 Returns : number
 Args    : n/a (new() or file() must allready have been supplied with a
                filename or filehandle)

=cut

sub avg_length {
    my $self = shift;
    return $self->_getter('avg_length');
}

=head2 max_length

 Title   : max_length
 Usage   : my $length_of_longest_sequence = $obj->max_length();
 Function: Get the length of the longest sequence.
 Returns : int
 Args    : n/a (new() or file() must allready have been supplied with a
                filename or filehandle)

=cut

sub max_length {
    my $self = shift;
    return $self->_getter('max_length');
}

=head2 standard_deviations

 Title   : standard_deviations
 Usage   : my ($total_sd, $per_base_sd) = $obj->standard_deviations();
 Function: Get the standard deviations at 0.25; total and per-base.
 Returns : list: two percentages
 Args    : n/a (new() or file() must allready have been supplied with a
                filename or filehandle)

=cut

sub standard_deviations {
    my $self = shift;
    return @{$self->_getter('standard_deviations')};
}

sub _getter {
    my ($self, $var) = @_;
    
    $self->_get_header;
    unless (defined $self->{"_$var"}) {
        $self->throw("$var unknown - did you supply a file and was it a fastqcheck file?");
    }
    return $self->{"_$var"};
}

sub _get_header {
    my $self = shift;
    
    my $fh = $self->fh() || return;
    return 1 if $self->_header_parsed();
    
    my $saw = 0;
    while (<$fh>) {
        if (/^(\d+) sequences, (\d+) total length, (\d+\.\d+) average, (\d+) max/) {
            $self->{_num_sequences} = $1;
            $self->{_total_length} = $2;
            $self->{_avg_length} = $3;
            $self->{_max_length} = $4;
            $saw++;
        }
        elsif (/^Standard deviations at 0.25:  total  (\d+\.\d+) %, per base  (\d+\.\d+) %/) {
            $self->{_standard_deviations} = [$1, $2];
            $saw++;
        }
        elsif (/^\s+A\s+C/) {
            $saw++;
            last;
        }
    }
    
    if ($saw == 3) {
        $self->_set_header_parsed();
        return 1;
    }
    return;
}

=head2 result_holder

 Title   : result_holder
 Usage   : my $result_holder = $obj->result_holder()
 Function: Get the data structure that will hold the last result requested by
           next_result()
 Returns : array ref, where the elements are:
           [0] base position (starting at position 0 == 'total')
           [1] A %
           [2] C %
           [3] G %
           [4] T %
           [5] N %
           [6] count of bases at position [0] with quality 0
           [7] count of bases at position [0] with quality 1
           ... etc. for quality up to 40 (or higher?)
           [-1] the 'AQ' which is NOT the average quality... but what is it?!
 Args    : n/a

=cut

=head2 next_result

 Title   : next_result
 Usage   : while ($obj->next_result()) { # look in result_holder }
 Function: Parse the next base position line from the fastqcheck output.
 Returns : boolean (false at end of output; check the result_holder for the
           actual result information)
 Args    : n/a

=cut

sub next_result {
    my $self = shift;
    
    my $fh = $self->fh() || return;
    
    $self->_get_header() || $self->throw("Unable to parse header before first result - is this a fastqcheck file?");
    
    # get the next line
    my $line = <$fh> || return;
    
    my @data = split(qr/\s+/, $line);
    @data || return;
    
    my $r = 0;
    for my $i (0..$#data) {
        my $datum = $data[$i];
        if ($i == 0) {
            next if $datum eq 'base';
            $datum = 0 if $datum eq 'Total';
        }
        $self->{_result_holder}->[$r++] = $datum;
    }
    
    return 1;
}

=head2 avg_base_quals

 Title   : avg_base_quals
 Usage   : my ($bases,$quals) = $pars->avg_base_quals();
 Function: Get the average qualities for each base.
 Returns : list with two array references
 Args    : n/a (new() or file() must allready have been supplied with a
                filename or filehandle)

=cut

sub avg_base_quals {
    my $self = shift;
    
    $self->_save_position || return;
    $self->_seek_first_result();
    my $rh = $self->result_holder();
    
    $self->next_result();
    
    my (@bases, @quals);
    while ($self->next_result()) {
        push @bases, $rh->[0];
        push @quals, $self->_avg_base_qual;
    }
    
    $self->_restore_position;
    
    return (\@bases, \@quals);
}

=head2 avg_qual

 Title   : avg_qual
 Usage   : my $avg_qual = $pars->avg_qual();
 Function: Get the average quality of all bases.
 Returns : float
 Args    : n/a (new() or file() must allready have been supplied with a
                filename or filehandle)

=cut

sub avg_qual {
    my $self = shift;
    
    $self->_save_position || return;
    $self->_seek_first_result();
    
    $self->next_result();
    
    my $avg_qual = $self->_avg_base_qual;
    
    $self->_restore_position;
    
    return $avg_qual;
}

sub _avg_base_qual {
    my $self = shift;
    my $rh = $self->result_holder();
    
    my $sum    = 0;
    my $nvals  = 0;
    my $nquals = $#{$rh} - 1;
    for my $i (6..$nquals) {
        $nvals += $rh->[$i];
        $sum   += $rh->[$i] * ($i-6);
    }
    my $avg_qual = $nvals ? $sum / $nvals : 0;
    
    return $avg_qual;
}

1;
