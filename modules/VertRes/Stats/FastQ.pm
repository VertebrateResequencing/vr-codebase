=head1 NAME

VertRes::Stats::FastQ - get information about fastq files

=head1 SYNOPSIS

use VertRes::Stats::FastQ;

my $sfq = VertRes::Stats::FastQ->new();

# run fastqcheck on your fastq file of interest, then you can supply
# the fastqcheck output file as the arg to these methods:
my $num_of_sequences = $sfq->num_sequences($fastqcheck_file);
my $total_length_of_all_sequences = $sfq->total_length($fastqcheck_file);
my $average_length_of_a_sequence = $sfq->avg_length($fastqcheck_file);
my $length_of_longest_sequence = $sfq->max_length($fastqcheck_file);
my ($total_sd, $per_base_sd) = $sfq->standard_deviations($fastqcheck_file);

my $average_quality = $sfq->avg_quality($fastqcheck_file);

my ($before, $after) = $sfq->changed_quality('qualmapBayesian.txt');

=head1 DESCRIPTION

A module that can get interesting stats about fastq files. Currently only
gets info from fastqcheck and qualmapBayesian.txt files, not from the fastq
files directly.

- could be extended to get some info from fastq files as well...

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Stats::FastQ;

use strict;
use warnings;
use Math::NumberCruncher;

use base qw(VertRes::Base);
use VertRes::Parser::fastqcheck;
use VertRes::Parser::qualmapBayesian;

=head2 new

 Title   : new
 Usage   : my $fqs = VertRes::Stats::FastQ->new();
 Function: Make a new VertRes::Stats::FastQ object.
 Returns : VertRes::Stats::FastQ object
 Args    : n/a

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    return $self;
}

=head2 avg_quality

 Title   : avg_quality
 Usage   : my $average_quality = $obj->avg_quality($fastqcheck_file);
 Function: Get the overall average quality of all sequences in the fastq
           as summarised by fastqcheck.
 Returns : number (1 decimal places)
 Args    : fastqcheck filename

=cut

sub avg_quality {
    my ($self, $fastqcheck_file) = @_;
    
    my $parser = VertRes::Parser::fastqcheck->new(file => $fastqcheck_file);
    $parser->next_result || $self->throw("Unable to get first row from your fastqcheck file");
    my $rh = $parser->result_holder();
    $rh->[0] == 0 or $self->throw("First row I looked at wasn't the Total row!");
    
    my @vals = @{$rh}[6..$#{$rh}];
    my $AQ = pop @vals; # NB, this is NOT average quality
    
    # the numbers in @vals are the counts (in millions I think) of bases
    # at the quality of the subscript (i.e. from 0..x)
    # want average of these
    my ($q_n, $q_tot);
    for (my $i=0; $i < scalar @vals; ++$i){
        $q_n += $vals[$i]; 
        $q_tot += $vals[$i] * $i;
    }
    
    if ($q_n){
        return sprintf("%.1f", $q_tot / $q_n);
    }
    else {
        $self->throw("Didn't manage to find the avg_quality! Don't know why...");
    }
}

=head2 changed_quality

 Title   : changed_quality
 Usage   : my ($before, $after) = $obj->changed_quality('qualmapBayesian.txt');
 Function: Get the average quality per position before and after recalibration.
 Returns : two hash refs. The first hash is for the 'before', the
           second hash is for 'after'. Each hash will have (up to) 2 keys which
           identify the fastq file looked at (fastq_1 or fastq_2). The values
           are array refs. Those arrays contains the average quality at each
           position (element 0 being base position 1 etc.).
 Args    : qualmapBayesian.txt found in sample subdir following a recalibration,
           optional name prefix to identify the fastq files the qualmapBayesian
           file is reporting on. (defaults to 'fastq')

=cut

sub changed_quality {
    my ($self, $qmb_file, $name) = @_;
    $name ||= 'fastq';
    
    my $parser = VertRes::Parser::qualmapBayesian->new(file => $qmb_file);
    my $rh = $parser->result_holder();
    
    my @quals_per_position;
    while ($parser->next_result) {
        # ($file is 0, 1 or 2 so we can just use an array index)
        my ($file, $pos, $old_qual, $new_qual, $count, $expected_error) = @{$rh};
        my $old_a = $quals_per_position[$file]->[0]->[$pos] || [];
        my $new_a = $quals_per_position[$file]->[1]->[$pos] || [];
        $old_a->[0] += $old_qual * $count;
        $old_a->[1] += $count;
        $new_a->[0] += $new_qual * $count;
        $new_a->[1] += $count;
        $quals_per_position[$file]->[0]->[$pos] = $old_a;
        $quals_per_position[$file]->[1]->[$pos] = $new_a;
    }
    
    # find the avg qual per position and transform to output data structures
    my (%before, %after);
    foreach my $file (0, 1, 2) {
        my $old_ref = $quals_per_position[$file]->[0] || next;
        my $new_ref = $quals_per_position[$file]->[1];
        my (@before, @after);
        foreach my $pos (1..$#{$old_ref}) {
            push(@before, sprintf("%0.2f", $old_ref->[$pos]->[0] / $old_ref->[$pos]->[1]));
            push(@after,  sprintf("%0.2f", $new_ref->[$pos]->[0] / $new_ref->[$pos]->[1]));
        }
        $before{"${name}_${file}"} = \@before;
        $after{"${name}_${file}"} = \@after;
    }
    
    return (\%before, \%after);
}

=head2 num_sequences

 Title   : num_sequences
 Usage   : my $num_of_sequences = $obj->num_sequences($fastqcheck_file);
 Function: Get the number of sequences that the fastqcheck file is summarising.
 Returns : int
 Args    : fastqcheck filename

=cut

sub num_sequences {
    my ($self, $file) = @_;
    return $self->_fastqcheck_getter('num_sequences', $file);
}

=head2 total_length

 Title   : total_length
 Usage   : my $total_length_of_all_sequences = $obj->total_length($fqc_file);
 Function: Get the total length of sequences that the fastqcheck file is
           summarising.
 Returns : int
 Args    : fastqcheck filename

=cut

sub total_length {
    my ($self, $file) = @_;
    return $self->_fastqcheck_getter('total_length', $file);
}

=head2 avg_length

 Title   : avg_length
 Usage   : my $average_length_of_a_sequence = $obj->avg_length($fqc_file);
 Function: Get the average length of a sequence.
 Returns : number
 Args    : fastqcheck filename

=cut

sub avg_length {
    my ($self, $file) = @_;
    return $self->_fastqcheck_getter('avg_length', $file);
}

=head2 max_length

 Title   : max_length
 Usage   : my $length_of_longest_sequence = $obj->max_length($fastqcheck_file);
 Function: Get the length of the longest sequence.
 Returns : int
 Args    : fastqcheck filename

=cut

sub max_length {
    my ($self, $file) = @_;
    return $self->_fastqcheck_getter('max_length', $file);
}

=head2 standard_deviations

 Title   : standard_deviations
 Usage   : my ($total_sd, $per_base_sd) = $obj->standard_deviations($fqc_file);
 Function: Get the standard deviations at 0.25; total and per-base.
 Returns : list: two percentages
 Args    : fastqcheck filename

=cut

sub standard_deviations {
    my ($self, $file) = @_;
    return $self->_fastqcheck_getter('standard_deviations', $file);
}

sub _fastqcheck_getter {
    my ($self, $var, $fastqcheck_file) = @_;
    my $parser = VertRes::Parser::fastqcheck->new(file => $fastqcheck_file);
    return $parser->$var();
}


1;
