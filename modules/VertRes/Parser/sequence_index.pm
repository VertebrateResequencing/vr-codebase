=head1 NAME

VertRes::Parser::sequence_index - parse sequence.index files from the DCC

=head1 SYNOPSIS

use VertRes::Parser::sequence_index;

# create object, supplying sequence.index file or filehandle
my $pars = VertRes::Parser::sequence_index->new(file => 'sequence.index');

# get the array reference that will hold the most recently requested result
my $result_holder = $pars->result_holder();

# loop through the output, getting results
while ($pars->next_result()) {
    my $individual = $result_holder->[0];
    # etc.
}

# or get specific information on a particular run id (lane)
my $individual = $pars->lane_info('SRR005864', 'sample_name');

# get all the lanes for a particular individual
my @lanes = $pars->get_lanes('NA11994');

=head1 DESCRIPTION

A straightforward parser for sequence.index filse.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Parser::sequence_index;

use strict;
use warnings;

use base qw(VertRes::Parser::ParserI);

our %field_to_index = (FASTQ_FILE => 0,
                       MD5 => 1,
                       RUN_ID => 2,
                       STUDY_ID => 3,
                       STUDY_NAME => 4,
                       CENTER_NAME => 5,
                       SUBMISSION_ID => 6,
                       SUBMISSION_DATE => 7,
                       SAMPLE_ID => 8,
                       SAMPLE_NAME => 9,
                       POPULATION => 10,
                       EXPERIMENT_ID => 11,
                       INSTRUMENT_PLATFORM => 12,
                       INSTRUMENT_MODEL => 13,
                       LIBRARY_NAME => 14,
                       RUN_NAME => 15,
                       RUN_BLOCK_NAME => 16,
                       INSERT_SIZE => 17,
                       LIBRARY_LAYOUT => 18,
                       PAIRED_FASTQ => 19,
                       WITHDRAWN => 20,
                       WITHDRAWN_DATE => 21,
                       COMMENT => 22,
                       READ_COUNT => 23,
                       BASE_COUNT => 24,
                       ANALYSIS_GROUP => 25);


=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Parser::sequence_index->new(file => 'filename');
 Function: Build a new VertRes::Parser::sequence_index object.
 Returns : VertRes::Parser::sequence_index object
 Args    : file => filename -or- fh => filehandle

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    return $self;
}

sub _get_header {
    my $self = shift;
    
    my $fh = $self->fh() || return;
    return 1 if $self->_header_parsed();
    
    my $saw = 0;
    while (<$fh>) {
        if (/^FASTQ_FILE\s+MD5\s+RUN_ID\s+STUDY_ID/) {
            $saw++;
            last;
        }
        else {
            # allow header line to not be present. lines can also end with an
            # extraneous \t\n, so remove that first... causes problems in
            # fake sequence.indexes when people leave the last few columns
            # empty...
            #s/\t\n$//;
            my @a = split("\t", $_);
            if (@a == 25 || @a == 26) {
                $self->{_first_line} = $_;
                my $tell = tell($fh);
                if ($tell == -1) {
                    $self->throw("sequence.index has no header line, and you've piped it in - can't cope!");
                }
                $self->seek(0, 0);
                $saw++;
                last;
            }
            else {
                $self->warn("This file has no header line, and has ".scalar(@a)." columns instead of 25 or 26");
                last;
            }
        }
    }
    
    if ($saw) {
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
           [0]  FASTQ_FILE
           [1]  MD5
           [2]  RUN_ID
           [3]  STUDY_ID
           [4]  STUDY_NAME
           [5]  CENTER_NAME
           [6]  SUBMISSION_ID
           [7]  SUBMISSION_DATE
           [8]  SAMPLE_ID
           [9]  SAMPLE_NAME
           [10] POPULATION
           [11] EXPERIMENT_ID
           [12] INSTRUMENT_PLATFORM
           [13] INSTRUMENT_MODEL
           [14] LIBRARY_NAME
           [15] RUN_NAME
           [16] RUN_BLOCK_NAME
           [17] INSERT_SIZE
           [18] LIBRARY_LAYOUT
           [19] PAIRED_FASTQ
           [20] WITHDRAWN
           [21] WITHDRAWN_DATE
           [22] COMMENT
           [23] READ_COUNT
           [24] BASE_COUNT
           [25] ANALYSIS_GROUP
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
    
    # make sure we've gotten our header first
    $self->_get_header() || $self->throw("Unable to parse header before first result - is this really a sequence.index file?");
    
    # get the next line
    my $tell = tell($fh);
    my $line = <$fh>;
    unless ($line) {
        $self->{"saw_last_line_$fh_id"} = 1;
        return;
    }
    
    # bug in sequence.index generation can result in lines ending in \t\n
    # $line =~ s/\t\n//;
    # ... but the above fix prevents us ending the sequence.index in empty
    # values!
    
    my @data = split(/\t/, $line);
    chomp($data[-1]);
    if (@data != 26 && @data != 25) {
        $self->warn("Expected 26 or 25 columns, got " . scalar @data . ": " . join(",", @data) ."\n");
        return;
    }
    
    if ($data[2] eq 'RUN_ID') {
        $self->throw("got RUN_ID!");
    }
    
    $self->{'lanes'.$fh_id}->{$data[2]}->{$tell} = 1;
    
    for my $i (0..$#data) {
        $self->{_result_holder}->[$i] = $data[$i];
    }
    
    return 1;
}

=head2 lane_info

 Title   : lane_info
 Usage   : my $sample_name = $obj->lane_info('SRR005864', 'sample_name');
 Function: Get a particular bit of info for a particular lane (run_id).
 Returns : in scalar context, the most common answer; in list context, all
           the different answers
 Args    : lane/run_id, field name (see result_holder() for the list of valid
           field names)

=cut

sub lane_info {
    my ($self, $lane, $field) = @_;
    $self->throw("lane (aka run_id) must be supplied") unless $lane;
    $self->throw("field must be supplied") unless $field;
    $field = uc($field);
    my $index = $field_to_index{$field};
    $self->throw("'$field' wasn't a valid field name") unless defined $index;
    
    my $fh = $self->fh() || return;
    my $fh_id = $self->_fh_id;
    $self->_save_position || return;
    
    # since a lane can appear multiple times in the file, we have to just
    # parse the whole file
    unless ($self->{"saw_last_line_$fh_id"}) {
        while ($self->next_result) {
            next;
        }
    }
    if (! defined $self->{'lanes'.$fh_id}->{$lane}) {
        $self->warn("'$lane' wasn't a run_id in your sequence.index file");
        return;
    }
    
    my %answers;
    foreach my $tell (keys %{$self->{'lanes'.$fh_id}->{$lane}}) {
        $self->seek($tell, 0);
        $self->next_result;
        $answers{$self->{_result_holder}->[$index]}++;
    }
    
    $self->_restore_position;
    
    if (wantarray) {
        my @answers = sort keys %answers;
        return @answers;
    }
    else {
        my $most_common_answer;
        my $highest_count = 0;
        while (my ($answer, $count) = each %answers) {
            if ($count > $highest_count) {
                $highest_count = $count;
                $most_common_answer = $answer;
            }
        }
        return $most_common_answer;
    }
}

=head2 get_lanes

 Title   : get_lanes
 Usage   : my @lanes = $obj->get_lanes(sample_name => 'NA11994');
 Function: Get all the lane ids with a given property.
 Returns : string
 Args    : none for all lanes, or a single key => value pair to get lanes under
           that key, eg. sample_name => 'NA11994' to get all lanes of that
           sample.

           optionally, you can add one or more key value pairs where keys are
           valid field names (see result_holder() for the list) prefixed with
           'ignore_' and values are what should should be ignored (a regex).
           Eg. to ignore all withdrawn lanes:
           ignore_withdrawn => 1

=cut

sub get_lanes {
    my ($self, %args) = @_;
    
    $self->_save_position || return;
    
    my %ignores;
    my ($field, $value);
    my $do_ignores = 0;
    while (my ($key, $val) = each %args) {
        if ($key =~ /^ignore_(\S+)/i) {
            my $ignore = uc($1);
            my $index = $field_to_index{$ignore};
            $self->throw("'$ignore' wasn't a valid field name") unless defined $index;
            $ignores{$index} = $val;
            $do_ignores = 1;
        }
        else {
            ($field, $value) = ($key, $val);
        }
    }
    
    my $index;
    if ($field && $value) {
        $field = uc($field);
        $index = $field_to_index{$field};
    }
    
    my $result_holder = $self->result_holder();
    my %lanes;
    
    # parse the whole file
    $self->_seek_first_result();
    RESULT: while ($self->next_result) {
        if ($do_ignores) {
            keys %ignores; # reset the iterator, since we may have nexted out
            while (my ($index, $regex) = each %ignores) {
                next RESULT if $result_holder->[$index] =~ /$regex/i;
            }
        }
        
        if (defined $index) {
            if ($result_holder->[$index] eq $value) {
                $lanes{$result_holder->[2]} = 1;
            }
        }
        else {
            $lanes{$result_holder->[2]} = 1;
        }
    }
    
    $self->_restore_position;
    
    return sort keys %lanes;
}

1;
