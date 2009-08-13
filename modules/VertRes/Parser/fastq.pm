=head1 NAME

VertRes::Parser::fastq - parse fastq files

=head1 SYNOPSIS

use VertRes::Parser::fastq;

# create object, supplying fastq file or filehandle
my $pars = VertRes::Parser::fastq->new(file => 'read.fastq.gz');

# get the array reference that will hold the most recently requested result
my $result_holder = $pars->result_holder();

# loop through the output, getting results
while ($pars->next_result()) {
    my $id = $result_holder->[0];
    my $seq_string = $result_holder->[1];
    my $qual_string = $result_holder->[2];
}

# get a list of all the sequence ids
my @ids = $pars->sequence_ids();

# get data for a specific sequence:
my $seq_string = $pars->seq($ids[0]);
my $qual_string = $pars->quality($ids[0]);

# ask if a given id is present in the fastq file:
if ($pars->exists('xyz')) {
    # ...
}

=head1 DESCRIPTION

An indexing parser for fastq files.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Parser::fastq;

use strict;
use warnings;

use base qw(VertRes::Parser::ParserI);

our %type_to_index = (sequence => 1, quality => 2);

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Parser::fastq->new(file => 'filename');
 Function: Build a new VertRes::Parser::fastq object.
 Returns : VertRes::Parser::fastq object
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
           [0] sequence id
           [1] sequence string
           [2] quality string
 Args    : n/a

=cut

=head2 next_result

 Title   : next_result
 Usage   : while ($obj->next_result()) { # look in result_holder }
 Function: Parse the next entry from the fastq file.
 Returns : boolean (false at end of output; check the result_holder for the
           actual result information)
 Args    : n/a

=cut

sub next_result {
    my $self = shift;
    
    # just return if no file set
    my $fh = $self->fh() || return;
    my $fh_id = $self->_fh_id;
    
    # get the next entry (4 lines)
    my $tell = tell($fh);
    my $line = <$fh>;
    unless ($line) {
        $self->{"saw_last_line_$fh_id"} = 1;
        return;
    }
    
    # http://maq.sourceforge.net/fastq.shtml
    # @<seqname>\n<seq>\n+[<seqname>]\n<qual>\n
    # <seqname>	:=	[A-Za-z0-9_.:-]+
    # <seq>	:=	[A-Za-z\n\.~]+
    # <qual>	:=	[!-~\n]+
    
    # seqname
    # for speed purposes, we allow anything for the seqname up to the first
    # space
    unless (index($line, '@') == 0) {
        $self->throw("fastq file was bad, seqname line didn't start as expected: $line");
    }
    my ($seq_name) = split(' ', $line);
    $seq_name = substr($seq_name, 1);
    
    # seq
    # apparently sequence can be split over multiple lines, but we'll take the
    # lazy way out and assume it will always be on one line.
    # For speed, we don't validate the line.
    my $seq = <$fh>;
    $seq || $self->throw("fastq file was truncated - no sequence for $seq_name");
    chomp($seq);
    
    # seqname repeated
    $line = <$fh>;
    $line || $self->throw("fastq file was truncated - no + line $seq_name");
    unless (index($line, '+') == 0) {
        $self->throw("fastq file was bad, + line didn't start as expected: $line");
    }
    
    # qual
    # for speed purposes, we don't validate the line
    my $qual = <$fh>;
    $qual || $self->throw("fastq file was truncated - no quality line for $seq_name");
    chomp($qual);
    
    # we assume that seqnames are unique in our fastqs...
    $self->{'seqnames'.$fh_id}->{$seq_name} = $tell;
    
    $self->{_result_holder}->[0] = $seq_name;
    $self->{_result_holder}->[1] = $seq;
    $self->{_result_holder}->[2] = $qual;
    
    return 1;
}

=head2 sequence_ids

 Title   : sequence_ids
 Usage   : my @ids = $obj->sequence_ids();
 Function: Get all the sequence ids in the fastq file. NB: only works on
           seekable input.
 Returns : list of ids
 Args    : n/a

=cut

sub sequence_ids {
    my $self = shift;
    
    my $fh = $self->fh() || return;
    my $fh_id = $self->_fh_id;
    $self->_save_position || return;
    
    unless ($self->{"saw_last_line_$fh_id"}) {
        while ($self->next_result) {
            next;
        }
    }
    
    $self->_restore_position;
    
    return keys %{$self->{'seqnames'.$fh_id}};
}

=head2 seq

 Title   : seq
 Usage   : my $seq = $obj->seq($id);
 Function: Get the sequence string of a particular sequence.
 Returns : string
 Args    : string (id)

=cut

sub seq {
    my $self = shift;
    return $self->_get_seq_data(@_, 'sequence');
}

=head2 quality

 Title   : quality
 Usage   : my $quality = $obj->quality($id);
 Function: Get the quality string of a particular sequence.
 Returns : string
 Args    : string (id)

=cut

sub quality {
    my $self = shift;
    return $self->_get_seq_data(@_, 'quality');
}

sub _get_seq_data {
    my ($self, $seq_name, $type) = @_;
    
    my $fh = $self->fh() || return;
    my $fh_id = $self->_fh_id;
    
    $self->_parse_to_seqname($seq_name) || return;
    
    $self->_save_position || return;
    
    my $tell = $self->{'seqnames'.$fh_id}->{$seq_name};
    seek($fh, $tell, 0);
    $self->next_result;
    my $index = $type_to_index{$type};
    my $result = $self->result_holder->[$index];
    
    $self->_restore_position;
    
    return $result;
}

sub _parse_to_seqname {
    my ($self, $seq_name, $no_warning) = @_;
    
    my $fh = $self->fh() || return;
    my $fh_id = $self->_fh_id;
    $self->_save_position || return;
    
    if (! defined $self->{'seqnames'.$fh_id}->{$seq_name}) {
        if (! $self->{"saw_last_line_$fh_id"}) {
            while (! defined $self->{'seqnames'.$fh_id}->{$seq_name}) {
                $self->next_result;
                last if $self->{"saw_last_line_$fh_id"};
            }
        }
    }
    if (! defined $self->{'seqnames'.$fh_id}->{$seq_name}) {
        $self->warn("'$seq_name' wasn't found in this fastq file") unless $no_warning;
        return;
    }
    
    $self->_restore_position;
    
    return 1;
}

=head2 exists

 Title   : exists
 Usage   : if ($obj->exists($id)) { ... }
 Function: Find out if a particular sequence name is in this fastq.
 Returns : boolean
 Args    : string (id)

=cut

sub exists {
    my ($self, $seq_name) = @_;
    
    my $fh = $self->fh() || return;
    my $fh_id = $self->_fh_id;
    
    $self->_parse_to_seqname($seq_name, 1);
    
    return exists $self->{'seqnames'.$fh_id}->{$seq_name} ? 1 : 0;
}

1;
