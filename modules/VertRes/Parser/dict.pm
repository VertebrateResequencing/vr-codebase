=head1 NAME

VertRes::Parser::dict - parse sequence dictionary files

=head1 SYNOPSIS

use VertRes::Parser::dict;

# create object, supplying .dict file or filehandle
my $pars = VertRes::Parser::dict->new(file => 'reference.dict');

# get the hash reference that will hold the most recently requested result
my $result_holder = $pars->result_holder();

# loop through the output, getting results
while ($pars->next_result()) {
    # check $result_holder for desired info, eg:
    my $md5 = $result_holder->{M5}
}

=head1 DESCRIPTION

A parser for sequence dictionary files. These are produced by picard-tools
CreateSequenceDictionary.jar and by my sequence_dicter.pl script.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Parser::dict;

use strict;
use warnings;

use base qw(VertRes::Parser::ParserI);

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Parser::dict->new(file => 'filename');
 Function: Build a new VertRes::Parser::dict object.
 Returns : VertRes::Parser::dict object
 Args    : file => filename -or- fh => filehandle

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    # unlike normal parsers, our result holder is a hash ref
    $self->{_result_holder} = {};
    
    return $self;
}

=head2 result_holder

 Title   : result_holder
 Usage   : my $result_holder = $obj->result_holder()
 Function: Get the data structure that will hold the last result requested by
           next_result().
 Returns : hash ref, with the keys:
           SN (sequence identifier)
           LN (sequence length)
           UR (file:/... location of the fasta file this a is a .dict of)
           M5 (md5 checksum of the uppercase-spaces-removed sequence)
           AS (the assembly name)
           SP (the species)
           NB: some may not exist depending upon the completeness of the dict
           file
 Args    : n/a

=cut

=head2 next_result

 Title   : next_result
 Usage   : while ($obj->next_result()) { # look in result_holder }
 Function: Parse the next line from the .dict file.
 Returns : boolean (false at end of output; check the result_holder for the
           actual result information)
 Args    : n/a

=cut

sub next_result {
    my $self = shift;
    
    # get the next line
    my $fh = $self->fh() || return;
    my $line = <$fh> || return;
    
    #@HD     VN:1.0  SO:unsorted
    #@SQ     SN:1    LN:247249719    UR:file:/lustre/scratch103/sanger/team145/g1k/ref/human_b36_male.fa     M5:9ebc6df9496613f373e73396d5b3b6b6
    #@SQ     SN:1    LN:249250621    AS:NCBI37       UR:ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz  M5:1b22b98cdeb4a9304cb5d48026a85128     SP:Human
    
    # ignore header and any other non SQ lines
    while (index($line, '@SQ') != 0) {
        $line = <$fh> || return;
    }
    
    chomp($line);
    my @data = split(qr/\t/, $line);
    @data || return;
    
    my $result_holder = $self->{_result_holder};
    
    foreach my $key (keys %{$result_holder}) {
        delete $result_holder->{$key};
    }
    
    shift(@data);
    foreach my $tag (@data) {
        my ($name, $value) = $tag =~ /^(\S\S):(\S+)/;
        $result_holder->{$name} = $value;
    }
    
    return 1;
}

=head2 seq_lengths

 Title   : seq_lengths
 Usage   : my %seq_lengths = $obj->seq_lengths;
 Function: Get all the sequence lengths.
 Returns : hash (keys as sequence ids, values as int lengths)
 Args    : n/a

=cut

sub seq_lengths {
    my $self = shift;
    
    my %seq_lengths;
    
    $self->_save_position;
    $self->_seek_first_result;
    
    my $rh = $self->result_holder;
    while ($self->next_result) {
        $seq_lengths{$rh->{SN}} = $rh->{LN};
    }
    
    $self->_restore_position();
    return %seq_lengths;
}

=head2 seq_ids

 Title   : seq_ids
 Usage   : my @seq_ids = $obj->seq_ids;
 Function: Get all the sequence ids.
 Returns : list of strings
 Args    : n/a

=cut

sub seq_ids {
    my $self = shift;
    
    my @seq_ids;
    
    $self->_save_position;
    $self->_seek_first_result;
    
    my $rh = $self->result_holder;
    while ($self->next_result) {
        push(@seq_ids, $rh->{SN});
    }
    
    $self->_restore_position();
    return @seq_ids;
}

1;
