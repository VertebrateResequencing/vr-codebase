=head1 NAME

VertRes::Parser::dict - parse sequence dictionary files

=head1 SYNOPSIS

use VertRes::Parser::dict;

# create object, supplying .dict file or filehandle
my $pars = VertRes::Parser::dict->new(file => 'reference.dict');

# get the array reference that will hold the most recently requested result
my $result_holder = $pars->result_holder();

# loop through the output, getting results
while ($pars->next_result()) {
    # check $result_holder for desired info, eg:
    my $md5 = $result_holder->[3]
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
    
    return $self;
}

=head2 result_holder

 Title   : result_holder
 Usage   : my $result_holder = $obj->result_holder()
 Function: Get the data structure that will hold the last result requested by
           next_result().
 Returns : array ref, with the values:
           [0] sequence identifier
           [1] sequence length
           [2] file:/... location of the fasta file this a is a .dict of
           [3] md5 checksum of the uppercase-spaces-removed sequence
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
    
    # ignore header and any other non SQ lines
    while (index($line, '@SQ') != 0) {
        $line = <$fh> || return;
    }
    
    chomp($line);
    my @data = split(qr/\t/, $line);
    @data || return;
    
    my $result_holder = $self->{_result_holder};
    foreach my $i (0..3) {
        $result_holder->[$i] = substr($data[$i + 1], 3);
    }
    
    return 1;
}

1;
