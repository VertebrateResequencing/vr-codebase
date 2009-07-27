=head1 NAME

VertRes::Parser::sam - parse the sam format

=head1 SYNOPSIS

use VertRes::Parser::sam;

# create object, supplying sam file or filehandle
my $pars = VertRes::Parser::sam->new(file => 'mapping.sam');

# get the hash reference that will hold the most recently requested result
my $result_holder = $pars->result_holder();

# loop through the output, getting results
while ($pars->next_result()) {
    # check $result_holder for desired info, eg:
    my $flag = $result_holder->{FLAG};
    
    # get info about a flag, eg:
    my $mapped = $pars->is_mapped($flag);
}

=head1 DESCRIPTION

A parser for sam files. Currently ignores header lines.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Parser::sam;

use strict;
use warnings;

use base qw(VertRes::Parser::ParserI);

my %col_to_name = (0  => 'QNAME',
                   1  => 'FLAG',
                   2  => 'RNAME',
                   3  => 'POS',
                   4  => 'MAPQ',
                   5  => 'CIGAR',
                   6  => 'MRNM',
                   7  => 'MPOS',
                   8  => 'ISIZE',
                   9  => 'SEQ',
                   10 => 'QUAL');

my %flags = (paired_tech    => 0x0001,
             paired_map     => 0x0002,
             self_unmapped  => 0x0004,
             mate_unmapped  => 0x0008,
             self_reverse   => 0x0010,
             mate_reverse   => 0x0020,
             '1st_in_pair'  => 0x0040,
             '2nd_in_pair'  => 0x0080,
             not_primary    => 0x0100,
             failed_qc      => 0x0200,
             duplicate      => 0x0400);

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Parser::sam->new(file => 'filename');
 Function: Build a new VertRes::Parser::sam object.
 Returns : VertRes::Parser::sam object
 Args    : file => filename -or- fh => filehandle

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    # unlike normal parsers, our result holder is a hash ref
    $self->{_result_holder} = {};
    
    return $self;
}

=head2 is_sequencing_paired

 Title   : is_sequencing_paired
 Usage   : if ($obj->is_sequencing_paired($flag)) { ... };
 Function: Ask if a given flag indicates the read was paired in sequencing.
 Returns : boolean
 Args    : int (the flag recieved from $result_holder->{FLAG})

=cut

sub is_sequencing_paired {
    my ($self, $flag) = @_;
    return ($flag & $flags{paired_tech}) > 0 ? 1 : 0;
}

=head2 is_mapped_paired

 Title   : is_mapped_paired
 Usage   : if ($obj->is_mapped_paired($flag)) { ... };
 Function: Ask if a given flag indicates the read was mapped in a proper pair.
 Returns : boolean
 Args    : int (the flag recieved from $result_holder->{FLAG})

=cut

sub is_mapped_paired {
    my ($self, $flag) = @_;
    return ($flag & $flags{paired_map}) > 0 ? 1 : 0;
}

=head2 is_mapped

 Title   : is_mapped
 Usage   : if ($obj->is_mapped($flag)) { ... };
 Function: Ask if a given flag indicates the read was itself mapped.
 Returns : boolean
 Args    : int (the flag recieved from $result_holder->{FLAG})

=cut

sub is_mapped {
    my ($self, $flag) = @_;
    return ($flag & $flags{self_unmapped}) == 0 ? 1 : 0;
}

=head2 is_mate_mapped

 Title   : is_mate_mapped
 Usage   : if ($obj->is_mate_mapped($flag)) { ... };
 Function: Ask if a given flag indicates the read's mate was mapped.
 Returns : boolean
 Args    : int (the flag recieved from $result_holder->{FLAG})

=cut

sub is_mate_mapped {
    my ($self, $flag) = @_;
    return ($flag & $flags{mate_unmapped}) == 0 ? 1 : 0;
}

=head2 is_reverse_strand

 Title   : is_reverse_strand
 Usage   : if ($obj->is_reverse_strand($flag)) { ... };
 Function: Ask if a given flag indicates the read is on the reverse stand.
 Returns : boolean
 Args    : int (the flag recieved from $result_holder->{FLAG})

=cut

sub is_reverse_strand {
    my ($self, $flag) = @_;
    return ($flag & $flags{self_reverse}) > 0 ? 1 : 0;
}

=head2 is_mate_reverse_strand

 Title   : is_mate_reverse_strand
 Usage   : if ($obj->is_mate_reverse_strand($flag)) { ... };
 Function: Ask if a given flag indicates the read's mate is on the reverse
           stand.
 Returns : boolean
 Args    : int (the flag recieved from $result_holder->{FLAG})

=cut

sub is_mate_reverse_strand {
    my ($self, $flag) = @_;
    return ($flag & $flags{mate_reverse}) > 0 ? 1 : 0;
}

=head2 is_first

 Title   : is_first
 Usage   : if ($obj->is_first($flag)) { ... };
 Function: Ask if a given flag indicates the read was the first of a pair.
 Returns : boolean
 Args    : int (the flag recieved from $result_holder->{FLAG})

=cut

sub is_first {
    my ($self, $flag) = @_;
    return ($flag & $flags{'1st_in_pair'}) > 0 ? 1 : 0;
}

=head2 is_second

 Title   : is_second
 Usage   : if ($obj->is_second($flag)) { ... };
 Function: Ask if a given flag indicates the read was the second of a pair.
 Returns : boolean
 Args    : int (the flag recieved from $result_holder->{FLAG})

=cut

sub is_second {
    my ($self, $flag) = @_;
    return ($flag & $flags{'2nd_in_pair'}) > 0 ? 1 : 0;
}

=head2 is_primary

 Title   : is_primary
 Usage   : if ($obj->is_primary($flag)) { ... };
 Function: Ask if a given flag indicates the read alignment was primary.
 Returns : boolean
 Args    : int (the flag recieved from $result_holder->{FLAG})

=cut

sub is_primary {
    my ($self, $flag) = @_;
    return ($flag & $flags{not_primary}) == 0 ? 1 : 0;
}

=head2 passes_qc

 Title   : passes_qc
 Usage   : if ($obj->passes_qc($flag)) { ... };
 Function: Ask if a given flag indicates the read passes quality checks.
 Returns : boolean
 Args    : int (the flag recieved from $result_holder->{FLAG})

=cut

sub passes_qc {
    my ($self, $flag) = @_;
    return ($flag & $flags{failed_qc}) == 0 ? 1 : 0;
}

=head2 is_duplicate

 Title   : is_duplicate
 Usage   : if ($obj->is_duplicate($flag)) { ... };
 Function: Ask if a given flag indicates the read was a duplicate.
 Returns : boolean
 Args    : int (the flag recieved from $result_holder->{FLAG})

=cut

sub is_duplicate {
    my ($self, $flag) = @_;
    return ($flag & $flags{duplicate}) > 0 ? 1 : 0;
}

=head2 result_holder

 Title   : result_holder
 Usage   : my $result_holder = $obj->result_holder()
 Function: Get the data structure that will hold the last result requested by
           next_result()
 Returns : hash ref, with the keys:
           QNAME
           FLAG
           RNAME
           POS
           MAPQ
           CIGAR
           MRNM
           MPOS
           ISIZE
           SEQ
           QUAL

           For optional tag fields, the hash ref will also contain corresponding
           keys for those if present, eg. RG. The value will be the value; the
           value type is ignored.
 Args    : n/a

=cut

=head2 next_result

 Title   : next_result
 Usage   : while ($obj->next_result()) { # look in result_holder }
 Function: Parse the next line from the sam file.
 Returns : boolean (false at end of output; check the result_holder for the
           actual result information)
 Args    : n/a

=cut

sub next_result {
    my $self = shift;
    
    # get the next line
    my $fh = $self->fh() || return;
    my $line = <$fh> || return;
    
    # ignore header lines for now
    while (index($line, '@') == 0) {
        $line = <$fh> || return;
    }
    
    my @data = split(qr/\t/, $line);
    @data || return;
    
    # clear data first, incase we don't overwrite all fields
    my $result_holder = $self->{_result_holder};
    foreach my $key (keys %{$result_holder}) {
        delete $result_holder->{$key};
    }
    
    for my $i (0..$#data) {
        chomp($data[$i]) if $i == $#data;
        
        my $name = $col_to_name{$i};
        if ($name) {
            $result_holder->{$name} = $data[$i];
        }
        else {
            ($name, undef, my $value) = split(":", $data[$i]);
            $name || $self->throw("Unable to parse sam line:\n$line");
            $result_holder->{$name} = $value;
        }
    }
    
    return 1;
}

1;
