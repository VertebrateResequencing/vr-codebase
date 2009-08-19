=head1 NAME

VertRes::Parser::sam - parse the sam format

=head1 SYNOPSIS

use VertRes::Parser::sam;

# create object, supplying sam file or filehandle
my $pars = VertRes::Parser::sam->new(file => 'mapping.sam');

# get header information


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

A parser for sam files.

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

our %flags = (paired_tech    => 0x0001,
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

=head2 sam_version

 Title   : sam_version
 Usage   : my $sam_version = $obj->sam_version();
 Function: Return the file format version of this sam file, as given in the
           header.
 Returns : number (undef if no header)
 Args    : n/a

=cut

sub sam_version {
    my $self = shift;
    return $self->_get_single_header_tag('HD', 'VN');
}

=head2 group_order

 Title   : group_order
 Usage   : my $group_order = $obj->group_order();
 Function: Return the group order of this sam file, as given in the header.
 Returns : string (undef if no header or not given in header)
 Args    : n/a

=cut

sub group_order {
    my $self = shift;
    return $self->_get_single_header_tag('HD', 'GO');
}

=head2 sort_order

 Title   : sort_order
 Usage   : my $sort_order = $obj->sort_order();
 Function: Return the sort order of this sam file, as given in the header.
 Returns : string (undef if no header or not given in header)
 Args    : n/a

=cut

sub sort_order {
    my $self = shift;
    return $self->_get_single_header_tag('HD', 'SO');
}

=head2 program

 Title   : program
 Usage   : my $program = $obj->program();
 Function: Return the program used to do the mapping, as given in the header.
 Returns : string (undef if no header or not given in header)
 Args    : n/a

=cut

sub program {
    my $self = shift;
    return $self->_get_single_header_tag('PG', 'ID');
}

=head2 program_version

 Title   : program_version
 Usage   : my $program_version = $obj->program_version();
 Function: Return the program version used to do the mapping, as given in the
           header.
 Returns : string (undef if no header or not given in header)
 Args    : n/a

=cut

sub program_version {
    my $self = shift;
    return $self->_get_single_header_tag('PG', 'VN');
}

=head2 command_line

 Title   : command_line
 Usage   : my $command_line = $obj->command_line();
 Function: Return the command line used to do the mapping, as given in the
           header.
 Returns : string (undef if no header or not given in header)
 Args    : n/a

=cut

sub command_line {
    my $self = shift;
    return $self->_get_single_header_tag('PG', 'CL');
}

=head2 sequence_info

 Title   : sequence_info
 Usage   : my %all_sequences_info = $obj->sequence_info();
           my %sequence_info = $obj->sequence_info('chr1');
           my $seq_length = $obj->sequence_info('chr1', 'LN');
 Function: Get information about the reference sequences, as reported in the
           header.
 Returns : undef if no SQ lines in header, else:
           with no args: hash (keys are sequence ids, values are hash refs with
                               keys as tags (like LN and M5))
           with just a sequence id: hash (keys as tags, like LN and M5)
           with a sequence and a tag: the value of that tag for that sequence
 Args    : none for all info,
           sequence id for all the info for just that sequence,
           sequence id and tag (like 'LN' or 'M5') for specific info

=cut

sub sequence_info {
    my $self = shift;
    return $self->_handle_multi_line_header_types('SQ', @_);
}

=head2 readgroup_info

 Title   : readgroup_info
 Usage   : my %all_rg_info = $obj->readgroup_info();
           my %rg_info = $obj->sequence_info('SRR00001');
           my $library = $obj->sequence_info('SRR00001', 'LB');
 Function: Get information about the read groups, as reported in the header.
 Returns : undef if no RG lines in header, else:
           with no args: hash (keys are sequence ids, values are hash refs with
                               keys as tags (like LB and SM))
           with just a readgroup id: hash (keys as tags, like LB and SM)
           with a readgroup and a tag: the value of that tag for that readgroup
 Args    : none for all info,
           readgroup id for all the info for just that readgroup,
           readgroup id and tag (like 'LB' or 'SM') for specific info

=cut

sub readgroup_info {
    my $self = shift;
    return $self->_handle_multi_line_header_types('RG', @_);
}

sub _handle_multi_line_header_types {
    my ($self, $type, $id, $tag) = @_;
    
    my $lines = $self->_get_header_type($type) || return;
    
    # organise the data into by-id hash
    my %all_info;
    foreach my $line (@{$lines}) {
        my %this_data = $self->_tags_to_hash(@{$line});
        my $this_id = $this_data{SN} || $this_data{ID};
        delete $this_data{SN};
        delete $this_data{ID};
        
        $all_info{$this_id} = \%this_data;
    }
    
    if ($id) {
        my $id_info = $all_info{$id} || return;
        if ($tag) {
            return $id_info->{$tag};
        }
        else {
            return %{$id_info};
        }
    }
    else {
        return %all_info;
    }
}

sub _get_single_header_tag {
    my ($self, $type, $tag) = @_;
    
    my $type_data = $self->_get_header_type($type) || return;
    
    my %data = $self->_tags_to_hash(@{$type_data});
    
    return $data{$tag};
}

sub _tags_to_hash {
    my ($self, @tags) = @_;
    
    my %hash;
    foreach my $tag (@tags) {
        my ($this_tag, $value) = $tag =~ /^(\w\w):(.+)/;
        $hash{$this_tag} = $value;
    }
    return %hash;
}

sub _get_header_type {
    my ($self, $type) = @_;
    
    my $fh = $self->fh() || return;
    my $fh_id = $self->_fh_id;
    
    $self->_get_header();
    
    if (defined $self->{'_header'.$fh_id} && defined $self->{'_header'.$fh_id}->{$type}) {
        return $self->{'_header'.$fh_id}->{$type};
    }
    
    return;
}

sub _get_header {
    my $self = shift;
    
    my $fh = $self->fh() || return;
    my $fh_id = $self->_fh_id;
    
    return if $self->{'_got_header'.$fh_id};
    
    my $non_header;
    while (<$fh>) {
        if (/^@/) {
            #@HD     VN:1.0  GO:none SO:coordinate
            #@SQ     SN:1    LN:247249719    AS:NCBI36       UR:file:/nfs/sf8/G1K/ref/human_b36_female.fa    M5:28f4ff5cf14f5931d0d531a901236378
            #@RG     ID:SRR003447    PL:ILLUMINA     PU:BI.PE1.080723_SL-XBH_0003_FC3044EAAXX.7    LB:Solexa-5453    PI:500  SM:NA11918      CN:BI
            #@PG     ID:xxxx    VN:xxx  CL:xxx
            my @tags = split("\t", $_);
            my $type = shift @tags;
            $type = substr($type, 1);
            
            if ($type eq 'HD' || $type eq 'PG') {
                # we only expect and handle one of these lines per file
                $self->{'_header'.$fh_id}->{$type} = \@tags;
            }
            else {
                push(@{$self->{'_header'.$fh_id}->{$type}}, \@tags);
            }
        }
        else {
            # allow header line to not be present
            $non_header = $_;
            last;
        }
    }
    
    $self->{'_got_header'.$fh_id} = 1;
    return $non_header;
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
    
    # make sure we've gotten our header first
    my $line = $self->_get_header();
    $line ||= <$fh> || return;
    
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
