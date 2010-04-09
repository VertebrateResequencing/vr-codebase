=head1 NAME

VertRes::Parser::sam - parse the sam format

=head1 SYNOPSIS

use VertRes::Parser::sam;

# create object, supplying sam file or filehandle
my $pars = VertRes::Parser::sam->new(file => 'mapping.sam');

# get header information
my $program = $pars->program();
my %readgroup_info = $pars->readgroup_info();
# etc.

# get the hash reference that will hold the most recently requested result
my $result_holder = $pars->result_holder();

# loop through the output, getting results
while ($pars->next_result()) {
    # check $result_holder for desired info, eg:
    my $flag = $result_holder->{FLAG};
    
    # get info about a flag, eg:
    my $mapped = $pars->is_mapped($flag);
}

# or for speed critical situations:
$pars = VertRes::Parser::sam->new(file => 'mapping.bam');
while (my @fields = $pars->get_fields('QNAME', 'FLAG', 'RG')) {
    # @fields contains the qname, flag and rg tag
}

=head1 DESCRIPTION

A parser for sam and bam files.

The environment variable SAMTOOLS must point to a directory where samtools
source has been compiled, so containing at least bam.h and libbam.a.
See http://cpansearch.perl.org/src/LDS/Bio-SamTools-1.06/README for advice on
gettings things to work. Specifically, you'll probably need to add -fPIC and
-m64 to the CFLAGS line in samtools's Makefile before compiling.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Parser::sam;

use strict;
use warnings;
use Cwd qw(abs_path);
use Inline C => Config => FILTERS => 'Strip_POD' =>
           INC => "-I$ENV{SAMTOOLS}" =>
           LIBS => "-L$ENV{SAMTOOLS} -lbam -lz" =>
           CCFLAGS => '-D_IOLIB=2 -D_FILE_OFFSET_BITS=64';

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

=head2 file

 Title   : file
 Usage   : $obj->file('filename.sam'); # open to read a sam file
           $obj->file('filename.bam'); # open to read a bam file
 Function: Get/set filename; when setting also opens the file and sets fh().
           There is also read support for remote files like
           'ftp://ftp..../file.bam' and it will be downloaded to a temporary
           location and opened.
           There is no support for writing to sam/bam.
 Returns : absolute path of file
 Args    : filename

=cut

sub file {
    my ($self, $filename) = @_;
    
    if ($filename) {
        if ($filename =~ /^ftp:|^http:/) {
            $filename = $self->get_remote_file($filename) || $self->throw("Could not download remote file '$filename'");
        }
        
        # avoid potential problems with caller changing dir and things being
        # relative; also more informative and explicit to throw with full path
        $filename = abs_path($filename);
        
        # set up the open command, handling bam files automatically
        my $open = $filename;
        if ($filename =~ /\.bam$/) {
            $open = "samtools view -h $filename |";
        }
        
        # go ahead and open it (3 arg form not working when middle is optional)
        open(my $fh, $open) || $self->throw("Couldn't open '$open': $!");
        
        $self->{_filename} = $filename;
        $self->fh($fh);
    }
    
    return $self->{_filename};
}

use Inline C => <<'END_C';

=head2 is_sequencing_paired

 Title   : is_sequencing_paired
 Usage   : if ($obj->is_sequencing_paired($flag)) { ... };
 Function: Ask if a given flag indicates the read was paired in sequencing.
 Returns : boolean
 Args    : int (the flag recieved from $result_holder->{FLAG})

=cut

int is_sequencing_paired(SV* self, int flag) {
    return (flag & 0x0001) > 0 ? 1 : 0;
}

=head2 is_mapped_paired

 Title   : is_mapped_paired
 Usage   : if ($obj->is_mapped_paired($flag)) { ... };
 Function: Ask if a given flag indicates the read was mapped in a proper pair.
 Returns : boolean
 Args    : int (the flag recieved from $result_holder->{FLAG})

=cut

int is_mapped_paired(SV* self, int flag) {
    return (flag & 0x0002) > 0 ? 1 : 0;
}

=head2 is_mapped

 Title   : is_mapped
 Usage   : if ($obj->is_mapped($flag)) { ... };
 Function: Ask if a given flag indicates the read was itself mapped.
 Returns : boolean
 Args    : int (the flag recieved from $result_holder->{FLAG})

=cut

int is_mapped(SV* self, int flag) {
    return (flag & 0x0004) == 0 ? 1 : 0;
}

=head2 is_mate_mapped

 Title   : is_mate_mapped
 Usage   : if ($obj->is_mate_mapped($flag)) { ... };
 Function: Ask if a given flag indicates the read's mate was mapped.
 Returns : boolean
 Args    : int (the flag recieved from $result_holder->{FLAG})

=cut

int is_mate_mapped(SV* self, int flag) {
    return (flag & 0x0008) == 0 ? 1 : 0;
}

=head2 is_reverse_strand

 Title   : is_reverse_strand
 Usage   : if ($obj->is_reverse_strand($flag)) { ... };
 Function: Ask if a given flag indicates the read is on the reverse stand.
 Returns : boolean
 Args    : int (the flag recieved from $result_holder->{FLAG})

=cut

int is_reverse_strand(SV* self, int flag) {
    return (flag & 0x0010) > 0 ? 1 : 0;
}

=head2 is_mate_reverse_strand

 Title   : is_mate_reverse_strand
 Usage   : if ($obj->is_mate_reverse_strand($flag)) { ... };
 Function: Ask if a given flag indicates the read's mate is on the reverse
           stand.
 Returns : boolean
 Args    : int (the flag recieved from $result_holder->{FLAG})

=cut

int is_mate_reverse_strand(SV* self, int flag) {
    return (flag & 0x0020) > 0 ? 1 : 0;
}

=head2 is_first

 Title   : is_first
 Usage   : if ($obj->is_first($flag)) { ... };
 Function: Ask if a given flag indicates the read was the first of a pair.
 Returns : boolean
 Args    : int (the flag recieved from $result_holder->{FLAG})

=cut

int is_first(SV* self, int flag) {
    return (flag & 0x0040) > 0 ? 1 : 0;
}

=head2 is_second

 Title   : is_second
 Usage   : if ($obj->is_second($flag)) { ... };
 Function: Ask if a given flag indicates the read was the second of a pair.
 Returns : boolean
 Args    : int (the flag recieved from $result_holder->{FLAG})

=cut

int is_second(SV* self, int flag) {
    return (flag & 0x0080) > 0 ? 1 : 0;
}

=head2 is_primary

 Title   : is_primary
 Usage   : if ($obj->is_primary($flag)) { ... };
 Function: Ask if a given flag indicates the read alignment was primary.
 Returns : boolean
 Args    : int (the flag recieved from $result_holder->{FLAG})

=cut

int is_primary(SV* self, int flag) {
    return (flag & 0x0100) == 0 ? 1 : 0;
}

=head2 passes_qc

 Title   : passes_qc
 Usage   : if ($obj->passes_qc($flag)) { ... };
 Function: Ask if a given flag indicates the read passes quality checks.
 Returns : boolean
 Args    : int (the flag recieved from $result_holder->{FLAG})

=cut

int passes_qc(SV* self, int flag) {
    return (flag & 0x0200) == 0 ? 1 : 0;
}

=head2 is_duplicate

 Title   : is_duplicate
 Usage   : if ($obj->is_duplicate($flag)) { ... };
 Function: Ask if a given flag indicates the read was a duplicate.
 Returns : boolean
 Args    : int (the flag recieved from $result_holder->{FLAG})

=cut

int is_duplicate(SV* self, int flag) {
    return (flag & 0x0400) > 0 ? 1 : 0;
}

END_C

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

=head2 program_info

 Title   : program_info
 Usage   : my %all_program_info = $obj->program_info();
 Function: Get information about the programs used to create/process this bam,
           as reported in the header.
 Returns : undef if no PG lines in header, else:
           with no args: hash (keys are program ids, values are hash refs with
                               keys as tags (like VN and CL))
           with just a program id: hash (keys as tags, like VN and CL)
           with a program and a tag: the value of that tag for that program
 Args    : none for all info,
           program id for all the info for just that program,
           program id and tag (like 'VN' or 'CL') for specific info

=cut

sub program_info {
    my $self = shift;
    return $self->_handle_multi_line_header_types('PG', @_);
}

=head2 program

 Title   : program
 Usage   : my $program = $obj->program();
 Function: Return the program used to do the mapping, as given in the header.
           If there is more than 1 PG header line, tries to guess which one is
           for the mapping program.
 Returns : string (undef if no header or not given in header)
 Args    : n/a

=cut

sub program {
    my $self = shift;
    return $self->_guess_mapping_program();
}

sub _guess_mapping_program {
    my $self = shift;
    
    my %info = $self->program_info();
    my @programs = keys %info;
    
    if (@programs == 1) {
        return $programs[0];
    }
    else {
        foreach my $program (@programs) {
            if ($program =~ /bwa|maq|ssha/) {
                return $program;
            }
        }
        
        # guess randomly
        return $programs[0];
    }
}

=head2 program_version

 Title   : program_version
 Usage   : my $program_version = $obj->program_version();
 Function: Return the program version used to do the mapping, as given in the
           header.
           If there is more than 1 PG header line, tries to guess which one is
           for the mapping program.
 Returns : string (undef if no header or not given in header)
 Args    : n/a

=cut

sub program_version {
    my $self = shift;
    my $program_id = $self->_guess_mapping_program();
    return $self->program_info($program_id, 'VN');
}

=head2 command_line

 Title   : command_line
 Usage   : my $command_line = $obj->command_line();
 Function: Return the command line used to do the mapping, as given in the
           header.
           If there is more than 1 PG header line, tries to guess which one is
           for the mapping program.
 Returns : string (undef if no header or not given in header)
 Args    : n/a

=cut

sub command_line {
    my $self = shift;
    my $program_id = $self->_guess_mapping_program();
    return $self->program_info($program_id, 'CL');
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
           my %rg_info = $obj->readgroup_info('SRR00001');
           my $library = $obj->readgroup_info('SRR00001', 'LB');
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
            
            if ($type eq 'HD') {
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
    $self->{'_first_record'.$fh_id} = $non_header;
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
    my $fh_id = $self->_fh_id;
    
    # make sure we've gotten our header first
    $self->_get_header();
    my $line;
    if ($self->{'_first_record'.$fh_id}) {
        $line = delete $self->{'_first_record'.$fh_id};
    }
    else {
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

=head2 get_fields

 Title   : get_fields
 Usage   : while (my @fields = $obj->get_fields('QNAME', 'FLAG', 'RG')) { #... }
 Function: From the next line in the bam file, get the values of certain
           fields/tags. This is much faster than using next_result().
           NB: this is incompatible with other non-flag methods, and only works
           on bam files (not sam files, or opened filehandles).
 Returns : list of desired values (if a desired tag isn't present, '*' will be
           returned in that slot)
 Args    : list of desired fields (see result_holder()) or tags (like 'RG').
           additionaly, there are the psuedo-fields 'SEQ_LENGTH' to get the
           raw length of the read (including hard/soft clipped bases) and
           'MAPPED_SEQ_LENGTH' (only bases that match or mismatch to the
           reference, ie. cigar operator M).

=cut

sub get_fields {
    my ($self, @fields) = @_;
    
    my $fh_id = $self->_fh_id || return;
    unless (defined $self->{"_getfields_$fh_id"}) {
        ($self->{_chead}, $self->{_cbam}, $self->{_cb}) = $self->_initialize_bam($self->file());
        $self->{"_getfields_$fh_id"} = 1;
    }
    
    return $self->_get_fields($self->{_cbam}, $self->{_cb}, $self->{_chead}, @fields);
}

use Inline C => <<'END_C';

#include "bam.h"

void _initialize_bam(SV* self, char* bamfile) {
    bamFile *bam;
    bam = bam_open(bamfile, "r");
    
    bam1_t *b;
    b = bam_init1();
    
    bam_header_t *bh;
    bgzf_seek(bam,0,0);
    bh = bam_header_read(bam);
    
    Inline_Stack_Vars;
    Inline_Stack_Reset;
    Inline_Stack_Push(newRV_noinc(newSViv(bh)));
    Inline_Stack_Push(newRV_noinc(newSViv(bam)));
    Inline_Stack_Push(newRV_noinc(newSViv(b)));
    Inline_Stack_Done;
}

void _get_fields(SV* self, SV* bam_ref, SV* b_ref, SV* header_ref, ...) {
    bamFile *bam;
    bam = (bamFile*)SvIV(SvRV(bam_ref));
    bam1_t *b;
    b = (bam1_t*)SvIV(SvRV(b_ref));
    
    int i;
    char *field;
    STRLEN field_length;
    uint8_t *tag_value;
    int type;
    
    int32_t tid;
    bam_header_t *header;
    uint32_t  *cigar;
    int cigar_loop;
    AV *cigar_avref;
    char *cigar_str;
    char *cigar_digits;
    int cigar_digits_length;
    int cigar_digits_i;
    int cigar_chars_total;
    
    char *cigar_op;
    int cigar_op_length;
    int raw_seq_length;
    int mapped_seq_length;
    
    char *seq;
    int seq_i;
    uint8_t *qual;
    int qual_i;
    char *qual_str;
    
    Inline_Stack_Vars;
    
    Inline_Stack_Reset;
    if (bam_read1(bam, b) >= 0) {
        for (i = 2; i < Inline_Stack_Items; i++) {
            field = SvPV(Inline_Stack_Item(i), field_length);
            
            if (field_length > 2) {
                if (strEQ(field, "QNAME")) {
                    Inline_Stack_Push(sv_2mortal(newSVpv(bam1_qname(b), 0)));
                }
                else if (strEQ(field, "FLAG")) {
                    Inline_Stack_Push(sv_2mortal(newSVuv(b->core.flag)));
                }
                else if (strEQ(field, "RNAME")) {
                    if (b->core.tid < 0) {
                        Inline_Stack_Push(sv_2mortal(newSVpv("*", 1)));
                    }
                    else {
                        header = (bam_header_t*)SvIV(SvRV(header_ref));
                        Inline_Stack_Push(sv_2mortal(newSVpv(header->target_name[b->core.tid], 0)));
                    }
                }
                else if (strEQ(field, "POS")) {
                    Inline_Stack_Push(sv_2mortal(newSVuv(b->core.pos + 1)));
                }
                else if (strEQ(field, "MAPQ")) {
                    Inline_Stack_Push(sv_2mortal(newSVuv(b->core.qual)));
                }
                else if (strEQ(field, "CIGAR")) {
                    if (b->core.n_cigar == 0) {
                        Inline_Stack_Push(sv_2mortal(newSVpv("*", 1)));
                    }
                    else {
                        cigar = bam1_cigar(b);
                        cigar_str = Newxz(cigar_str, b->core.n_cigar * 5, char);
                        cigar_chars_total = 0;
                        cigar_digits = Newxz(cigar_digits, 3, char);
                        for (cigar_loop = 0; cigar_loop < b->core.n_cigar; ++cigar_loop) {
                            Renew(cigar_digits, 3, char);
                            cigar_digits_length = sprintf(cigar_digits, "%i", cigar[cigar_loop]>>BAM_CIGAR_SHIFT);
                            for (cigar_digits_i = 0; cigar_digits_i < cigar_digits_length; ++cigar_digits_i) {
                                cigar_str[cigar_chars_total] = cigar_digits[cigar_digits_i];
                                cigar_chars_total++;
                            }
                            
                            cigar_str[cigar_chars_total] =  "MIDNSHP"[cigar[cigar_loop]&BAM_CIGAR_MASK];
                            cigar_chars_total++;
                        }
                        
                        Inline_Stack_Push(sv_2mortal(newSVpv(cigar_str, cigar_chars_total)));
                        
                        Safefree(cigar_str);
                        Safefree(cigar_digits);
                    }
                }
                else if (strEQ(field, "MRNM")) {
                    if (b->core.mtid < 0) {
                        Inline_Stack_Push(sv_2mortal(newSVpv("*", 1)));
                    }
                    else {
                        header = (bam_header_t*)SvIV(SvRV(header_ref));
                        Inline_Stack_Push(sv_2mortal(newSVpv(header->target_name[b->core.mtid], 0)));
                    }
                }
                else if (strEQ(field, "MPOS")) {
                    Inline_Stack_Push(sv_2mortal(newSVuv(b->core.mpos + 1)));
                }
                else if (strEQ(field, "ISIZE")) {
                    Inline_Stack_Push(sv_2mortal(newSViv((int*)b->core.isize)));
                }
                else if (strEQ(field, "SEQ_LENGTH") || strEQ(field, "MAPPED_SEQ_LENGTH")) {
                    if (b->core.n_cigar == 0) {
                        if (b->core.l_qseq) {
                            Inline_Stack_Push(sv_2mortal(newSVuv(b->core.l_qseq)));
                        }
                        else {
                            Inline_Stack_Push(sv_2mortal(newSVuv(0)));
                        }
                    }
                    else {
                        cigar = bam1_cigar(b);
                        raw_seq_length = 0;
                        mapped_seq_length = 0;
                        for (cigar_loop = 0; cigar_loop < b->core.n_cigar; ++cigar_loop) {
                            cigar_op_length = cigar[cigar_loop]>>BAM_CIGAR_SHIFT;
                            cigar_op = "MIDNSHP"[cigar[cigar_loop]&BAM_CIGAR_MASK];
                            
                            if (cigar_op == 'S' || cigar_op == 'H' || cigar_op == 'I') {
                                raw_seq_length = raw_seq_length + cigar_op_length;
                            }
                            else if (cigar_op == 'M') {
                                raw_seq_length = raw_seq_length + cigar_op_length;
                                mapped_seq_length = mapped_seq_length + cigar_op_length;
                            }
                        }
                        
                        if (strEQ(field, "SEQ_LENGTH")) {
                            Inline_Stack_Push(sv_2mortal(newSVuv(raw_seq_length)));
                        }
                        else {
                            Inline_Stack_Push(sv_2mortal(newSVuv(mapped_seq_length)));
                        }
                    }
                }
                else if (strEQ(field, "SEQ")) {
                    if (b->core.l_qseq) {
                        seq = Newxz(seq, b->core.l_qseq + 1, char);
                        for (seq_i = 0; seq_i < b->core.l_qseq; ++seq_i) {
                            seq[seq_i] = bam_nt16_rev_table[bam1_seqi(bam1_seq(b), seq_i)];
                        }
                        Inline_Stack_Push(sv_2mortal(newSVpv(seq, b->core.l_qseq)));
                        Safefree(seq);
                    }
                    else {
                        Inline_Stack_Push(sv_2mortal(newSVpv("*", 1)));
                    }
                }
                else if (strEQ(field, "QUAL")) {
                    if (b->core.l_qseq) {
                        qual = bam1_qual(b);
                        if (qual[0] != 0xff) {
                            qual_str = Newxz(qual_str, b->core.l_qseq + 1, char);
                            for (qual_i = 0; qual_i < b->core.l_qseq; ++qual_i) {
                                qual_str[qual_i] = qual[qual_i] + 33;
                            }
                            Inline_Stack_Push(sv_2mortal(newSVpv(qual_str, b->core.l_qseq)));
                            Safefree(qual_str);
                        }
                        else {
                            Inline_Stack_Push(sv_2mortal(newSVpv("*", 1)));
                        }
                        
                    }
                    else {
                        Inline_Stack_Push(sv_2mortal(newSVpv("*", 1)));
                    }
                }
            }
            else {
                tag_value = bam_aux_get(b, field);
                
                if (tag_value != 0) {
                    type = *tag_value++;
                    switch (type) {
                        case 'c':
                            Inline_Stack_Push(sv_2mortal(newSViv((int32_t)*(int8_t*)tag_value)));
                            break;
                        case 'C':
                            Inline_Stack_Push(sv_2mortal(newSViv((int32_t)*(uint8_t*)tag_value)));
                            break;
                        case 's':
                            Inline_Stack_Push(sv_2mortal(newSViv((int32_t)*(int16_t*)tag_value)));
                            break;
                        case 'S':
                            Inline_Stack_Push(sv_2mortal(newSViv((int32_t)*(uint16_t*)tag_value)));
                            break;
                        case 'i':
                            Inline_Stack_Push(sv_2mortal(newSViv(*(int32_t*)tag_value)));
                            break;
                        case 'I':
                            Inline_Stack_Push(sv_2mortal(newSViv((int32_t)*(uint32_t*)tag_value)));
                            break;
                        case 'f':
                            Inline_Stack_Push(sv_2mortal(newSVnv(*(float*)tag_value)));
                            break;
                        case 'A':
                            Inline_Stack_Push(sv_2mortal(newSVpv((char*)tag_value, 1)));
                            break;
                        case 'Z':
                        case 'H':
                            Inline_Stack_Push(sv_2mortal(newSVpv((char*)tag_value, 0)));
                            break;
                    }
                }
                else {
                    Inline_Stack_Push(sv_2mortal(newSVpv("*", 1)));
                }
            }
        }
    }
    Inline_Stack_Done;
}

END_C

1;
