=head1 NAME

VertRes::Wrapper::bowtie2 - wrapper for bowtie2

=head1 SYNOPSIS

[stub]

=head1 DESCRIPTION

bowtie2-build ref.fa out_base
bowtie2 -e aspermaq -X maxinsert --sam --time index_basename -1 read1.fastq -2 read2.fastq output

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Wrapper::bowtie2;

use strict;
use warnings;
use File::Copy;
use VertRes::IO;
use VertRes::Parser::fastqcheck;

use base qw(VertRes::Wrapper::MapperI);

my %e_parameter = (
    37 	=> 70,
    63	=> 120,
    76  => 140,
    92 	=> 160,
    108 => 160,
    123	=> 200,
    156	=> 240
);

=head2 new

 Title   : new
 Usage   : my $wrapper = VertRes::Wrapper::bowtie2->new();
 Function: Create a VertRes::Wrapper::bowtie2 object.
 Returns : VertRes::Wrapper::bowtie2 object
 Args    : quiet   => boolean

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args, exe => '/software/pathogen/external/apps/usr/bin/bowtie2');
    $self->{orig_exe} = $self->exe;
    
    return $self;
}

=head2 version

 Title   : version
 Usage   : my $version = $obj->version();
 Function: Returns the program version.
 Returns : string representing version of the program 
 Args    : n/a

=cut

sub version {
    my $self = shift;
    
    my $exe = $self->{orig_exe};
    open(my $fh, "$exe --version 2>&1 |") || $self->throw("Could not start $exe");
    my $version = 2;
    while (<$fh>) {
        if (/bowtie2-align\s+version\s+(\S+)/) {
            $version = $1;
            last;
        }
    }
    close($fh);
    
    return $version;
}

=head2 setup_reference

 Title   : setup_reference
 Usage   : $obj->setup_reference($ref_fasta);
 Function: Do whatever needs to be done with the reference to allow mapping.
 Returns : boolean
 Args    : n/a

=cut

sub setup_reference {
    my ($self, $ref) = @_;
    
    #  .1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, and .rev.2.bt2
    my $existing_bt2s = 0;
    foreach my $suffix (qw(.1.bt2 .2.bt2 .3.bt2 .4.bt2 .rev.1.bt2 .rev.2.bt2)) {
        my $bt2 = $ref.$suffix;
        $existing_bt2s += -s $bt2 ? 1 : 0;
    }
    
    unless ($existing_bt2s == 6) {
        my $orig_exe = $self->exe;
        $self->exe($orig_exe.'-build');
        $self->simple_run("$ref $ref");
        $self->exe($orig_exe);
    }
    
    return 1;
}

=head2 setup_fastqs

 Title   : setup_fastqs
 Usage   : $obj->setup_fastqs($ref_fasta, @fastqs);
 Function: Do whatever needs to be done with the fastqs to allow mapping.
 Returns : boolean
 Args    : n/a

=cut

sub setup_fastqs {
    my ($self, $ref, @fqs) = @_;
    
    foreach my $fq (@fqs) {
        if ($fq =~ /\.gz$/) {
            my $fq_new = $fq;
            $fq_new =~ s/\.gz$//;
            
            unless (-s $fq_new) {
                my $i = VertRes::IO->new(file => $fq);
                my $o = VertRes::IO->new(file => ">$fq_new");
                my $ifh = $i->fh;
                my $ofh = $o->fh;
                while (<$ifh>) {
                    print $ofh $_;
                }
                $i->close;
                $o->close;
            }
        }
    }
    
    return 1;
}

=head2 generate_sam

 Title   : generate_sam
 Usage   : $obj->generate_sam($out_sam, $ref_fasta, @fastqs);
 Function: Do whatever needs to be done with the reference and fastqs to
           complete mapping and generate a sam/bam file.
 Returns : boolean
 Args    : n/a

=cut

sub generate_sam {
    my ($self, $out, $ref, @fqs) = @_;
    
    unless (-s $out) {
        foreach my $fq (@fqs) {
            $fq =~ s/\.gz$//;
        }
        
        my $X = 1000;
        $self->simple_run("-X $X --time $ref -1 $fqs[0] -2 $fqs[1] $out");
    }
    
    return -s $out ? 1 : 0;
}

=head2 add_unmapped

 Title   : add_unmapped
 Usage   : $obj->add_unmapped($sam_file, $ref_fasta, @fastqs);
 Function: Do whatever needs to be done with the sam file to add in unmapped
           reads.
 Returns : boolean
 Args    : n/a

=cut

sub add_unmapped {
    my ($self, $sam, $ref, @fqs) = @_;
    return 1;
}

=head2 do_mapping

 Title   : do_mapping
 Usage   : $wrapper->do_mapping(ref => 'ref.fa',
                                read1 => 'reads_1.fastq',
                                read2 => 'reads_2.fastq',
                                output => 'output.sam');
 Function: Run mapper on the supplied files, generating a sam file of the
           mapping. Checks the sam file isn't truncated.
 Returns : n/a
 Args    : required options:
           ref => 'ref.fa'
           output => 'output.sam'

           read1 => 'reads_1.fastq', read2 => 'reads_2.fastq'
           -or-
           read0 => 'reads.fastq'

=cut

=head2 run

 Title   : run
 Usage   : Do not call directly: use one of the other methods instead.
 Function: n/a
 Returns : n/a
 Args    : paths to input/output files

=cut

1;
