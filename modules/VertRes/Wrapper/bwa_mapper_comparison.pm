=head1 NAME

VertRes::Wrapper::bwa - wrapper for bwa for mapper comparison

=head1 SYNOPSIS

[stub]

=head1 DESCRIPTION

aln -q 15

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Wrapper::bwa_mapper_comparison;

use strict;
use warnings;
use File::Copy;
use VertRes::IO;

use base qw(VertRes::Wrapper::MapperI);


=head2 new

 Title   : new
 Usage   : my $wrapper = VertRes::Wrapper::bwa_mapper_comparison->new();
 Function: Create a VertRes::Wrapper::bwa_mapper_comparison object.
 Returns : VertRes::Wrapper::srprism object
 Args    : quiet   => boolean

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args, exe => 'bwa');
    
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
    return 0;
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
    
    my @suffixes = qw(amb ann bwt pac rbwt rpac rsa sa);
    my $indexed = 0;
    foreach my $suffix (@suffixes) {
        if (-s "$ref.$suffix") {
            $indexed++;
        }
    }
    
    unless ($indexed == @suffixes) {
        $self->simple_run("index -a bwtsw $ref");
        
        $indexed = 0;
        foreach my $suffix (@suffixes) {
            if (-s "$ref.$suffix") {
                $indexed++;
            }
        }
    }
    
    return $indexed == @suffixes ? 1 : 0;
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
        my $sai = $fq.'.sai';
        unless (-s $sai) {
            $self->simple_run("aln -q 15 $ref $fq > $sai");
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
        my @sais;
        foreach my $fq (@fqs) {
            my $sai = $fq.'.sai';
            push(@sais, $sai);
        }
        
        $self->simple_run("sampe $ref @sais @fqs > $out");
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
    my $self = shift;
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
