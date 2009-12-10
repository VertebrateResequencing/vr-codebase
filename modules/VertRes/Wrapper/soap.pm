=head1 NAME

VertRes::Wrapper::soap - wrapper for soap

=head1 SYNOPSIS

[stub]

=head1 DESCRIPTION

37:
    default
54:
    -l 35 -v 1
76:
    -l 35 -v 2
108:
    -l 35 -v 3

2bwt-builder ref.fa
soap -a 1.fq -b 2.fq -D ref.fa.index -o out.txt -m 100
soap2sam.pl out.txt > out.sam ??

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Wrapper::soap;

use strict;
use warnings;
use File::Copy;
use VertRes::IO;

use base qw(VertRes::Wrapper::MapperI);


=head2 new

 Title   : new
 Usage   : my $wrapper = VertRes::Wrapper::soap->new();
 Function: Create a VertRes::Wrapper::soap object.
 Returns : VertRes::Wrapper::soap object
 Args    : quiet   => boolean

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args, exe => '/lustre/scratch102/user/sb10/mapper_comparisons/mappers/soap2.20release/');
    
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
    
    my $index = $ref.'.index';
    my @suffixes = qw(amb ann lkt pac rev.lkt rev.pac);
    my $indexed = 0;
    foreach my $suffix (@suffixes) {
        if (-s "$index.$suffix") {
            $indexed++;
        }
    }
    
    unless ($indexed == @suffixes) {
        my $orig_exe = $self->exe;
        $self->exe($orig_exe.'2bwt-builder');
        $self->simple_run($ref);
        $self->exe($orig_exe);
        
        $indexed = 0;
        foreach my $suffix (@suffixes) {
            if (-s "$index.$suffix") {
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
    
    my $soap_out = $out.'.soap';
    unless (-s $soap_out) {
        my $orig_exe = $self->exe;
        $self->exe($orig_exe.'soap');
        
        # settings depend on read length
        my $longest_read = 0;
        foreach my $fq (@fqs) {
            my $pars = VertRes::Parser::fastqcheck->new(file => "$fq.fastqcheck");
            my $length = $pars->max_length();
            if ($length > $longest_read) {
                $longest_read = $length;
            }
        }
        
        my $extra_args = '';
        if ($longest_read >= 108) {
            $extra_args = ' -l 35 -v 3';
        }
        elsif ($longest_read >= 76) {
            $extra_args = ' -l 35 -v 2';
        }
        elsif ($longest_read >= 54) {
            $extra_args = ' -l 35 -v 1';
        }
        
        $self->simple_run("-a $fqs[0] -b $fqs[1] -D $ref.index -o $soap_out -m 100$extra_args");
        
        $self->exe($orig_exe);
    }
    
    unless (-s $out) {
        my $orig_exe = $self->exe;
        $self->exe($orig_exe.'soap2sam.pl');
        $self->simple_run("$soap_out > $out");
        $self->exe($orig_exe);
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
