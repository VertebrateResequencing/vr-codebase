=head1 NAME

VertRes::Wrapper::bowtie - wrapper for bowtie

=head1 SYNOPSIS

[stub]

=head1 DESCRIPTION

bowtie-build ref.fa out_base
bowtie -e aspermaq -X maxinsert --sam --time index_basename -1 read1.fastq -2 read2.fastq output

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Wrapper::bowtie;

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
 Usage   : my $wrapper = VertRes::Wrapper::bowtie->new();
 Function: Create a VertRes::Wrapper::bowtie object.
 Returns : VertRes::Wrapper::bowtie object
 Args    : quiet   => boolean

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args, exe => '/software/pathogen/external/apps/usr/bin/bowtie');
    
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
    
    #  .1.ebwt, .2.ebwt, .3.ebwt, .4.ebwt, .rev.1.ebwt, and .rev.2.ebwt
    my $existing_ebwts = 0;
    foreach my $suffix (qw(.1.ebwt .2.ebwt .3.ebwt .4.ebwt .rev.1.ebwt .rev.2.ebwt)) {
        my $ebwt = $ref.$suffix;
        $existing_ebwts += -s $ebwt ? 1 : 0;
    }
    
    unless ($existing_ebwts == 6) {
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
        my $e = 70;
        
        my $longest_read = 0;
        foreach my $fq (@fqs) {
            my $pars = VertRes::Parser::fastqcheck->new(file => "$fq.fastqcheck");
            my $length = $pars->max_length();
            if ($length > $longest_read) {
                $longest_read = $length;
            }
        }
        
        foreach (sort { $a <=> $b } keys %e_parameter) {
            if ($longest_read <= $_) {
                $e = $e_parameter{$_};
            }
        }
        
        foreach my $fq (@fqs) {
            $fq =~ s/\.gz$//;
        }
        
        my $X = 1000;
        $self->simple_run("-e $e -X $X --sam --time $ref -1 $fqs[0] -2 $fqs[1] $out");
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
