=head1 NAME

VertRes::Wrapper::srprism - wrapper for srprism

=head1 SYNOPSIS

[stub]

=head1 DESCRIPTION

# supposed to run in under 8GB, but using more than 10?

srprism mkindex -i reference.fa -o reference -M 7168 
srprism search --trace-level info -I reference -i "1.fq,2.fq" -F fastq \
-n 2 -R 0 -r 2 -M 7168 -p true -s 300 -f 300 -o out -O sam 

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Wrapper::srprism;

use strict;
use warnings;
use File::Copy;
use VertRes::IO;

use base qw(VertRes::Wrapper::MapperI);


=head2 new

 Title   : new
 Usage   : my $wrapper = VertRes::Wrapper::srprism->new();
 Function: Create a VertRes::Wrapper::srprism object.
 Returns : VertRes::Wrapper::srprism object
 Args    : quiet   => boolean

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args, exe => '/lustre/scratch102/user/sb10/mapper_comparisons/mappers/srprism/app/srprism');
    
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
    
    my @suffixes = qw(idx imp map pmp ss);
    my $indexed = 0;
    foreach my $suffix (@suffixes) {
        if (-s "$ref.$suffix") {
            $indexed++;
        }
    }
    
    unless ($indexed == @suffixes) {
        $self->simple_run("mkindex -i $ref -o $ref -M 7168");
        
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
        
        $self->simple_run("search -i \"$fqs[0],$fqs[1]\" -F fastq -I $ref --trace-level info -o $out -O sam -n 2 -R 0 -r 2 -M 7168 -p true -s 300 -f 300");
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
    
    # could filter out non unique best hits?
    # gawk '{if($5==100){print}}' the.sam > the.sam.uniq
    
    return $self->SUPER::add_unmapped(@_);
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
