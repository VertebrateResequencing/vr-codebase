=head1 NAME

VertRes::Wrapper::karma - wrapper for karma

=head1 SYNOPSIS

[stub]

=head1 DESCRIPTION

Needs 24GB, or 16GB with --occurrenceCutoff < about 100

karma ÐcreateIndex Ðreference 1kgref.fa ÐwordSize <readsize/4>
karma Ðreference 1kgref.fa ÐpairedReads _1.fastq.gz _2.fastq.gz

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Wrapper::karma;

use strict;
use warnings;
use File::Copy;
use File::Basename;
use VertRes::IO;

use base qw(VertRes::Wrapper::MapperI);


=head2 new

 Title   : new
 Usage   : my $wrapper = VertRes::Wrapper::karma->new();
 Function: Create a VertRes::Wrapper::karma object.
 Returns : VertRes::Wrapper::karma object
 Args    : quiet   => boolean

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args, exe => '/lustre/scratch102/user/sb10/mapper_comparisons/mappers/karma-0.8.8/karma/karma');
    
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
    
    my $wordSize = 15; # readsize / 4 ... author suggested using 15 for all read sizes
    
    my $no_fa = $ref;
    $no_fa =~ s/\.fa.*$//;
    my @indexes = qw(umfa umwhl umwhr umwihi umwiwp);
    my $indexed  = 0;
    foreach my $suffix (@indexes) {
        my $index = $no_fa.'.'.$suffix;
        if (-s $index || -l $index) {
            $indexed++;
        }
    }
    
    unless ($indexed == @indexes) {
        # new recommendation is to use no wordSize or occurrenceCutoff;
        # this step will need about 72GB of ram
        $self->simple_run("--createIndex --reference $ref");
        
        $indexed  = 0;
        foreach my $suffix (@indexes) {
            my $index = $no_fa.'.'.$suffix;
            if (-s $index) {
                $indexed++;
            }
        }
    }
    
    return $indexed == @indexes ? 1 : 0;
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
    
    unless (-s $out) {
        $self->simple_run("--reference $ref --pairedReads @fqs");
        
        # it generates output files in the current working directory named after
        # the first fastq
        my $fq = basename($fqs[0]);
        $fq =~ s/\.gz//;
        
        my ($out_path) = $out =~ /^(.+)$fq/;
        foreach my $suffix ('R', 'sam', 'stats') {
            move($fq.'.'.$suffix, $out_path);
        }
        
        my $sam = $out_path.$fq.'.sam';
        if (-s $sam) {
            move($sam, $out);
        }
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
