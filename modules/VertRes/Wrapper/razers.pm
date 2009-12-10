=head1 NAME

VertRes::Wrapper::razers - wrapper for razers

=head1 SYNOPSIS

[stub]

=head1 DESCRIPTION

Author never suggested any settings, so will use defaults.

Example given in readme:
razers example/genome.fa example/reads.fa example/reads2.fa -id -mN -of 4 -o out.sam

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Wrapper::razers;

use strict;
use warnings;
use File::Copy;
use VertRes::IO;
use VertRes::Parser::fastq;

use base qw(VertRes::Wrapper::MapperI);


=head2 new

 Title   : new
 Usage   : my $wrapper = VertRes::Wrapper::razers->new();
 Function: Create a VertRes::Wrapper::razers object.
 Returns : VertRes::Wrapper::razers object
 Args    : quiet => boolean

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args, exe => '/lustre/scratch102/user/sb10/mapper_comparisons/mappers/seqan/projects/library/apps/razers2/razers');
    
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
    
    # convert to fasta
    my $dones = 0;
    foreach my $fq (@fqs) {
        my $fa = $fq.'.fa';
        
        if (-s $fa) {
            $dones++;
        }
        else {
            open(my $fh, '>', $fa);
            my $pars = VertRes::Parser::fastq->new(file => $fq);
            my $rh = $pars->result_holder;
            while ($pars->next_result()) {
                print $fh '>', $rh->[0], "\n", $rh->[1], "\n";
            }
            close($fh);
            
            if (-s $fa) {
                $dones++;
            }
        }
    }
    
    return $dones == @fqs ? 1 : 0;
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
        $self->simple_run("$ref $fqs[0].fa $fqs[1].fa -id -mN -of 4 -o $out");
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
