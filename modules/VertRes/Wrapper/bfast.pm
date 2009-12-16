=head1 NAME

VertRes::Wrapper::bfast - wrapper for bfast

=head1 SYNOPSIS

[stub]

=head1 DESCRIPTION

Convert the reads: 
$perl bfast-0.6.1c/scripts/qseq2fastq.pl s <N>

Convert the reference: 
$bfast-0.6.1c/bfast/bfast fasta2brg -f hg18.fa

Create the indexes: 
$bfast-0.6.1c/bfast index -f hg18.fa -m <mask> -w 14 -i <index number>

Search the indexes: 
$bfast-0.6.1c/bfast match -f hg18.fa -r reads.s <N>.fastq > bfast.matches.file.s <N>.bmf

Perform local alignment: 
$bfast-0.6.1c/bfast localalign -f hg18.fa bfast.rg.file.hg18.0.brg -m bfast.matches.file.s <N>.bmf > bfast.aligned.file.s <N>.baf

Filter alignments: 
$bfast-0.6.1c/bfast postprocess -f hg18.fa -i bfast.aligned.file.s <N>.baf -a 3 -O 3 > bfast.reported.file.s <N>.sam 



To reduce RAM consumption at the cost of speed, see the ”-d” option in “bfast index”.
bfast fasta2brg -f <ref.fa> -t
bfast index -f <ref.fa> -m <mask string> -d 0 -w 14 -i <index number> -n <num threads> -T <tmp dir> -t   # since they said 24GB, I'll use -d 1 to split into 4 parts?
bfast match -f <ref.fa> -r <reads.fastq> -n <num threads> -T <tmp dir> -t > <out.bmf>
bfast localalign -f <ref.fa> -m <out.bmf> -n <num threads> -t > <out.baf>
bfast postprocess -f <ref.fa> -i <out.baf> -t > <out.sam>


=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Wrapper::bfast;

use strict;
use warnings;
use File::Basename;
use File::Copy;
use VertRes::IO;
use VertRes::Parser::fastqcheck;

use base qw(VertRes::Wrapper::MapperI);


=head2 new

 Title   : new
 Usage   : my $wrapper = VertRes::Wrapper::bfast->new();
 Function: Create a VertRes::Wrapper::bfast object.
 Returns : VertRes::Wrapper::bfast object
 Args    : quiet   => boolean

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args, exe => '/lustre/scratch102/user/sb10/mapper_comparisons/mappers/bfast-0.6.1c/bfast/bfast');
    
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
    
    unless (-s "$ref.nt.brg") {
        $self->simple_run("fasta2brg -f $ref -t");
    }
    
    # masks 1-10 are for reads <40bp, 11-20 are for >40bp
    my @masks = qw(111111111111111111
                   111010001110001110100011011111
                   11110100110111101010101111
                   11111111111111001111
                   1111011101100101001111111
                   11110111000101010000010101110111
                   1011001101011110100110010010111
                   1110110010100001000101100111001111
                   1111011111111111111
                   11011111100010110111101101
                   1111111111111111111111
                   1111101110111010100101011011111
                   1011110101101001011000011010001111111
                   10111001101001100100111101010001011111
                   11111011011101111011111111
                   111111100101001000101111101110111
                   11110101110010100010101101010111111
                   111101101011011001100000101101001011101
                   1111011010001000110101100101100110100111
                   1111010010110110101110010110111011);
    
    my $all_bifs = 0;
    for my $i (1..20) {
        my $made_bifs = 0;
        for my $j (1..4) {
            my $bif = $ref.".nt.$i.$j.bif";
            if (-s $bif) {
                $made_bifs++;
            }
        }
        if ($made_bifs == 4) {
            $all_bifs += 4;
            next;
        }
        for my $j (1..4) {
            my $bif = $ref.".nt.$i.$j.bif";
            unlink($bif);
        }
        
        my $mask = $masks[$i - 1];
        $self->simple_run("index -f $ref -m $mask -d 1 -w 14 -i $i -n 8 -T /tmp/ -t");
        
        for my $j (1..4) {
            my $bif = $ref.".nt.$i.$j.bif";
            if (-s $bif) {
                $all_bifs++;
            }
        }
    }
    
    return $all_bifs == 80 ? 1 : 0;
}

=head2 setup_fastqs

 Title   : setup_fastqs
 Usage   : $obj->setup_fastqs(@fastqs);
 Function: Do whatever needs to be done with the fastqs to allow mapping.
 Returns : boolean
 Args    : n/a

=cut

sub setup_fastqs {
    my ($self, $ref, @fqs) = @_;
    
    # fastqs must be merged into one file
    my $merged_fq = $self->_merged_fastq_file(@fqs);
    unless (@fqs > 1 && -s $merged_fq) {
        my $o = VertRes::IO->new(file => ">$merged_fq");
        my $ofh = $o->fh;
        foreach my $fq (@fqs) {
            my $i = VertRes::IO->new(file => $fq);
            my $ifh = $i->fh;
            while (<$ifh>) {
                print $ofh $_;
            }
            $i->close;
        }
        $o->close;
    }
    
    my $out = "$merged_fq.bmf";
    return if -s $out;
    
    # settings change depending on read length
    my $max_length = 0;
    foreach my $fq (@fqs) {
        my $pars = VertRes::Parser::fastqcheck->new(file => "$fq.fastqcheck");
        my $length = $pars->max_length();
        if ($length > $max_length) {
            $max_length = $length;
        }
    }
    
    my $i = $max_length >= 40 ? '11-20' : '1-10';
    $self->simple_run("match -f $ref -r $merged_fq -n 8 -i $i -T /tmp/ -t > $out");
    
    return 1;
}

sub _merged_fastq_file {
    my ($self, @fqs) = @_;
    
    my $merged_fq;
    if (@fqs > 1) {
        foreach my $fq (@fqs) {
            my ($name) = $fq =~ /^(.+)\.fastq\.gz$/;
            
            unless($merged_fq) {
                $merged_fq = $name;
            }
            else {
                $merged_fq .= '-'.basename($name);
            }
        }
        
        $merged_fq .= '.fastq';
    }
    else {
        $merged_fq = $fqs[0];
    }
    
    return $merged_fq;
}

=head2 generate_sam

 Title   : generate_sam
 Usage   : $obj->generate_sam($ref_fasta, @fastqs);
 Function: Do whatever needs to be done with the reference and fastqs to
           complete mapping and generate a sam/bam file.
 Returns : boolean
 Args    : n/a

=cut

sub generate_sam {
    my ($self, $out, $ref, @fqs) = @_;
    
    my $merged_fq = $self->_merged_fastq_file(@fqs);
    
    my $baf = "$out.baf";
    unless (-s $baf) {
        $self->simple_run("localalign -f $ref -m $merged_fq.bmf -n 1 -t > $baf");
    }
    
    unless (-s $out) {
        $self->simple_run("postprocess -f $ref -i $baf -t > $out");
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
