=head1 NAME

VertRes::Wrapper::smalt - wrapper for smalt for mapper comparison

=head1 SYNOPSIS

[stub]

=head1 DESCRIPTION

> smalt index -k 13 -s 4 <hash name> <fasta file>
This produces two files <hash name>.sma and <hash name>.smi
For > 70bp Illumina data of the human genome, '-s 6' should be sufficiently sensitive when generating the hash index. 

-k 13 -s 13 for ~1k sequences


> smalt map -f samsoft -o <out> <hash name> <fastq file A> <fastq file B>


=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Wrapper::smalt;

use strict;
use warnings;
use File::Copy;
use VertRes::IO;

use base qw(VertRes::Wrapper::MapperI);


=head2 new

 Title   : new
 Usage   : my $wrapper = VertRes::Wrapper::smalt->new();
 Function: Create a VertRes::Wrapper::smalt object.
 Returns : VertRes::Wrapper::smalt object
 Args    : quiet       => boolean
           read_length => int (the maximum read_length of your sequences; if
                               fastqcheck files are present, this is auto-set)

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args, exe => 'smalt');
    
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
    
    my $exe = $self->exe;
    open(my $fh, "$exe version 2>&1 |") || $self->throw("Could not start $exe");
    my $version = 0;
    while (<$fh>) {
        if (/Version: (\S+)/) {
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
    
    my @suffixes = qw(small.sma small.smi large.sma large.smi);
    my $indexed = 0;
    foreach my $suffix (@suffixes) {
        if (-s "$ref.$suffix") {
            $indexed++;
        }
    }
    
    unless ($indexed == @suffixes) {
        # we produce multiple sets of hashes, one for <70bp reads, one for >70bp,
        # one for >=100bp and one for >500bp
        $self->simple_run("index -k 13 -s 4 $ref.small $ref");
        $self->simple_run("index -k 13 -s 6 $ref.medium $ref");
        $self->simple_run("index -k 20 -s 13 $ref.large $ref");
        
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
                my $o = VertRes::IO->new(file => ">$fq_new.tmp");
                my $ifh = $i->fh;
                my $ofh = $o->fh;
                my $lines = 0;
                while (<$ifh>) {
                    $lines++;
                    print $ofh $_;
                }
                $i->close;
                $o->close;
                
                # check the decompressed fastq isn't truncated
                $i = VertRes::IO->new(file => "$fq_new.tmp");
                my $actual_lines = $i->num_lines;
                if ($actual_lines == $lines) {
                    move("$fq_new.tmp", $fq_new);
                }
                else {
                    $self->throw("Made $fq_new.tmp, but it only had $actual_lines instead of $lines lines");
                }
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
    my ($self, $out, $ref, $fq1, $fq2, %other_args) = @_;
    
    unless (-s $out) {
        # settings change depending on read length
        my $max_length = 0;
        foreach my $fq ($fq1, $fq2) {
            $fq || next;
            if (-s "$fq.fastqcheck") {
                my $pars = VertRes::Parser::fastqcheck->new(file => "$fq.fastqcheck");
                my $length = $pars->max_length();
                if ($length > $max_length) {
                    $max_length = $length;
                }
            }
            else {
                $max_length = $self->{read_length};
            }
        }
        my $hash_name;
        if ($max_length < 70) {
            $hash_name = $ref.'.small';
        }
        elsif ($max_length >= 100) {
            $hash_name = $ref.'.large';
        }
        else {
            $hash_name = $ref.'.medium';
        }
        
        foreach my $fq ($fq1, $fq2) {
            $fq || next;
            $fq =~ s/\.gz$//;
        }
        
        my $insert_size_arg = '';
        if (defined $other_args{i}) {
            $insert_size_arg = " -i $other_args{i}";
        }
        
        $self->simple_run("map -f samsoft$insert_size_arg -o $out $hash_name $fq1 $fq2");
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
