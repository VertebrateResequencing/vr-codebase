=head1 NAME

VertRes::Wrapper::mosaik - wrapper for mosaik

=head1 SYNOPSIS

[stub]

=head1 DESCRIPTION

Below are examples from example scripts supplied with release.

MosaikBuild -fr reference/c.elegans_chr2.fasta -oa reference/c.elegans_chr2.dat
# MosaikJump -ia h.sapiens.dat -out h.sapiens_15 -hs 15
MosaikBuild -q fastq/c_elegans_chr2_test.fastq -out sequence_archives/c_elegans_chr2_test.dat -st illumina
MosaikAligner -in sequence_archives/c_elegans_chr2_test.dat -out sequence_archives/c_elegans_chr2_test_aligned.dat -ia reference/c.elegans_chr2.dat -hs 14 -act 17 -mm 2 -m unique
MosaikSort -in sequence_archives/c_elegans_chr2_test_aligned.dat -out sequence_archives/c_elegans_chr2_test_sorted.dat
# MosaikAssembler -in sequence_archives/c_elegans_chr2_test_sorted.dat -out assembly/c.elegans_chr2_test -ia reference/c.elegans_chr2.dat -f ace
MosaikText -in yeast_aligned.dat -sam yeast_aligned.sam

Author suggests MosaikAligner settings of:
37: -mm 4 -act 20 -bw 13 -mhp 100 -ls 100
54: -mm 6 -act 25 -bw 17 -mhp 100 -ls 100
76: -mm 12 -act 35 -bw 29 -mhp 100 -ls 100
108: -mm 15 -act 40 -bw 35 -mhp 100 -ls 100
"-mfl 200" is recommended when using *MosaiBuild*

jump database, which is strongly recommended. MosaikJump is a tool to generate jump database. Please refer Pages 16 and 17 in the attached document. If the jump database is used, the parameters go to
=============================
37: -mm 4 -act 20 -bw 13 -mhp 100 -ls 100 -j jumpFileName
54: -mm 6 -act 25 -bw 17 -mhp 100 -ls 100 -j jumpFileName
76: -mm 12 -act 35 -bw 29 -mhp 100 -ls 100 -j jumpFileName
108: -mm 15 -act 40 -bw 35 -mhp 100 -ls 100 -j jumpFileName

-p 8 for multiprocessors


# latest version setting, dealing with including unmapped reads, Wang Ping used:
MosaikBuild:	-mfl 200 -st illumina
MosaikAligner:	-mm 15 -act 40 -bw 35 -mhp 100 -ls 250
MosaikSort:	-rmm
appendUnmappedRead2Sam fastq1 fastq2 sam updatedSam duplicatedReads


=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Wrapper::mosaik;

use strict;
use warnings;
use File::Copy;
use File::Basename;
use VertRes::IO;
use VertRes::Parser::fastqcheck;

use base qw(VertRes::Wrapper::MapperI);

my $append_exe = '/lustre/scratch102/user/sb10/mapper_comparisons/mappers/AppendUnmappedRead2Sam/appendUnmappedRead2Sam';

=head2 new

 Title   : new
 Usage   : my $wrapper = VertRes::Wrapper::mosaik->new();
 Function: Create a VertRes::Wrapper::mosaik object.
 Returns : VertRes::Wrapper::mosaik object
 Args    : quiet   => boolean

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args, exe => '/lustre/scratch102/user/sb10/mapper_comparisons/mappers/mosaik-source/bin/');
    
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
    
    my $out = $ref.'.dat';
    unless (-s $out) {
        my $orig_exe = $self->exe;
        $self->exe($orig_exe.'MosaikBuild');
        $self->simple_run("-fr $ref -oa $out");
        $self->exe($orig_exe);
    }
    
    if (-s $out) {
        my @jump_suffixes = ('_keys.jmp', '_meta.jmp', '_positions.jmp');
        my $jumped = 1;
        foreach my $suffix (@jump_suffixes) {
            unless (-e "$ref$suffix") {
                $jumped = 0;
                last;
            }
        }
        
        unless ($jumped) {
            my $orig_exe = $self->exe;
            $self->exe($orig_exe.'MosaikJump');
            $self->simple_run("-ia $out -out $ref -hs 15 -mhp 100");
            $self->exe($orig_exe);
            
            $jumped = 1;
            foreach my $suffix (@jump_suffixes) {
                unless (-e "$ref$suffix") {
                    $jumped = 0;
                    last;
                }
            }
        }
        
        return $jumped ? 1 : 0;
    }
    
    return 0;
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
    
    my $out = $self->_fastqs_dat(@fqs);
    
    unless (-s $out) {
        my $orig_exe = $self->exe;
        $self->exe($orig_exe.'MosaikBuild');
        $self->simple_run("-q $fqs[0] -q2 $fqs[1] -out $out -mfl 200 -st illumina");
        $self->exe($orig_exe);
    }
    
    return -s $out ? 1 : 0;
}

sub _fastqs_dat {
    my ($self, @fqs) = @_;
    
    my $merged_fq;
    if (@fqs > 1) {
        foreach my $fq (@fqs) {
            my ($name) = $fq =~ /^(.+)\.fastq\.gz$/;
            
            unless($merged_fq) {
                $merged_fq = $name;
            }
            else {
                $merged_fq .= basename($name);
            }
        }
        
        $merged_fq .= '.fastq.gz';
    }
    else {
        $merged_fq = $fqs[0];
    }
    
    $merged_fq =~ s/\.fastq\.gz$/.dat/;
    
    return $merged_fq;
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
    
    my $orig_exe = $self->exe;
    
    # align
    my $align_out = $out.'.align';
    unless (-s $align_out) {
        # settings change depending on read length
        my $max_length = 0;
        foreach my $fq (@fqs) {
            my $pars = VertRes::Parser::fastqcheck->new(file => "$fq.fastqcheck");
            my $length = $pars->max_length();
            if ($length > $max_length) {
                $max_length = $length;
            }
        }
        my $extra_settings;
        if ($max_length >= 108) {
            $extra_settings = "-mm 15 -act 40 -bw 35 -mhp 100 -ls 100 -j $ref";
        }
        elsif ($max_length >= 76) {
            $extra_settings = "-mm 12 -act 35 -bw 29 -mhp 100 -ls 100 -j $ref";
        }
        elsif ($max_length >= 54) {
            $extra_settings = "-mm 6 -act 25 -bw 17 -mhp 100 -ls 100 -j $ref";
        }
        else {
            $extra_settings = "-mm 4 -act 20 -bw 13 -mhp 100 -ls 100 -j $ref";
        }
        
        $self->exe($orig_exe.'MosaikAligner');
        my $fastqs_dat = $self->_fastqs_dat(@fqs);
        $self->simple_run("-in $fastqs_dat -out $align_out -rur $out.unmapped.fastq -ia $ref.dat -p 8 $extra_settings");
        $self->exe($orig_exe);
    }
    unless (-s $align_out) {
        die "failed during the align step\n";
    }
    
    # sort
    my $sort_out = $out.'.sort';
    unless (-s $sort_out) {
        $self->exe($orig_exe.'MosaikSort');
        $self->simple_run("-in $align_out -out $sort_out -rmm");
        $self->exe($orig_exe);
    }
    unless (-s $sort_out) {
        die "failed during the sort step\n";
    }
    
    # assemble??
    # MosaikAssembler -in sequence_archives/c_elegans_chr2_test_sorted.dat -out assembly/c.elegans_chr2_test -ia reference/c.elegans_chr2.dat -f ace
    
    # convert to sam
    my $sam_out = $out.'.gz';
    unless (-s $sam_out) {
        $self->exe($orig_exe.'MosaikText');
        $self->simple_run("-in $sort_out -sam $out");
        $self->exe($orig_exe);
    }
    
    unless (-s $out) {
        system("gunzip $sam_out");
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
    
    my $all_sam = "$sam.all_reads.sam";
    my $failed = system("$append_exe @fqs $sam $all_sam $sam.duplicated_reads");
    unless ($failed) {
        copy($sam, "$sam.orig");
        move($all_sam, $sam);
        return 1;
    }
    else {
        $self->warn("failed [$append_exe @fqs $sam $all_sam $sam.duplicated_reads]: $!");
        return 0;
    }
    
    #my $unmapped = $sam.'.unmapped.fastq';
    
    # mosaik renames the fastq ids, stripping off the /1 or /2 and doing like:
    # @9:107904153+67M200D9M:R:-161 (mate 2, length=76)
    # but for mapper comparison purposes, it doesn't really matter; the add
    # method still works
    
    #return $self->add_unmapped_from_fastq($sam, $unmapped);
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
