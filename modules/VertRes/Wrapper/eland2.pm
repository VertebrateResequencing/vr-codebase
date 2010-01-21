=head1 NAME

VertRes::Wrapper::eland2 - wrapper for eland2

=head1 SYNOPSIS

[stub]

=head1 DESCRIPTION

*** waiting on statically linked version

squashGenome /human_genome_squashed /human_genome/*.fasta
?? buildSeqFromFastq.pl input.fastq batch_size convert
make -j 1 paired CORES=1 INPUT1=<left.fastq> INPUT2=<right.fastq> SAM_OUTPUT=<sam_output.sam> GENOME=<squash_genome_directory>

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Wrapper::eland2;

use strict;
use warnings;
use File::Copy;
use File::Spec;
use File::Basename;
use Cwd qw(abs_path cwd);
use VertRes::IO;
use Bio::SeqIO;

use base qw(VertRes::Wrapper::MapperI);


=head2 new

 Title   : new
 Usage   : my $wrapper = VertRes::Wrapper::eland2->new();
 Function: Create a VertRes::Wrapper::eland2 object.
 Returns : VertRes::Wrapper::eland2 object
 Args    : quiet   => boolean

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args, exe => '/lustre/scratch102/user/sb10/mapper_comparisons/mappers/eland2_standalone_distribution');
    
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
    
    my $squash_dir = $self->_squash_dir($ref);
    
    # 2GB file size limit; split ref by chromosome, but keep all NT/NC sequences
    # in one file to avoid overflow problem when too many ref files
    my $done_splits = File::Spec->catfile($squash_dir, '.done_splits');
    unless (-e $done_splits) {
        my @expected;
        my $i = Bio::SeqIO->new(-file => $ref);
        while (my $s = $i->next_seq) {
            my $id = $s->id;
            if ($id =~ /^N/) {
                $id = 'N';
            }
            my $o_file = File::Spec->catfile($squash_dir, $id.".fa");
            push(@expected, $o_file);
            
            my $o;
            if ($id eq 'N') {
                if (-s $o_file) {
                    my $ni = Bio::SeqIO->new(-file => $o_file);
                    my $done_sequence = 0;
                    while (my $n = $ni->next_seq) {
                        if ($n->id eq $s->id) {
                            $done_sequence = 1;
                            last;
                        }
                    }
                    $ni->close;
                    next if $done_sequence;
                }
                $o = Bio::SeqIO->new(-file => ">>$o_file");
            }
            else {
                next if -s $o_file;
                $o = Bio::SeqIO->new(-file => ">$o_file");
            }
            
            $o->write_seq($s);
            $o->close;
        }
        
        my $exist = 0;
        foreach my $fa (@expected) {
            if (-s $fa) {
                $exist++;
            }
        }
        if ($exist == @expected) {
            open(my $fh, '>', $done_splits);
            foreach my $fa (@expected) {
                print $fh $fa, "\n";
            }
            close($fh);
        }
    }
    
    my $indexed = 0;
    my @expected;
    open(my $fh, '<', $done_splits);
    while (<$fh>) {
        chomp;
        my $fa = $_;
        $fa =~ /fa/ || next;
        foreach my $suffix (qw(2bpb vld)) {
            my $file = $fa.".$suffix";
            push(@expected, $file);
            if (-s $file) {
                $indexed++;
            }
        }
    }
    close($fh);
    
    return 1 if $indexed == @expected;
    
    my $orig_exe = $self->exe;
    $self->exe($orig_exe.'/squashGenome');
    $self->simple_run("$squash_dir $squash_dir/*.fa");
    $self->exe($orig_exe);
    
    $indexed = 0;
    foreach my $file (@expected) {
        if (-s $file) {
            $indexed++;
        }
        else {
            warn "$file is missing\n";
        }
    }
    
    return $indexed == @expected ? 1 : 0;
}

sub _squash_dir {
    my ($self, $ref) = @_;
    my $basename = basename($ref);
    my ($path) = $ref =~ /^(.+)\/$basename$/;
    my $squash_dir = File::Spec->catdir($path, 'ref_squashed');
    mkdir($squash_dir);
    return $squash_dir;
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
        my $orig_exe = $self->exe;
        $self->exe('make');
        my $squash = $self->_squash_dir($ref);
        
        my $fq_dir;
        my $base_dir;
        foreach my $fq (@fqs) {
            $fq_dir .= basename($fq);
            unless ($base_dir) {
                $base_dir = $fq;
                $base_dir =~ s/$fq_dir//;
            }
            $fq =~ s/\.gz$//;
        }
        
        # fqs and the whole eland2 distributions need to be in their own
        # directory or multiple runs will override same-named temp files
        $fq_dir = File::Spec->catdir($base_dir, $fq_dir);
        mkdir($fq_dir);
        my @fq_links;
        foreach my $fq (@fqs) {
            my $dest = File::Spec->catfile($fq_dir, basename($fq));
            symlink($fq, $dest);
            push(@fq_links, $dest);
        }
        
        my $cwd = cwd();
        chdir($fq_dir);
        system("cp -r $orig_exe/* .");
        $self->simple_run("-j 1 paired CORES=1 INPUT1=$fq_links[0] INPUT2=$fq_links[1] SAM_OUTPUT=$out GENOME=$squash");
        $self->exe($orig_exe);
        chdir($cwd);
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
