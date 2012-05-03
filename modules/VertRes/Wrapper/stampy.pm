=head1 NAME

VertRes::Wrapper::stampy - wrapper for stampy

=head1 SYNOPSIS

[stub]

=head1 DESCRIPTION

stampy.py -G hg18 /data/genomes/hg18/*.fa.gz (hg18.stidx)
stampy.py -g hg18 -H hg18 (hg18.sthash)
stampy.py --bwaoptions="-q15 /lustre/scratch102/user/sb10/mapper_comparisons/test_runs/bwa/ref.fa" -g hg18 -h hg18 -M solexareads_1.fastq,solexareads_2.fastq

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Wrapper::stampy;

use strict;
use warnings;
use File::Copy;
use VertRes::IO;
use VertRes::Utils::FileSystem;

use base qw(VertRes::Wrapper::MapperI);


=head2 new

 Title   : new
 Usage   : my $wrapper = VertRes::Wrapper::stampy->new();
 Function: Create a VertRes::Wrapper::stampy object.
 Returns : VertRes::Wrapper::stampy object
 Args    : quiet => boolean

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(exe => 'stampy', @args);
    
    $self->{orig_exe} = $self->exe;
    
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
    
    my $exe = $self->{orig_exe};
    open(my $fh, "$exe |");
    my $line = <$fh>;
    while (<$fh>) {
        next;
    }
    close($fh);
    
    my @nums = $line =~ /v(\d+)\.(\d+).+?(\d+)/;
    
    return join('.', @nums);
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
    
    my $stidx = $ref.'.stidx';
    unless (-s $stidx) {
        $self->simple_run("-G $ref $ref");
    }
    
    my $sthash = $ref.'.sthash';
    unless (-s $sthash) {
        $self->simple_run("-g $ref -H $ref");
    }
    
    return (-s $stidx && -s $sthash) ? 1 : 0;
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
    my ($self, $out, $ref, $fq1, $fq2, %other_args) = @_;
    my @fqs = ($fq1);
    if ($fq2) {
        push(@fqs, $fq2);
    }
    my $fqs_string = join(',', @fqs);
    
    unless (-s $out) {
        # reference files must be copied to /tmp
        #1 sleep for a short random amount of time
        #2 check if a sentinel file /tmp/reference.copying exists.  If yes, sleep for 10 seconds and go to 1
        #3 check if /tmp/reference.stidx and /tmp/reference.sthash exist.  If yes, done & start stampy
        #4 create the sentinel file
        #5 copy hash and index to /tmp
        #6 remove sentinel file, done & start stampy
        my $g = $ref.'.stidx';
        my $h = $ref.'.sthash';
        my $fsu = VertRes::Utils::FileSystem->new();
        my $md5 = $fsu->calculate_md5($ref);
        my $local_ref = '/tmp/VertRes_'.$md5.'_stampy_ref';
        my $local_g = $local_ref.'.stidx';
        my $local_h = $local_ref.'.sthash';
        my $sentinal = '/tmp/VertRes_'.$md5.'_copying_stampy_ref_files';
        my $in_use_base = $local_ref.'.inuse';
        my $in_use = $in_use_base.'.'.$$;
        
        my $max_checks = 400;
        sleep(int(rand(14)) + 1);
        my $checks = 0;
        while (-e $sentinal) {
            $checks++;
            sleep(int(rand(14)) + 1);
            $self->throw("waited for another job to finish coping reference files, but now giving up") if $checks >= $max_checks;
        }
        
        my $start_stampy = 0;
        if (-s $local_g && -s $local_h) {
            $start_stampy = 1;
        }
        else {
            $self->register_for_unlinking($sentinal);
            open(my $sfh, '>', $sentinal) || $self->throw("could not create sentinal file");
            close($sfh);
            copy($g, $local_g) || $self->throw("could not copy $g to $local_g");
            copy($h, $local_h) || $self->throw("could not copy $h to $local_h");
            $start_stampy = 1;
            unlink($sentinal);
        }
        
        if ($start_stampy) {
            $self->register_for_unlinking($in_use);
            open(my $iufh, '>', $in_use) || $self->throw("Could not write to $in_use");
            close($iufh);
            my $bwa_options = '';
            my $aln_q = delete $other_args{aln_q};
            unless (defined $aln_q) {
                $aln_q = 15;
            }
            my $sampe_a = delete $other_args{sampe_a};
            #if ($sampe_a) {
            #    $bwa_options .= "-a $sampe_a ";
            #}
            $bwa_options .= "-q $aln_q $ref";
            $self->simple_run("--bwaoptions=\"$bwa_options\" -g $local_ref -h $local_ref -M $fqs_string > $out");
            unlink($in_use);
            
            # if no other pid is using the ref files, delete them
            open(my $lsfh, 'ls $in_use_base.* 2> /dev/null |');
            my $others = 0;
            while (<$lsfh>) {
                $others++;
            }
            close($lsfh);
            
            unless ($others) {
                unlink($local_g);
                unlink($local_h);
            }
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
