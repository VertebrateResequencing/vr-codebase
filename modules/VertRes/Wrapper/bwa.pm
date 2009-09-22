=head1 NAME

VertRes::Wrapper::bwa - wrapper for bwa

=head1 SYNOPSIS

use VertRes::Wrapper::bwa;

my $wrapper = VertRes::Wrapper::bwa->new();

# Now you can call a method corrseponding to a bwa command.
# All methods take the filenames as a list to start with, followed by a hash
# of any options the bwa command understands.

$wrapper->index('ref.fa', a => 'bwtsw');
$wrapper->aln('ref.fa', 'reads.fastq', 'output.sai');
$wrapper->sampe('ref.fa', 'aln1.sai', 'aln2.sai',
                'reads1.fastq', 'reads2.fastq', 'output.sam', a => 2000);
# etc.

# Or just run the full set of commands needed to do mapping:
$wrapper->do_mapping(ref => 'ref.fa',
                     read1 => 'reads_1.fastq',
                     read2 => 'reads_2.fastq',
                     output => 'output.sam',
                     index_a => 'bwtsw',
                     sampe_a => 2000);

=head1 DESCRIPTION

Runs bwa in a nice way. Encapsulates the series of commands that must be
run into just one command. (and let's you run individual commands, but not all
have been wrapped yet...)

First creates BWT index of fasta reference sequence, using bwtsw algorithm since
only it can handle human genome:
bwa index -a bwtsw ref.fa
Then generates alignments in suffix array (SA) coordinates:
bwa aln ref.fa reads_1.fastq > aln_1.sai
bwa aln ref.fa reads_2.fastq > aln_2.sai
Finally converts SA coordinates to chromosomal coords and does pairing:
bwa sampe ref.fa aln_1.sai aln_2.sai reads_1.fastq reads_2.fastq > aln.sam

For subsequent analysis you'd convert the .sam to a .bam and sort it:
use VertRes::Wrapper::samtools;
my $sam_wrapper = VertRes::Wrapper::samtools->new();
$sam_wrapper->faidx('ref.fa');
$sam_wrapper->view('aln.sam', 'aln.bam', b => 1, t => 'ref.fa.fai');
$sam_wrapper->sort('aln.bam' 'aln_sorted');

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Wrapper::bwa;

use strict;
use warnings;
use File::Copy;
use VertRes::IO;

use base qw(VertRes::Wrapper::WrapperI);


=head2 new

 Title   : new
 Usage   : my $wrapper = VertRes::Wrapper::bwa->new();
 Function: Create a VertRes::Wrapper::bwa object.
 Returns : VertRes::Wrapper::bwa object
 Args    : quiet   => boolean

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args, exe => 'bwa');
    
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
    open(my $fh, "$exe 2>&1 |") || $self->throw("Could not start $exe");
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

=head2 index

 Title   : index
 Usage   : $wrapper->index('in.fa', %opts);
 Function: index...
 Returns : n/a (generates a bunch of files with different suffixes to your
                input fasta)
 Args    : list of file paths, options as a hash. does a => 'bwtsw' by default.
           (the alternative is a => 'is' for small genomes)

=cut

sub index {
    my ($self, $in_fa, %opts) = @_;
    
    $self->exe('bwa index');
    
    $self->switches([qw(c)]);
    $self->params([qw(a p)]);
    $self->_set_params_and_switches_from_args(a => 'bwtsw', %opts);
    
    my $out_file = $in_fa;
    $out_file .= '.bwt';
    $self->register_output_file_to_check($out_file);
    
    return $self->run($in_fa);
}

=head2 aln

 Title   : aln
 Usage   : $wrapper->aln('ref.fa', 'reads.fastq', 'output.sai', %opts);
 Function: aln...
 Returns : n/a
 Args    : list of file paths (the first being the reference fasta, requiring
           that you've already indexed it with index()), options as a hash.

=cut

sub aln {
    my ($self, $ref, $fastq, $out_sai, %opts) = @_;
    
    $self->exe('bwa aln');
    
    $self->switches([qw(c L N)]);
    $self->params([qw(n o e i d l k m t M O E R)]);
    $self->_set_params_and_switches_from_args(%opts);
    
    $self->register_output_file_to_check($out_sai);
    
    return $self->run($ref, $fastq, ' > '.$out_sai);
}

=head2 sampe

 Title   : sampe
 Usage   : $wrapper->sampe('ref.fa', 'aln1.sai', 'aln2.sai', 'reads1.fastq',
                           'reads2.fastq, 'output.sam', %opts);
 Function: sampe...
 Returns : n/a
 Args    : list of file paths (the first being the reference fasta, requiring
           that you've already indexed it with index()), options as a hash. NB:
           it is probably important to change the 'a' option, as by default this
           is a max insert size of just 500.

=cut

sub sampe {
    my ($self, $ref, $sai1, $sai2, $fq1, $fq2, $out_sam, %opts) = @_;
    
    $self->exe('bwa sampe');
    
    $self->switches([qw(s)]);
    $self->params([qw(a o)]);
    
    $self->_set_params_and_switches_from_args(%opts);
    
    $self->register_output_file_to_check($out_sam);
    
    return $self->run($ref, $sai1, $sai2, $fq1, $fq2, ' > '.$out_sam);
}

=head2 samse

 Title   : samse
 Usage   : $wrapper->sampe('ref.fa', 'aln1.sai', 'reads1.fastq',
                           'output.sam', %opts);
 Function: samse...
 Returns : n/a
 Args    : list of file paths (the first being the reference fasta, requiring
           that you've already indexed it with index()), options as a hash.

=cut

sub samse {
    my ($self, $ref, $sai1, $fq1, $out_sam, %opts) = @_;
    
    $self->exe('bwa samse');
    
    $self->switches([]);
    $self->params([qw(n)]);
    
    $self->_set_params_and_switches_from_args(%opts);
    
    $self->register_output_file_to_check($out_sam);
    
    return $self->run($ref, $sai1, $fq1, ' > '.$out_sam);
}

=head2 do_mapping

 Title   : do_mapping
 Usage   : $wrapper->do_mapping(ref => 'ref.fa',
                                read1 => 'reads_1.fastq',
                                read2 => 'reads_2.fastq',
                                output => 'output.sam',
                                index_a => 'bwtsw',
                                sampe_a => 2000);
 Function: Run bwa on the supplied files, generating a sam file of the mapping.
           Checks the sam file isn't truncated.
 Returns : n/a
 Args    : required options:
           ref => 'ref.fa'
           output => 'output.sam'

           read1 => 'reads_1.fastq', read2 => 'reads_2.fastq'
           -or-
           read0 => 'reads.fastq'

           Optionally, also supply opts as understood by index(), aln(), samse()
           or sampe(), but prefix the option name with the name of the command,
           eg. index_a => 'is', sampe_a => 500.

=cut

sub do_mapping {
    my ($self, %args) = @_;
    
    my $orig_run_method = $self->run_method;
    $self->run_method('system');
    
    my $ref_fa = delete $args{ref} || $self->throw("ref is required");
    my @fqs;
    if (defined $args{read1} && defined $args{read2} && ! defined $args{read0}) {
        push(@fqs, delete $args{read1});
        push(@fqs, delete $args{read2});
    }
    elsif (defined $args{read0} && ! defined $args{read1} && ! defined $args{read2}) {
        push(@fqs, delete $args{read0});
    }
    else {
        $self->throw("Bad read args: need read1 and read2, or read0");
    }
    my $out_sam = delete $args{output};
    
    my %command_opts;
    while (my ($key, $val) = each %args) {
        if ($key =~ /^(index|aln|sampe|samse)_(\w+)/) {
            $command_opts{$1}->{$2} = $val;
        }
    }
    
    # index
    my $bwt = $ref_fa.'.bwt';
    unless (-s $bwt) {
        $self->index($ref_fa, %{$command_opts{index} || {}});
        $self->throw("failed during the index step, giving up for now") unless $self->run_status >= 1;
    }
    
    # aln
    my ($sai1, $sai2) = @fqs;
    $sai1 =~ s/\.f[^.]+(?:\.gz)?$/.sai/;
    $sai2 =~ s/\.f[^.]+(?:\.gz)?$/.sai/ if $sai2;
    foreach my $fs ([$fqs[0], $sai1], $fqs[1] ? [$fqs[1], $sai2] : []) {
        my ($fq, $sai) = @{$fs};
        $fq || next;
        
        unless (-s $sai) {
            $self->aln($ref_fa, $fq, $sai, %{$command_opts{aln} || {}});
            $self->throw("failed during the aln step, giving up for now") unless $self->run_status >= 1;
        }
    }
    
    # sampe or sampse
    unless (-s $out_sam) {
        my $tmp_sam = $out_sam.'_tmp';
        
        if (@fqs == 2) {
            $self->sampe($ref_fa, $sai1, $sai2, @fqs, $tmp_sam, %{$command_opts{sampe} || {}});
        }
        else {
            $self->samse($ref_fa, $sai1, @fqs, $tmp_sam, %{$command_opts{samse} || {}});
        }
        $self->throw("failed during the sam[sp]e step, giving up for now") unless $self->run_status >= 1;
        
        # check the sam file isn't truncted
        my $num_reads = 0;
        my $io = VertRes::IO->new();
        foreach my $fastq_file (@fqs) {
            $io->file($fastq_file);
            my $fastq_lines = $io->num_lines();
            $num_reads += $fastq_lines / 4;
        }
        $io->file($tmp_sam);
        my $sam_lines = $io->num_lines();
        
        if ($sam_lines >= $num_reads) {
            move($tmp_sam, $out_sam) || $self->throw("Failed to move $tmp_sam to $out_sam: $!");
        }
        else {
            system("cp $tmp_sam tmp.sam");
            $self->warn("a sam file was made by do_mapping(), but it only had $sam_lines lines, not the expected $num_reads - will unlink it");
            $self->_set_run_status(-1);
            unlink("$tmp_sam");
            exit;
        }
    }
    else {
        $self->_set_run_status(1);
    }
    
    $self->run_method($orig_run_method);
    return;
}

=head2 run

 Title   : run
 Usage   : Do not call directly: use one of the other methods like index()
           instead.
 Function: Run your chosen bwa command on the supplied file(s).
 Returns : n/a
 Args    : paths to input/output files

=cut

1;
