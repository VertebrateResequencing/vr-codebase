=head1 NAME

VertRes::Wrapper::ssaha - wrapper for ssaha2

=head1 SYNOPSIS

use VertRes::Wrapper::ssaha;

my $wrapper = VertRes::Wrapper::ssaha->new();

...

# Or just run the full set of commands needed to do mapping:
$wrapper->do_mapping(ref => 'ref.fa',
                     read0 => 'reads.fastq'
                     read1 => 'reads_1.fastq',
                     read2 => 'reads_2.fastq',
                     output => 'output.sam',
                     index_a => 'bwtsw',
                     sampe_a => 2000);

=head1 DESCRIPTION

Runs ssaha2 in a nice way. Encapsulates the series of commands that must be
run into just one command. (and let's you run individual commands, but not all
have been wrapped yet...)


For subsequent analysis you'd convert the output to a .sam and then to a .bam
and sort it:

use VertRes::Wrapper::samtools;
my $sam_wrapper = VertRes::Wrapper::samtools->new();
$sam_wrapper->faidx('ref.fa');
$sam_wrapper->view('aln.sam', 'aln.bam', b => 1, t => 'ref.fa.fai');
$sam_wrapper->sort('aln.bam' 'aln_sorted');

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Wrapper::ssaha;

use strict;
use warnings;
use VertRes::IO;
use VertRes::Utils::FileSystem;
use File::Basename;
use File::Spec;
use File::Copy;
use Cwd 'abs_path';
use VertRes::Utils::Cigar;
use VertRes::Utils::FastQ;
use SamTools;

use base qw(VertRes::Wrapper::WrapperI);

our @ref_hash_suffixes = ('head', 'body', 'name', 'base', 'size');

=head2 new

 Title   : new
 Usage   : my $wrapper = VertRes::Wrapper::ssaha->new();
 Function: Create a VertRes::Wrapper::ssaha object.
 Returns : VertRes::Wrapper::ssaha object
 Args    : quiet => boolean

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args, exe => 'ssaha2');
    
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
    open(my $fh, "$exe -v |") || $self->throw("Could not start $exe");
    my $version = 0;
    while (<$fh>) {
        if (/ SSAHA2 version (\S+)/) {
            $version = $1;
            last;
        }
    }
    close($fh);
    
    return $version;
}

=head2 ssaha2Build

 Title   : ssaha2Build
 Usage   : $wrapper->ssaha2Build('ref.fa', 'output_basename', skip => 3);
 Function: Index the reference genome.
 Returns : n/a
 Args    : input reference fasta, output basename, hash of options understood
           by ssaha2Build (currently limited to kmer and skip. NB: skip must
           match the skip parameter you later plan to use when mapping)

=cut

sub ssaha2Build {
    my ($self, $in_fa, $out_basename, %opts) = @_;
    $in_fa = abs_path($in_fa);
    my (undef, $path) = fileparse($in_fa);
    $out_basename = File::Spec->catfile($path, basename($out_basename));
    
    $self->exe('ssaha2Build');
    
    $self->switches([]);
    $self->params([qw(kmer skip save)]);
    $self->_set_params_and_switches_from_args(%opts, save => $out_basename);
    
    foreach my $suffix (@ref_hash_suffixes) {
        $self->register_output_file_to_check($out_basename.'.'.$suffix);
    }
    
    $self->{silence_stdout} = 1;
    my @results = $self->run($in_fa);
    $self->{silence_stdout} = 0;
    return @results;
}

=head2 ssaha2

 Title   : ssaha2
 Usage   : $wrapper->ssaha2(\@inputs, $output, skip => 3);
 Function: Align the @inputs.
 Returns : n/a
 Args    : input sequence file name(s) in an array reference, output name (undef
           if using run_method('open')), hash of options understood by ssaha2
           (NB: skip must match the skip parameter you set in ssaha2Build if
           using the save option)

=cut

sub ssaha2 {
    my ($self, $seqs, $out, %opts) = @_;
    
    $self->exe('ssaha2');
    
    $self->switches([qw(sense best 454 NQS tags name fix solexa)]);
    $self->params([qw(kmer skip save ckmer cmatch cut seeds depth memory score
                      identity port align edge array start end quality output
                      diff udiff disk weight rtype pair outfile mthresh)]);
    $self->_set_params_and_switches_from_args(%opts);
    
    # we pipe via tee because some sort of buffer issue prevents any output
    # to file or straight pipe to perl normally
    my @args = (@{$seqs}, ' | tee');
    if ($out) {
        $self->register_output_file_to_check($out);
        push(@args, " > $out");
    }
    
    return $self->run(@args);
}

=head2 do_mapping

 Title   : do_mapping
 Usage   : $wrapper->do_mapping(ref => 'ref.fa',
                                read0 => 'reads.fastq'
                                read1 => 'reads_1.fastq',
                                read2 => 'reads_2.fastq',
                                output => 'output.sam',
                                local_cache => '/local/space/for/files',
                                read_group => 'SRR000000');
 Function: Run ssaha2 on the supplied files, generating a sam file of the
           mapping. Checks the sam file isn't truncated. Optimises for 454.
           Filters short reads. Does custom read pair matching. Indexes the
           reference with ssaha2Build if not already done.
 Returns : n/a
 Args    : required options:
           ref => 'ref.fa'
           output => 'output.sam'
           insert_size => int (default 2000)
           local_cache => '/path' (defaults to standard tmp space)
           read_group => string (defaults to undef, so your sam will have no
                                 RG tag!)

           read1 => 'reads_1.fastq', read2 => 'reads_2.fastq'
           -or-
           read0 => 'reads.fastq'

           Optionally, also supply opts as understood by ssaha2Build() and
           ssaha2(), but prefix the option name with the name of the command,
           eg. ssaha2Build_skip => 3, ssaha2_454 => 1. Many sensible defaults
           are set - best not to supply any of your own, really.

=cut

sub do_mapping {
    my ($self, %args) = @_;
    
    my $io = VertRes::IO->new();
    my $fsu = VertRes::Utils::FileSystem->new();
    my $cigar_util = VertRes::Utils::Cigar->new();
    
    my $orig_run_method = $self->run_method;
    $self->run_method('system');
    
    my $read_group = delete $args{read_group};
    
    my $ref_fa = delete $args{ref} || $self->throw("ref is required");
    my $ref_fa_hash_base = $ref_fa;
    #$ref_fa_hash_base =~ s/\.fa[^.]*$//;
    my $insert_size = delete $args{insert_size} || 2000;
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
        if ($key =~ /^(ssaha2Build|ssaha2)_(\w+)/) {
            $command_opts{$1}->{$2} = $val;
        }
    }
    
    # copy reference hash files to local cache; create hash files first if
    # necessary
    my $local_cache = delete $args{local_cache};
    unless ($local_cache) {
        $local_cache = File::Spec->tmpdir();
        # == /tmp : ideally we don't want these files deleted when we're done
        # with them, so we don't clean them up by using VertRes::Utils::FileSystem->tempdir
    }
    my $hash_files_found = 0;
    foreach my $suffix (@ref_hash_suffixes) {
        $hash_files_found++ if -s "$ref_fa_hash_base.$suffix";
    }
    unless ($hash_files_found == 5) {
        $self->ssaha2Build($ref_fa, $ref_fa_hash_base, skip => 3, %{$command_opts{ssah2Build}});
        $self->throw("failed during the ref hash index build step, giving up for now") unless $self->run_status >= 1;
    }
    my $hash_basename = basename($ref_fa_hash_base);
    foreach my $suffix (@ref_hash_suffixes) {
        my $source = "$ref_fa_hash_base.$suffix";
        my $dest = $fsu->catfile($local_cache, "$hash_basename.$suffix");
        
        # if dest already exists, check it is really our file
        my $do_copy = 1;
        if (-s $dest) {
            my $diff = `diff $source $dest`;
            unless ($diff) {
                $do_copy = 0;
            }
        }
        
        # copy if it doesn't exist or isn't correct
        if ($do_copy) {
            my $ok = $fsu->copy($source, $dest);
            $ok || $self->throw("failed prior to attempting the mapping (could not copy $source -> $dest)");
        }
    }
    $ref_fa_hash_base = $fsu->catfile($local_cache, $hash_basename);
    
    unless (-s $out_sam) {
        my $tmp_sam = $out_sam.'_tmp';
        
        my (undef, $out_dir) = fileparse($tmp_sam);
        my (@filtered_fastqs, @cigar_outputs);
        $self->run_method('open');
        
        foreach my $fastq (@fqs) {
            my $temp_dir = $fsu->tempdir();
            
            my $fq_basename = basename($fastq);
            $fq_basename =~ s/\.gz$//;
            my $tmp_fastq = $fsu->catfile($temp_dir, $fq_basename);
            
            my $cigar_name = $fq_basename;
            $cigar_name =~ s/\.fastq$//;
            my $cigar_out = $fsu->catfile($out_dir, $cigar_name.'.cigar.gz');
            
            unless (-s $cigar_out) {
                # filter out short reads from fastq
                my $fu = VertRes::Utils::FastQ->new();
                $fu->filter_reads($fastq, $tmp_fastq, min_length => 30);
                
                # run ssaha2, filtering the output to get the top 10 hits per read,
                # grouping by readname, and compressing it
                my $sfh = $self->ssaha2([$tmp_fastq], undef, disk => 1, '454' => 1, output => 'cigar', diff => 10, save => $ref_fa_hash_base);
                
                my $tmp_cigar = $cigar_out;
                $tmp_cigar =~ s/cigar\.gz$/cigar.tmp.gz/;
                open(my $cfh, "| gzip -c > $tmp_cigar");
                $cigar_util->top_10_hits_per_read($sfh, $cfh);
                
                # check for completion
                my $done = `zcat $tmp_cigar | tail -1 | grep "SSAHA2 finished" | wc -l`;
                ($done) = $done =~ /^(\d+)/;
                if ($done) {
                    move($tmp_cigar, $cigar_out) || $self->throw("Failed to move $tmp_cigar to $cigar_out: $!");
                }
                else {
                    $self->throw("ssah2 failed to finish for fastq $fastq");
                }
            }
            
            push(@filtered_fastqs, $tmp_fastq);
            push(@cigar_outputs, $cigar_out);
        }
        
        # convert to sam
        my $expected_sam_lines = 0;
        if (@cigar_outputs == 1 || @cigar_outputs == 2) {
            my $print_group = $read_group ? $read_group : 'undef';
            $self->debug("will do cigar_to_sam([@fqs], [@cigar_outputs], $insert_size, $print_group, $tmp_sam)");
            $expected_sam_lines = $cigar_util->cigar_to_sam(\@fqs, \@cigar_outputs, $insert_size, $read_group, $tmp_sam);
        }
        else {
            $self->throw("Something went wrong, ended up with ([@filtered_fastqs], [@cigar_outputs])");
        }
        
        # sanity check the sam
        $io->file($tmp_sam);
        my $actual_sam_lines = $io->num_lines;
        if ($actual_sam_lines == $expected_sam_lines) {
            move($tmp_sam, $out_sam) || $self->throw("Failed to move $tmp_sam to $out_sam: $!");
            $self->_set_run_status(1);
        }
        else {
            $self->warn("a sam file was made by do_mapping(), but it only had $actual_sam_lines lines, not the expected $expected_sam_lines - will unlink it");
            $self->_set_run_status(-1);
            unlink("$tmp_sam");
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
 Function: Run your chosen ssaha command on the supplied file(s).
 Returns : n/a
 Args    : paths to input/output files

=cut

1;
