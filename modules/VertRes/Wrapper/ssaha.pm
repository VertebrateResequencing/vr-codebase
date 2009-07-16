=head1 NAME

VertRes::Wrapper::ssaha - wrapper for ssaha

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

Runs ssaha in a nice way. Encapsulates the series of commands that must be
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
use File::Basename;
use AssemblyTools;
use SamTools;

use base qw(VertRes::Wrapper::WrapperI);


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

=head2 do_mapping

 Title   : do_mapping
 Usage   : $wrapper->do_mapping(ref => 'ref.fa',
                                read0 => 'reads.fastq'
                                read1 => 'reads_1.fastq',
                                read2 => 'reads_2.fastq',
                                output => 'output.sam');
 Function: Run bwa on the supplied files, generating a sam file of the mapping.
           Checks the sam file isn't truncated.
 Returns : n/a
 Args    : required options:
           ref => 'ref.fa'
           output => 'output.sam'
           insert_size => int (default 2000)

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
    $ref_fa =~ s/\.fa$//;
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
        if ($key =~ /^(index|aln|sampe|samse)_(\w+)/) {
            $command_opts{$1}->{$2} = $val;
        }
    }
    
    my $io = VertRes::IO->new;
    my @files;
    foreach my $fastq (@fqs) {
        my $temp_dir = $io->tempdir();
        
        # filter out short reads from fastq
        my $fq_basename = basename($fastq);
        $fq_basename =~ s/\.gz$//;
        my $tmp_fastq = $io->catfile($temp_dir, $fq_basename);
        AssemblyTools::filterOutShortReads($fastq, 30, $tmp_fastq);
        
        # run ssaha2, filtering the output to get the top 10 hits per read,
        # and compressing it
        my $cigar_out = $io->catfile($temp_dir, 'ssaha.cigar.gz');
        system("ssaha2 -disk 1 -454 -output cigar -diff 10 -save $ref_fa $tmp_fastq | /nfs/users/nfs_s/sb10/src/vert_reseq/user/tk2/miscScripts/filterCigarStreamTop10.pl | gzip -c > $cigar_out");
        
        # check for completion
        my $done = `zcat $cigar_out | tail -50 | grep "SSAHA2 finished" | wc -l`;
        $done || $self->throw("ssah2 failed to finish for fastq $fastq");
        
        #$cmd .= qq{; perl -w -e "use Mapping_454_ssaha;Mapping_454_ssaha::cigarStat( \\"$currentDir/$cigarName\\", \\"$currentDir/$cigarName.mapstat\\");"'};
        
        push(@files, ($tmp_fastq, $cigar_out));
    }
    
    # convert to sam
    if (@files == 2) {
        $self->debug("will do SamTools::ssaha2samUnpaired(@files, undef, $out_sam.'.gz')");
        SamTools::ssaha2samUnpaired(@files, undef, $out_sam.'.gz');
    }
    elsif (@files == 4) {
        $self->debug("will do SamTools::ssaha2samPaired(@files, $insert_size, undef, $out_sam.'.gz')");
        SamTools::ssaha2samPaired(@files, $insert_size, undef, $out_sam.'.gz');
    }
    else {
        $self->throw("Something went wrong, ended up with (@files)");
    }
    
    system("gunzip -f $out_sam.gz");
    
    $self->_set_run_status(1);
    
    $self->run_method($orig_run_method);
    return;
}

=head2 filter_cigar_stream_top10

 Title   : filter_cigar_stream_top10
 Usage   : $wrapper->filter_cigar_stream_top10($infh, $outfh);
 Function: Filters the cigar output of ssaha2, outputting only the top 10 hits
           per read.
 Returns : n/a
 Args    : input and output filehandles

=cut

sub filter_cigar_stream_top10 {
    my ($self, $infh, $outfh) = @_;
    
    my $currentRead = '';
    my @currentHits;
    while (<$infh>) {
        chomp;
        
        if (/^cigar::/) {
            my @s = split(/\s+/, $_);
            
            if (length($currentRead) == 0) {
                $currentRead = $s[1];
                push(@currentHits, $_);
            }
            elsif ($currentRead ne $s[1]) {
                foreach my $i (0..$#currentHits) {
                    print $outfh $currentHits[$i], "\n";
                }
                
                $currentRead = $s[1];
                undef(@currentHits);
                push(@currentHits, $_);
            }
            else {
                push(@currentHits, $_);
            }
        }
        elsif (/^SSAHA2 finished/) {
            print $outfh $_, "\n";
        }
    }
    
    foreach my $i (0..$#currentHits) {
        print $outfh $currentHits[$i], "\n";
    }
    
    close($infh);
    close($outfh);
    
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
