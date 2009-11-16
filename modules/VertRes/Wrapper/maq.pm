=head1 NAME

VertRes::Wrapper::maq - wrapper for maq

=head1 SYNOPSIS

use VertRes::Wrapper::maq;

my $wrapper = VertRes::Wrapper::maq->new();

# Now you can call a method corrseponding to a maq command.
# All methods take the filenames as a list to start with, followed by a hash
# of any options the maq command understands.

$wrapper->fasta2bfa('ref.fa', 'ref.fa.bfa');
$wrapper->fastq2bfq('reads.fastq', 'reads.fastq.bfq');
$wrapper->map('output.map', 'ref.fa.bfa', ['reads.fastq.bfq'], a => 200);

# Or just run the full set of commands needed to do mapping:
$wrapper->do_mapping(ref => 'ref.fa',
                     read1 => 'reads_1.fastq',
                     read2 => 'reads_2.fastq',
                     output => 'output.sam',
                     a => 200);

=head1 DESCRIPTION

Runs maq in a nice way. Encapsulates the series of commands that must be
run into just one command. (and let's you run individual commands, but not all
have been wrapped yet...)

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Wrapper::maq;

use strict;
use warnings;
use File::Copy;
use File::Basename;
use VertRes::IO;
use VertRes::Wrapper::samtools;
use VertRes::Wrapper::fastqcheck;
use VertRes::Parser::fastqcheck;
use VertRes::Utils::FastQ;
use VertRes::Utils::Sam;

use base qw(VertRes::Wrapper::WrapperI);

#read length ranges and what to set the -e parameter in maq to
my %e_parameter = (
    37 	=> 80, 	#2 * 40
    63	=> 120,	#3 * 40
    92 	=> 160, #4 * 40
    123	=> 200, #5 * 40
    156	=> 240 	#6 * 40
);

=head2 new

 Title   : new
 Usage   : my $wrapper = VertRes::Wrapper::maq->new();
 Function: Create a VertRes::Wrapper::maq object.
 Returns : VertRes::Wrapper::maq object
 Args    : quiet   => boolean

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(exe => 'maq', @args);
    
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

=head2 fasta2bfa

 Title   : fasta2bfa
 Usage   : $wrapper->fasta2bfa('ref.fa', 'ref.fa.bfa');
 Function: Create an index of your reference fasta.
 Returns : n/a
 Args    : input fasta, output bfq

=cut

sub fasta2bfa {
    my ($self, $in_fa, $out_bfa) = @_;
    
    $self->exe($self->{orig_exe}.' fasta2bfa');
    
    $self->switches([]);
    $self->params([]);
    $self->_set_params_and_switches_from_args();
    
    $self->register_output_file_to_check($out_bfa);
    
    return $self->run($in_fa, $out_bfa);
}


=head2 fastq2bfq

 Title   : fastq2bfq
 Usage   : $wrapper->fastq2bfq('fastq.fa', 'fastq.fa.bfq');
 Function: Create an index of your fastq sequence reads.
 Returns : n/a
 Args    : input fasta, output bfq, optionally n => \d+ to set the nreads
           option

=cut

sub fastq2bfq {
    my ($self, $in_fa, $out_bfq, %args) = @_;
    
    unless ($out_bfq =~ /\.bfq$/) {
        $out_bfq .= '.bfq';
    }
    
    $self->exe($self->{orig_exe}.' fastq2bfq');
    
    $self->switches([]);
    $self->params(['n']);
    $self->_set_params_and_switches_from_args(%args);
    
    $self->register_output_file_to_check($out_bfq);
    
    return $self->run($in_fa, $out_bfq);
}

=head2 mapstat

 Title   : mapstat
 Usage   : $wrapper->mapstat('mapping.map', 'mapping.map.mapstat');
 Function: Generate some statistics describing a mapping carried out with map().
 Returns : n/a
 Args    : input map file, output stat file

=cut

sub mapstat {
    my ($self, $in, $out) = @_;
    
    $self->exe($self->{orig_exe}.' mapstat');
    
    $self->switches([]);
    $self->params([]);
    $self->_set_params_and_switches_from_args();
    
    $self->register_output_file_to_check($out);
    
    return $self->run($in, ' > '.$out);
}

=head2 map

 Title   : map
 Usage   : $wrapper->map('output.map', 'ref.bfa', ['reads1_bfq', 'reads2_bfq'],
                         %opts);
 Function: map...
 Returns : n/a
 Args    : list of file paths (the first being the output file, requiring
           that you've already indexed it with index()),
           options can be supplied as a hash, with keys corresponding to maq
           map option names, and values the value you wish to use, eg.
           1 => 37, 2 => 37, a => 200

=cut

sub map {
    my ($self, $out_map, $ref, $bfqs, %opts) = @_;
    
    $self->exe($self->{orig_exe}.' map');
    
    $self->switches([qw(W t c)]);
    $self->params([qw(1 2 m e d a A n M u H C s)]);
    
    $self->_set_params_and_switches_from_args(%opts);
    
    $self->register_output_file_to_check($out_map);
    
    return $self->run($out_map, $ref, @{$bfqs});
}

=head2 do_mapping

 Title   : do_mapping
 Usage   : $wrapper->do_mapping(ref => 'ref.fa',
                                read1 => 'reads_1.fastq',
                                read2 => 'reads_2.fastq',
                                output => 'output.sam',
                                a => 1000);
 Function: Run maq on the supplied files, generating a sam file of the mapping.
           Checks the sam file isn't truncated. The sam file will contain
           all reads, including those there were not mapped.
           maq map options 1,2, and e are determined automatically from your
           reads and set unless supplied by you.
 Returns : n/a
 Args    : required options:
           ref => 'ref.fa'
           output => 'output.sam'

           read1 => 'reads_1.fastq', read2 => 'reads_2.fastq'
           -or-
           read0 => 'reads.fastq'

           Optionally, also supply opts as understood by map(). Note that if you
           know an expected insert size, provide this as the 'a' option.
           Internally a and A will be set appropriately based on that. '1' and
           '2' are also set automatically unless supplied, by calculating a
           suitable clip point with VertRes::Utils::FastQ->clip_point.

=cut

sub do_mapping {
    my ($self, %args) = @_;
    
    my $orig_run_method = $self->run_method;
    $self->run_method('system');
    
    my $ref_fa = delete $args{ref} || $self->throw("ref is required");
    my @fqs;
    my $paired = 0;
    if (defined $args{read1} && defined $args{read2} && ! defined $args{read0}) {
        push(@fqs, delete $args{read1});
        push(@fqs, delete $args{read2});
        $paired = 1;
    }
    elsif (defined $args{read0} && ! defined $args{read1} && ! defined $args{read2}) {
        push(@fqs, delete $args{read0});
    }
    else {
        $self->throw("Bad read args: need read1 and read2, or read0");
    }
    my $out_sam = delete $args{output};
    
    # fasta2bfa for reference
    my $bfa = $ref_fa.'.bfa';
    unless (-s $bfa) {
        $self->fasta2bfa($ref_fa, $bfa);
        $self->throw("failed during the fasta2bfa step, giving up for now") unless $self->run_status >= 1;
    }
    
    # fastq2bfq for the fastqs
    my ($bfq1, $bfq2) = @fqs;
    $bfq1 =~ s/\.f[^.]+(?:\.gz)?$/.bfq/;
    $bfq2 =~ s/\.f[^.]+(?:\.gz)?$/.bfq/ if $bfq2;
    foreach my $fs ([$fqs[0], $bfq1], $fqs[1] ? [$fqs[1], $bfq2] : []) {
        my ($fq, $bfq) = @{$fs};
        $fq || next;
        
        unless (-s $bfq) {
            $self->fastq2bfq($fq, $bfq);
            $self->throw("failed during the fasta2bfq step for $fq, giving up for now") unless $self->run_status >= 1;
        }
    }
    
    # map
    my $out_map = $out_sam;
    $out_map =~ s/\.sam$/.map/;
    
    # our output sam must contain all reads, including the unmapped, so we
    # force 'u' on
    unless (defined $args{u}) {
        $args{u} = $out_map.'.unmapped';
    }
    
    # sometimes we won't use all input fastqs, so we need to make a note of that
    my $skipped_fqs_file = $out_map.'.skipped_fastqs';
    
    unless (-s $out_map) {
        my @bfqs;
        if (@fqs == 2) {
            @bfqs = ($bfq1, $bfq2);
        }
        else {
            @bfqs = ($bfq1);
        }
        
        # work out 'e', '1' and '2'. If these turn out to be less than maq's
        # minimum size of 12, don't try to map
        my @final_bfqs = @bfqs;
        unless (defined $args{e} && defined $args{1} && (@fqs > 1 ? defined $args{2} : 1)) {
            my @fqcps = $self->_get_fastq_details(@fqs);
            my @read_lengths;
            my $fqu = VertRes::Utils::FastQ->new();
            my $i = 0;
            foreach my $fqcp (@fqcps) {
                # calculate the clip point
                my $length = $fqu->clip_point($fqs[$i]);
                if ($length == -1) {
                    $length = int($fqcp->avg_length);
                }
                else {
                    $self->warn("Set clip point for $fqs[$i] to $length");
                }
                
                push(@read_lengths, $length);
                $i++;
            }
            foreach (sort { $a <=> $b } keys %e_parameter) {
                if ($read_lengths[0] < $_ && ($read_lengths[1] ? $read_lengths[1] < $_ : 1)) {
                    $args{e} = $e_parameter{$_} unless defined $args{e};
                    
                    @final_bfqs = ();
                    if ($read_lengths[0] <= 12) {
                        $self->warn("fastq file $fqs[0] will be skipped, since after hard-clipping to position $read_lengths[0] it is too small to map");
                        open(my $skipped_fh, '>>', $skipped_fqs_file) || $self->throw("Couldn't write to $skipped_fqs_file");
                        print $skipped_fh 0, "\n";
                        close($skipped_fh);
                    }
                    else {
                        $args{1} = $read_lengths[0];
                        push(@final_bfqs, $bfqs[0]);
                    }
                    
                    if ($read_lengths[1]) {
                        if ($read_lengths[1] <= 12) {
                            $self->warn("fastq file $fqs[1] will be skipped since, after hard-clipping to position $read_lengths[1], sequences are too short to map");
                            open(my $skipped_fh, '>>', $skipped_fqs_file) || $self->throw("Couldn't write to $skipped_fqs_file");
                            print $skipped_fh 1, "\n";
                            close($skipped_fh);
                        }
                        else {
                            $args{2} = $read_lengths[1];
                            push(@final_bfqs, $bfqs[1]);
                        }
                    }
                    
                    last;
                }
            }
        }
        
        # work out 'a' and 'A'
        if (@final_bfqs > 1) {
            my $a = $args{a};
            # the aim with insert size settings is to capture everything within
            # a reasonable distance as being paired. So in the first instance
            # we just set a to 1000
            $args{a} = 1000;
            
            # now, if we expect an insert size of greater than 1000, we probably
            # have big inserts and reads pointing outward, so we set A
            if ($a && $a > 1000) {
                $args{A} = $a * 2;
            }
        }
        else {
            delete $args{a};
            delete $args{A};
        }
        
        if (@final_bfqs) {
            $self->map($out_map, $bfa, \@final_bfqs, %args);
            $self->throw("failed during the map, giving up for now") unless $self->run_status >= 1;
        }
        else {
            $self->warn("none of the input fastqs had long enough sequences after clipping to map; will generate a pure unmapped sam file");
        }
    }
    
    my $io = VertRes::IO->new();
    
    unless (-s $out_sam) {
        my $tmp_sam = $out_sam.'_tmp';
        
        # convert to sam *** should create a maq2sam wrapper...
        system("maq2sam-long $out_map > $tmp_sam");
        
        # append unmapped reads
        my $unmapped = $args{u};
        if (-s $unmapped) {
            open(my $sfh, '>>', $tmp_sam) || $self->throw("Could not append to $tmp_sam");
            open(my $ufh, $unmapped) || $self->throw("Could not open $unmapped");
            my $su = VertRes::Utils::Sam->new();
            my %seen_reads;
            while (<$ufh>) {
                my ($id, undef, $seq, $qual) = split;
                
                # work out the SAM flags
                my $flags;
                if (exists $seen_reads{$id}) {
                    $flags = $su->calculate_flag(paired_tech => $paired,
                                                 $paired ? (mate_unmapped => 1) : (),
                                                 self_unmapped => 1,
                                                 $paired ? ('2nd_in_pair' => 1) : ());
                }
                else {
                    $flags = $su->calculate_flag(paired_tech => $paired,
                                                 $paired ? (mate_unmapped => 1) : (),
                                                 self_unmapped => 1,
                                                 $paired ? ('1st_in_pair' => 1) : ());
                    $seen_reads{$id} = 1;
                }
                
                print $sfh "$id\t$flags\t*\t0\t0\t*\t*\t0\t0\t$seq\t$qual\n";
            }
        }
        if (-s $skipped_fqs_file) {
            open(my $skipped_fh, '<', $skipped_fqs_file) || $self->throw("Couldn't open $skipped_fqs_file");
            my @skipped_fqs;
            while (<$skipped_fh>) {
                chomp;
                next unless /^\d+$/;
                push(@skipped_fqs, $fqs[$_]);
            }
            
            open(my $sfh, '>>', $tmp_sam) || $self->throw("Could not append to $tmp_sam");
            my $su = VertRes::Utils::Sam->new();
            my %seen_reads;
            foreach my $fq (@skipped_fqs) {
                my $fp = VertRes::Parser::fastq->new(file => $fq);
                my $rh = $fp->result_holder;
                
                while ($fp->next_result) {
                    my ($id, $seq, $qual) = @{$rh};
                    
                    # work out the SAM flags
                    my $flags;
                    if (exists $seen_reads{$id}) {
                        $flags = $su->calculate_flag(paired_tech => $paired,
                                                     mate_unmapped => 1,
                                                     self_unmapped => 1,
                                                     '2nd_in_pair' => 1);
                    }
                    else {
                        $flags = $su->calculate_flag(paired_tech => $paired,
                                                     mate_unmapped => 1,
                                                     self_unmapped => 1,
                                                     '1st_in_pair' => 1);
                        $seen_reads{$id} = 1;
                    }
                    
                    print $sfh "$id\t$flags\t*\t0\t0\t*\t*\t0\t0\t$seq\t$qual\n";
                }
            }
        }
        
        # check the sam file isn't truncated
        my $num_reads = 0;
        my @fqcps = $self->_get_fastq_details(@fqs);
        foreach my $fqcp (@fqcps) {
            $num_reads += $fqcp->num_sequences;
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
        }
    }
    else {
        $self->_set_run_status(1);
    }
    
    $self->run_method($orig_run_method);
    return;
}

sub _get_fastq_details {
    my ($self, @fqs) = @_;
    
    my $io = VertRes::IO->new();
    
    # make temp fastqcheck files if they don't already exist
    my @fqcs;
    foreach my $fq (@fqs) {
        my $fqc = $fq.'.fastqcheck';
        my $older = $self->_is_older($fqc, $fq);
        
        unless (-s $fqc && $older) {
            my $fqcw = VertRes::Wrapper::fastqcheck->new(quiet => 1);
            my $tempdir = $io->tempdir;
            $fqc = $io->catfile($tempdir, basename($fq).'.fastqcheck');
            if ($self->_is_older($fq, $fqc)) {
                unlink($fqc);
            }
            $fqcw->run($fq, $fqc) unless -s $fqc;
        }
        
        push(@fqcs, $fqc);
    }
    
    # parse the fastqcheck files, returning the parser objects
    my @fqcps;
    foreach my $fqc (@fqcs) {
        -s $fqc || $self->throw("Expected a fastqcheck file '$fqc', but it wasn't there!");
        my $fqcp = VertRes::Parser::fastqcheck->new(file => $fqc);
        push(@fqcps, $fqcp);
    }
    
    @fqs == @fqcps || $self->throw("Didn't get as many fastqcheck parsers as input fastqs!");
    
    return @fqcps;
}

=head2 run

 Title   : run
 Usage   : Do not call directly: use one of the other methods like index()
           instead.
 Function: Run your chosen maq command on the supplied file(s).
 Returns : n/a
 Args    : paths to input/output files

=cut

1;
