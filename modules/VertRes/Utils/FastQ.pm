=head1 NAME

VertRes::Utils::FastQ - fastq utility functions

=head1 SYNOPSIS

use VertRes::Utils::FastQ;

my $fastq_util = VertRes::Utils::FastQ->new();

$fastq_util->split(\@fastq_files, split_dir => 'split');

=head1 DESCRIPTION

General utility functions for working on or with fastq files.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Utils::FastQ;

use strict;
use warnings;
use POSIX;
use VertRes::IO;
use VertRes::Utils::FileSystem;
use File::Basename;
use VertRes::Parser::fastqcheck;
use VertRes::Parser::fastq;
use Cwd 'abs_path';
use Inline C => Config => FILTERS => 'Strip_POD';

use base qw(VertRes::Base);

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Utils::FastQ->new();
 Function: Create a new VertRes::Utils::FastQ object.
 Returns : VertRes::Utils::FastQ object
 Args    : n/a

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    return $self;
}

=head2 split

 Title   : split
 Usage   : $obj->split(\@fastqs,
                       split_dir => '/path/to/desired/split_dir',
                       chunk_size => 1000000);
 Function: Split the fastq(s) into multiple smaller files. If more than one
           fastq file is supplied, they are treated as a unit with regard to
           the chunk size. So if you had two fastq files where all sequences
           were 100bp long, and you supplied only one of them with a chunk_size
           of 1000, you'd end up with chunks containing 10 sequences each. But
           with both fastq files supplied, you'd end up with chunks containing
           5 sequences each. The idea being, you'll then use the Nth chunk of
           both fastqs at the same time, bringing the total bases to the
           chunk size.
 Returns : int (the number of splits created)
 Args    : array ref of fastq files,
           split_dir => '/path/to/desired/split_dir' (location to store the
                                                      resulting files)
           chunk_size => int (max number of bases per chunk, default 10000000)

=cut

sub split {
    my ($self, $fastqs, %args) = @_;
    my $chunk_size = $args{chunk_size} || 10000000;
    my $split_dir = $args{split_dir} || $self->throw("split_dir must be supplied");
    
    mkdir($split_dir);
    
    my @ins;
    my @outs;
    my $split_num = 1;
    my $io = VertRes::IO->new();
    my $fsu = VertRes::Utils::FileSystem->new();
    my $fqc = VertRes::Parser::fastqcheck->new();
    my $total_bases = 0;
    my $num_fqcs = 0;
    foreach my $fastq_file (@{$fastqs}) {
        my $basename = basename($fastq_file);
        my $prefix = $basename;
        $prefix =~ s/\.f[^.]+(?:\.gz)?$//;
        
        my $in = VertRes::IO->new(file => $fastq_file, verbose => $self->verbose);
        push(@ins, [$fastq_file, $in]);
        
        # if there are corresponding fastqcheck files we'll check them for the
        # total bases in each fastq, and if we'd only generate one chunk, simply
        # symlink the input to the ouput
        my $fastqcheck = $fastq_file.'.fastqcheck';
        if (-s $fastqcheck) {
            $num_fqcs++;
            $fqc->file($fastqcheck);
            $total_bases += $fqc->total_length();
        }
        
        my $split_file = $fsu->catfile($split_dir, "$prefix.$split_num.fastq.gz");
        my $out = VertRes::IO->new(file => ">$split_file", verbose => $self->verbose);
        push(@outs, [$prefix, $out, $split_file]);
    }
    
    if ($num_fqcs == @ins) {
        my $splits = ceil($total_bases / $chunk_size);
        if ($splits == 1) {
            foreach my $i (0..$#ins) {
                my $fastq_file = $ins[$i]->[0];
                my $out = $outs[$i]->[1];
                $out->close();
                my $split_file = $outs[$i]->[2];
                unlink($split_file);
                unless ($fastq_file =~ /\.gz$/) {
                    $split_file =~ s/\.gz$//;
                }
                symlink(abs_path($fastq_file), $split_file);
            }
            return 1;
        }
    }
    
    my $num_bases = 0;
    my $expected_lines = @ins * 4;
    my $count = 0;
    while (1) {
        # get the next entry (4 lines) from each input fastq
        my @seqs;
        my $lines = 0;
        my $these_bases = 0;
        foreach my $i (0..$#ins) {
            my $in_fh = $ins[$i]->[1]->fh;
            
            for (1..4) {
                my $line = <$in_fh>;
                defined $line || next;
                $lines++;
                
                push(@{$seqs[$i]}, $line);
                
                if ($_ == 2) {
                    my $seq = $seqs[$i]->[1];
                    chomp($seq);
                    $these_bases += length($seq);
                }
            }
        }
        $count++;
        
        # check for truncation/ eof
        if ($lines == 0) {
            last;
        }
        elsif ($lines != $expected_lines) {
            $self->throw("one of the fastq files ended early");
        }
        
        # start a new chunk if necessary
        $num_bases += $these_bases;
        if ($num_bases > $chunk_size) {
            $split_num++;
            $num_bases = $these_bases;
            
            foreach my $ref (@outs) {
                my ($prefix, $old) = @{$ref};
                $old->close;
                
                my $split_file = $fsu->catfile($split_dir, "$prefix.$split_num.fastq.gz");
                my $out = VertRes::IO->new(file => ">$split_file", verbose => $self->verbose);
                $ref->[1] = $out;
            }
        }
        
        # print out the entries
        foreach my $i (0..$#seqs) {
            my @lines = @{$seqs[$i]};
            my $out_fh = $outs[$i]->[1]->fh;
            foreach (@lines) {
                print $out_fh $_;
            }
        }
    }
    foreach my $ref (@ins, @outs) {
        $ref->[1]->close;
    }
    
    # check the chunks seem fine
    foreach my $i (0..$#ins) {
        my $fastq_file = $ins[$i]->[0];
        $io->file($fastq_file);
        my $in_lines = $io->num_lines;
        
        my $prefix = $outs[$i]->[0];
        my $out_lines = 0;
        foreach my $test_split_num (1..$split_num) {
            my $split_file = $fsu->catfile($split_dir, "$prefix.$test_split_num.fastq.gz");
            $io->file($split_file);
            $out_lines += $io->num_lines;
        }
        
        unless ($out_lines == $in_lines) {
            $self->throw("$fastq_file had $in_lines lines, but the split files ended up with only $out_lines!");
        }
    }
    
    return $split_num;
}

=head2 filter_reads

 Title   : filter_reads
 Usage   : $obj->filter_reads($fastq_file, $filt_file, min_length => 30);
 Function: Filter a fastq file so that certain sequences are excluded.
 Returns : int (the number of sequences output)
 Args    : input fastq (can be gzip compressed), output fastq (can be
           automatically compressed if it is named ...gz), filtering options
           hash:
           min_length => int (only sequences this length or longer are output)

=cut

sub filter_reads {
    my ($self, $in, $out, %filt) = @_;
    
    my $i = VertRes::IO->new(file => $in);
    my $ifh = $i->fh;
    my $o = VertRes::IO->new(file => ">$out");
    my $ofh = $o->fh;
    
    my $min_length = $filt{min_length};
    
    # we only support one option right now, so die if it wasn't supplied
    $min_length || $self->throw("Since min_length is the only filtering option right now, it is required");
    
    my $count = 0;
    while (<$ifh>) {
        chomp;
        my $name = $_;
        
        my $seq = <$ifh>;
        chomp($seq);
        
        my $qname = <$ifh>;
        chomp($qname);
        
        my $quals = <$ifh>;
        chomp($quals);
        
        unless ($name && $seq && $qname && $quals) {
            $self->throw("Fastq file '$in' is bad: truncated sequence entry");
        }
        
        if (defined $min_length && length($seq) >= $min_length) {
            $count++;
            print $ofh "$name\n$seq\n$qname\n$quals\n";
        }
    }
    close($ifh);
    close($ofh);
    
    return $count;
}

=head2 clip_point

 Title   : clip_point
 Usage   : my $clip_point = $obj->clip_point($fastq_file);
 Function: Find the base position where 90% of reads have fewer than 2 Ns.
           This can be used as a suitable hard clipping point for MAQ mapping.
 Returns : int (-1 if no clipping is necessary or could satisfy the 90%)
 Args    : input fastq (can be gzip compressed)

=cut

sub clip_point {
    my ($self, $fastq) = @_;
    
    my $fp = VertRes::Parser::fastq->new(file => $fastq);
    my $rh = $fp->result_holder;
    
    my $totalReads = 0;
    my %clipLengths;
    while ($fp->next_result) {
        my $seq = uc($rh->[1]);
        
        # we'll look for the position of the second N in the sequence
        my $ns = $seq =~ tr/N//;
        if ($ns >= 2) {
            my @s = split(//, $seq);
            my $nCount = 0;
            my $clipLength = @s;
            foreach my $i (0..$#s) {
                if ($s[$i] eq 'N') {
                    $nCount++;
                    
                    if ($nCount == 2) {
                        # *** shouldn't this be $i + 1? Keeping it as it was
                        #     originally...
                        $clipLengths{$i}++;
                        last;
                    }
                }
            }
        }
        
        $totalReads++;
    }
    
    # determine the point where 90% of the reads are ok
    my $numReads = 0;
    foreach (sort { $a <=> $b } keys %clipLengths) {
        # *** does this make sense? Could investigate...
        $numReads += $clipLengths{$_};
        
        if ($numReads > $totalReads * 0.9) {
            return ($_ - 1);
        }
    }
    
    return -1;
}

use Inline C => <<'END_C';

=head2 qual_to_ints

 Title   : qual_to_ints
 Usage   : my @qualities = $obj->qual_to_ints($quality_string);
 Function: Convert the quality string of a fastq sequence into quality integers.
           NB: this currently only works correctly for sanger (phred) quality
           strings, as found in sam files.
 Returns : list of int
 Args    : quality string

=cut

void qual_to_ints(SV* obj, char* str) {
    Inline_Stack_Vars;
    
    Inline_Stack_Reset;
    
    char c;
    int i = 0;
    while (c = str[i++]) {
        Inline_Stack_Push(sv_2mortal(newSViv(c - 33)));
    }
    
    Inline_Stack_Done;
}

END_C

1;
