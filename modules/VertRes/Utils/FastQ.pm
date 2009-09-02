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
use File::Basename;
use VertRes::Parser::fastqcheck;
use Cwd 'abs_path';

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
        
        my $split_file = $io->catfile($split_dir, "$prefix.$split_num.fastq.gz");
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
                
                my $split_file = $io->catfile($split_dir, "$prefix.$split_num.fastq.gz");
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
            my $split_file = $io->catfile($split_dir, "$prefix.$test_split_num.fastq.gz");
            $io->file($split_file);
            $out_lines += $io->num_lines;
        }
        
        unless ($out_lines == $in_lines) {
            $self->throw("$fastq_file had $in_lines lines, but the split files ended up with only $out_lines!");
        }
    }
    
    return $split_num;
}


use Inline C => DATA => FILTERS => 'Strip_POD';

1;

__DATA__
__C__

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
