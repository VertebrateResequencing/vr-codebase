=head1 NAME

VertRes::Utils::Sam - sam utility functions

=head1 SYNOPSIS

use VertRes::Utils::Sam;

my $sam_util = VertRes::Utils::Sam->new();

$sam_util->diff();

=head1 DESCRIPTION

General utility functions for working on or with sam/bam files.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Utils::Sam;

use strict;
use warnings;
use VertRes::IO;
use VertRes::Wrapper::samtools;

use base qw(VertRes::Base);


=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Utils::Sam->new();
 Function: Create a new VertRes::Utils::Sam object.
 Returns : VertRes::Utils::Sam object
 Args    : n/a

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    return $self;
}

=head2 bams_are_similar

 Title   : bams_are_similar
 Usage   : my $are_similar = $obj->bams_are_similar(@bam_files);
 Function: Check if multiple name-sorted bam files are the same, ignoring
           mapping quality 0 reads.
 Returns : boolean (false if they differ)
 Args    : list of bam files

=cut

sub bams_are_similar {
    my ($self, @files) = @_;
    
    @files >= 2 or $self->throw("At least 2 bam files are needed");
    
    my $sw = VertRes::Wrapper::samtools->new(verbose => $self->verbose,
                                             run_method => 'open',
                                             quiet => 1);
    
    # open
    my @in_fhs;
    foreach my $file (@files) {
        my $bamfh = $sw->view($file, undef, h => 1);
        push(@in_fhs, $bamfh);
    }
    
    my $first_fh = $in_fhs[0];
    while (! eof($first_fh)) {
        # if a line in one of the input bams has mapping quality 0, we skip it
        # along with the corresponding line in all other input bams, even if
        # they have better mapping qualities
        my @lines;
        while (1) {
            @lines = ();
            my $good_mapping_qs = 0;
            foreach my $fh (@in_fhs) {
                my $line = <$fh>;
                push(@lines, $line);
                
                if ($line =~ /^\@/) {
                    $good_mapping_qs++;
                    next;
                }
                
                my (undef, undef, undef, undef, $mapping_q) = split("\t", $line);
                $good_mapping_qs++ if $mapping_q > 0;
            }
            
            last if $good_mapping_qs == @in_fhs;
            last if eof($first_fh);
        }
        
        my %lines;
        foreach my $line (@lines) {
            if ($line =~ /^\@/) {
                $lines{$line}++;
                next;
            }
            
            my ($col1, $col2, $col3, $col4, $mapping_q) = split("\t", $line);
            $lines{"$col1\t[...]\t$col3\t$col4\t[...]"}++;
        }
        
        unless (keys %lines == 1) {
            $self->debug("Files started differing here:\n".join("\n", keys %lines));
            return 0;
        }
        
        my ($line, $count) = each %lines;
        return 0 unless $count == @in_fhs;
    }
    
    foreach my $fh (@in_fhs) {
        return 0 unless eof($fh);
    }
    
    return 1;
}

1;
