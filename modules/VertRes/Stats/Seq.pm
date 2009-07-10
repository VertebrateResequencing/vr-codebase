# perl module to calculate various things about sequences

=head1 NAME

SeqStats - calculate various things about sequences

=head1 SYNOPSIS

use SeqStats;

# create object
my $seq_stats_obj = SeqStats->new();

# find out something about a sequence string
my $gc_percentage = $seq_stats_obj->gc_percent('atgc');

=head1 DESCRIPTION

A module that can contain various methods related to sequence stats. Currently
implements a single method for calculating the GC percentage of a DNA sequence.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package SeqStats;

use strict;
use warnings;
use Carp;


=head2 new

 Title   : new
 Usage   : my $obj = SeqStats->new()
 Function: Build a new SeqStats object
 Returns : SeqStats object
 Args    : n/a

=cut

sub new {
    my $class = shift;
    
    $class = ref($class) || $class;
    my $self = {};
    bless($self, $class);
    
    return $self;
}

=head2 gc_percent

 Title   : gc_percent
 Usage   : my $gc = $obj->gc_percent('atgc');
 Function: calculate the GC percentage of a DNA sequence
 Returns : real number (the percent)
 Args    : DNA sequence string, optionally the number of decimal places to
           calculate percentage to (default 2)

=cut

sub gc_percent {
    my ($self, $dna, $places) = @_;
    $dna || return 0;
    unless (defined $places) {
        $places = 2;
    }
    
    # we'll just assume at this point our $dna is only atgc; if it might not be
    # we could add sequence checking etc. to this module
    
    # count the gc quickly
    my $gc_count = $dna =~ tr/gcGC//;
    
    return sprintf("%0.${places}f", ((100 / length($dna)) * $gc_count));
}

1;
