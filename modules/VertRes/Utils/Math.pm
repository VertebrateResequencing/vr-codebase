=head1 NAME

VertRes::Utils::Math - math utility functions

=head1 SYNOPSIS

use VertRes::Utils::Math;

my $math_util = VertRes::Utils::Math->new();

my $median = $math_util->histogram_median({1 => 5, 2 => 10, 3 => 5});

=head1 DESCRIPTION

General utility functions for doing math/stats stuff.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Utils::Math;

use strict;
use warnings;

use base qw(VertRes::Base);


=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Utils::Math->new();
 Function: Create a new VertRes::Utils::Math object.
 Returns : VertRes::Utils::Math object
 Args    : n/a

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    return $self;
}

=head2 histogram_median

 Title   : histogram_median
 Usage   : my $median = $obj->histogram_median($histogram_hash_ref);
 Function: Get the median bin of a histogram.
 Returns : int
 Args    : hash ref where keys are bins and values are frequencies
           NB: assumes equal sized bins

=cut

sub histogram_median {
    my ($self, $hash) = @_;
    
    # find the half-way frequency
    my $total = 0;
    foreach my $freq (values %{$hash}) {
        $total += $freq;
    }
    my $half = sprintf("%0.0f", $total / 2);
    
    # find the corresponding bin
    my $median = 0; 
    my $current = 0;
    foreach my $bin (sort { $a <=> $b } keys %{$hash}) {
        $current += $hash->{$bin};
        if ($current >= $half) {
            $median = $bin;
            last;
        }
    }
    
    return $median;
}

1;
