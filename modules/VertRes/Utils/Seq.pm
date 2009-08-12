=head1 NAME

VertRes::Utils::Seq - sequence utility functions

=head1 SYNOPSIS

use VertRes::Utils::Seq;

my $seq_util = VertRes::Utils::Seq->new();

$seq_util->rev_com('atgc');

=head1 DESCRIPTION

General utility functions for working on or with sequence strings and files.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Utils::Seq;

use strict;
use warnings;

use base qw(VertRes::Base);


=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Utils::Seq->new();
 Function: Create a new VertRes::Utils::Seq object.
 Returns : VertRes::Utils::Seq object
 Args    : n/a

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    return $self;
}

=head2 rev_com

 Title   : rev_com
 Usage   : my $rev_com = $obj->rev_com('atgc');
 Function: Reverse compliment a nucleotide string.
 Returns : string
 Args    : string

=cut

sub rev_com {
    my ($self, $seq) = @_;
    
    $seq = uc $seq;
    $seq = reverse $seq;
    $seq =~ tr/ACGTN/TGCAN/;
    
    return $seq;
}

1;
