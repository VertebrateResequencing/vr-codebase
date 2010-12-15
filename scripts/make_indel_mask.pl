#!/usr/bin/env perl
#
# usage: make_indel_mask.pl indel.calls.bed mask_size indels.mask.bed
#
# Author: Sendu Bala <bix@sendu.me.uk>

use strict;
use warnings;

use File::Spec;
use File::Basename;

# get user input
my ($in_bed, $mask_size, $out_bed) = @ARGV;
chomp($out_bed);

($in_bed && -e $in_bed && $mask_size && $out_bed) or die <<USAGE;
Creates a mask version of a bed file.

Usage: $0 indel.calls.bed mask_size indels.mask.bed

where indel.calls.bed is the input bed file, mask_size is the size in bases
(eg. 10), and indels.mask.bed is the output file

USAGE

open(my $ifh, $in_bed) || die "Could not open $in_bed\n";
open(my $ofh, '>', $out_bed) || die "Could not write to $out_bed\n";
$mask_size = int($mask_size);
$mask_size >= 1 || die "Invalid mask_size '$mask_size' (should be int >= 1)\n";


# covert to mask
while (<$ifh>) {
    # chr1    71996   72005   -AAAAAAAAA:4/6
    my ($chr, $start, $stop, $notes) = split;
    
    $start -= $mask_size;
    $stop += $mask_size;
    $notes .= ":+/-$mask_size" ;
    
    print $ofh join("\t", $chr, $start, $stop, $notes), "\n";
}
close($ifh);
close($ofh);

exit;
