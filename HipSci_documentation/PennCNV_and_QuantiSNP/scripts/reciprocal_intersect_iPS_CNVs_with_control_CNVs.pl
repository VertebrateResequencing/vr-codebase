#!/bin/perl -w

use strict;

if($#ARGV != 1){
    print "Usage: $0 <CNV BED FILE> <CONTROL CNV BED FILE>\n";
    exit 0;
}

my(
    $CNV_FILE, 
    $CONTROL_CNV_FILE, 
    $CNVCount, 
    $controlOverlapCount, 
    $noOverlapCNVs, 
    $cutOff
);

$CNV_FILE = $ARGV[0];
$CONTROL_CNV_FILE = $ARGV[1];

$cutOff = 0.5;

#bedtools intersect -u -f 0.5 -r -a 283163_A02_qc1hip5533821-batch02-02.bed_10_130_10 -b 283163_D02_qc1hip5533824-batch02-11.bed_10_130_10 | wc -l

$CNVCount = 0;
open(CNV_FILE, "$CNV_FILE") || die "cannot open $CNV_FILE ($!)";
while(<CNV_FILE>){
    $CNVCount++;
}
close(CNV_FILE);

$controlOverlapCount = 0;
chomp($controlOverlapCount = `bedtools intersect -u -f $cutOff -r -a $CNV_FILE -b $CONTROL_CNV_FILE | wc -l`);
$noOverlapCNVs = $CNVCount - $controlOverlapCount;

print "Number of CNVs:            $CNVCount\n";
print "Overlap with control CNVs: $controlOverlapCount\n";
print "CNVs minus control:        $noOverlapCNVs\n";

