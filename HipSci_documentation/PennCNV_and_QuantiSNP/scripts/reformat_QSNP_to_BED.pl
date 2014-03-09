#!/bin/perl -w

use strict;

if($#ARGV != 0){
    print "Usage: $0 <QSNP DIR>\n";
    exit 0;
}

my(
    @CNVInfo, 
    $QUANTISNP_DIR, 
    $CNV_FILE, 
    $BED_FILE, 
    $CNVChr, 
    $CNVCoordOne, 
    $CNVCoordTwo, 
    @inFiles
);

$QUANTISNP_DIR = $ARGV[0];

@inFiles = <$QUANTISNP_DIR/*.cnv>;

foreach(@inFiles){
    $CNV_FILE = $_; 
    /^(.*)\.cnv/;
    $BED_FILE = $1 . '.bed';
    open(CNV_FILE, "$CNV_FILE") || die "cannot open $CNV_FILE ($!)";
    open(BED_FILE, ">$BED_FILE") || die "cannot open $BED_FILE ($!)";
    <CNV_FILE>;
    while(<CNV_FILE>){
        chomp;
        @CNVInfo = split /\s+/;
        $CNVChr = $CNVInfo[1];
        $CNVCoordOne = $CNVInfo[2];
        $CNVCoordTwo = $CNVInfo[3];
        print BED_FILE "$CNVChr\t$CNVCoordOne\t$CNVCoordTwo\n";
    }
    close(BED_FILE);
    close(CNV_FILE);
    print "Converted $CNV_FILE to $BED_FILE\n";
}

