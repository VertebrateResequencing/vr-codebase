#!/bin/perl -w

use strict;

if($#ARGV != 0){
    print "Usage: $0 <PENNCNV DIR>\n";
    exit 0;
}

my(
    @CNVInfo, 
    $PENNCNV_DIR, 
    $CNV_FILE, 
    $BED_FILE, 
    $CNVChr, 
    $CNVCoordOne, 
    $CNVCoordTwo, 
    @inFiles
);

$PENNCNV_DIR = $ARGV[0];

@inFiles = <$PENNCNV_DIR/*.filtered>;

foreach(@inFiles){
    $CNV_FILE = $_; 
    /^(.*)\.filtered/;
    $BED_FILE = $1 . '.bed';
    open(CNV_FILE, "$CNV_FILE") || die "cannot open $CNV_FILE ($!)";
    open(BED_FILE, ">$BED_FILE") || die "cannot open $BED_FILE ($!)";
    while(<CNV_FILE>){
        chomp;
        @CNVInfo = split /\s+/;
        $_ = $CNVInfo[0];
        if(!/^chr(.*)\:([0-9]+)\-([0-9]+)/){
            print "Problem with the format of $_\n";
            exit;
        }
        $_ = $1;
        /^(.*)\:(.*)-(.*)/;
        $CNVChr = $1;
        $CNVCoordOne = $2;
        $CNVCoordTwo = $3;
        print BED_FILE "$CNVChr\t$CNVCoordOne\t$CNVCoordTwo\n";
    }
    close(BED_FILE);
    close(CNV_FILE);
    print "Converted $CNV_FILE to $BED_FILE\n";
}
