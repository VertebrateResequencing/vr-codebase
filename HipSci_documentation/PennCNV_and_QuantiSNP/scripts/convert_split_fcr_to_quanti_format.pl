#!/bin/perl -w

use strict;

if($#ARGV != 0){
    print "Usage: $0 <SPLIT FCR FILE>\n";
    exit 0;
}

#EXAMPLE USAGE:
#perl convert_split_fcr_to_quanti_format.pl split_files/271298_A01_hipscigt5466706 HumanCoreExome-12v1-0_A.csv

my(
    $i, 
    $manifestInfoCount, 
    @tmpArr, 
    $splitGTFile, 
    $manifestFile, 
    @manifestSNPInfo, 
    $found, 
    $chr, 
    $position, 
    $name, 
    $bAlleleFreq, 
    $logRRatio
);

$splitGTFile = $ARGV[0];
$manifestFile = '/lustre/scratch105/vrpipe/refs/hipsci/resources/genotyping/archive/HumanCoreExome-12v1-0_A.csv';

$manifestInfoCount = 0;
open(MANIFEST_FILE, "$manifestFile") || die "cannot open $manifestFile ($!)";
<MANIFEST_FILE>;
while(<MANIFEST_FILE>){
    chomp;
    @tmpArr = split /,/;
    $manifestSNPInfo[$manifestInfoCount][0] = $tmpArr[1];
    $manifestSNPInfo[$manifestInfoCount][1] = $tmpArr[9];
    $manifestSNPInfo[$manifestInfoCount][2] = $tmpArr[10];
    $manifestInfoCount++;
}
close(MANIFEST_FILE);

open(SPLIT_FCR_FILE, "$splitGTFile") || die "cannot open $splitGTFile ($!)";
$_ = <SPLIT_FCR_FILE>;
print "Name\tSample ID\tChr\tPosition\tB Allele Freq\tLog R Ratio\n";
while(<SPLIT_FCR_FILE>){
    chomp;
    @tmpArr = split /\t/;
    $name = $bAlleleFreq = $logRRatio = $chr = $position = '';
    $name = $tmpArr[0];
    $bAlleleFreq = $tmpArr[11];
    $logRRatio = $tmpArr[12];
    $found = 0;
    for($i = 0; $i < $manifestInfoCount; $i++){
        if($manifestSNPInfo[$i][0] eq $name){
            $chr = $manifestSNPInfo[$i][1];
            $position = $manifestSNPInfo[$i][2];
            $found = 1;
            last;
        }
    }
    if($found == 0){
        print "#Couldn't find $name\n";
        exit;
    }
    if($chr eq 'X'){
        $chr = 23;
    }
    elsif($chr eq 'Y'){
        $chr = 24;
    }
    elsif($chr eq 'XY'){
        $chr = 23;
        next;
    }
    elsif($chr eq 'MT'){
        $chr = 'M';
        next;
    }
    print "$name\t$chr\t$position\t$bAlleleFreq\t$logRRatio\n";
}
close(SPLIT_FCR_FILE);
