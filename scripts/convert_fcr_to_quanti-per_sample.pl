#!/bin/perl -w

use strict;

if($#ARGV != 1){
    print "Usage: $0 <SPLIT GT (FCR) FILE> <MANIFEST FILE>\n";
    exit 0;
}

#EXAMPLE USAGE:
#perl convert_fcr_to_quanti-per_sample.pl split_files/271298_A01_hipscigt5466706 HumanCoreExome-12v1-0_A.csv
my(
    $i, 
    @tmpArr, 
    $splitGTFile, 
    $manifestFile, 
    $sampleName, 
    %manifestSNPInfo, 
    @chrArr,
    $found, 
    $chr, 
    $position, 
    $name, 
    $bAlleleFreq, 
    $logRRatio,
    $linenumber,
);

$splitGTFile = $ARGV[0];
$manifestFile =  $ARGV[1];
$linenumber = 0;

open(MANIFEST_FILE, "$manifestFile") || die "cannot open $manifestFile ($!)";
<MANIFEST_FILE>;
while(<MANIFEST_FILE>){
    chomp;
    @tmpArr = split /,/;
    if ( $tmpArr[9] && $tmpArr[10] ) {
		push @{ $manifestSNPInfo{ $tmpArr[1] } }, $tmpArr[9];
		push @{ $manifestSNPInfo{ $tmpArr[1] } }, $tmpArr[10];
	}
}
close(MANIFEST_FILE);

open(SPLIT_FCR_FILE, "$splitGTFile") || die "cannot open $splitGTFile ($!)";
$linenumber = -1;
print "Name\tChr\tPosition\tB Allele Freq\tLog R Ratio\n";
while(<SPLIT_FCR_FILE>){
	$linenumber++;
	next if $linenumber == 0;
    chomp;
    @tmpArr = split /\t/;
    $name = $bAlleleFreq = $logRRatio = $chr = $position = '';
    $name = $tmpArr[0];
    $bAlleleFreq = $tmpArr[11];
    $logRRatio = $tmpArr[12];
    $found = 0;
    if( $manifestSNPInfo{$name} ){
	    @chrArr = @{ $manifestSNPInfo{$name} };
		$chr = $chrArr[0];
		$position = $chrArr[1];
		$found = 1;
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
