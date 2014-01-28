#!/usr/bin/env perl

use strict;
use warnings;
no warnings 'uninitialized';

if($#ARGV != 3){
    print "Usage: $0 <PARSE> <ANNO> <ID MAPPING> <FIRST N SAMPLES>\n";
    exit 0;
}

my(
    $i, $j, $k, 
    @VALS, 
    $ANNOS_FILE, 
    $PROFILE_FILE, 
    $MAPPING_FILE, 
    $IDOne,
    $IDTwo,
    $IDThree, 
    $annoCount, 
    $found, 
    @arr, 
    @IDMap, 
    $numberOfSamples
);

$PROFILE_FILE = $ARGV[0];
$ANNOS_FILE = $ARGV[1];
$MAPPING_FILE = $ARGV[2];
$numberOfSamples = $ARGV[3];

$i = 0;

open (ANNOS, "$ANNOS_FILE") || die "cannot open $ANNOS_FILE ($!)";
<ANNOS>;
while(<ANNOS>){
    chomp;
    @VALS = split /\t/;

    $arr[$i][0] = $VALS[0]; # TargetID e.g. 7A5
    $arr[$i][1] = $VALS[1]; # ProbeID e.g. 6450255
    $arr[$i++][2] = $VALS[4]; # PROBE_ID e.g. ILMN_1762337
}
close (ANNOS);

$annoCount = $i;

$i = 0;
open (IDS, "$MAPPING_FILE") || die "cannot open $MAPPING_FILE ($!)";
while(<IDS>){
    #print;
    chomp;
    @VALS = split /\t/;
    $IDMap[$i][0] = $VALS[0]; # A
    $_ = $VALS[1]; 
    s/ //g; # RV Bob -> RVBob
    $IDMap[$i++][1] = $_; # RVBob
}
close (IDS);

open (PROFILE, "$PROFILE_FILE") || die "cannot open $PROFILE_FILE ($!)";
$_ = <PROFILE>;
chomp;
@VALS = split /\t/;

my $valsCount = $#VALS;
my $sampleCount = (($valsCount - 2) + 1) / 8;

my $indexStart;
my $positionOne;
my $positionTwo;
my $positionThree;
my $positionFour;

#1
#2
#3#MIN_Signal-9257622034_A	
#4#AVG_Signal-9257622034_A	
#5#MAX_Signal-9257622034_A	
#6#NARRAYS-9257622034_A	
#7#ARRAY_STDEV-9257622034_ABEAD_STDEV-9257622034_A	
#8#Avg_NBEADS-9257622034_A	
#9#Detection-9257622034_A
#3  7  8  9
#11 15 16 17

#PRINT OUT THE HEADER:

$sampleCount = $numberOfSamples;

print "$VALS[1]\t";
for($k = 0; $k < $sampleCount; $k++){
    $indexStart = ($k * 8) + 2;
    $positionOne = ($indexStart + 1);
    $positionTwo = ($indexStart + 5);
    $positionThree = ($indexStart + 6);
    $positionFour = ($indexStart + 7);
    my $headerOne = $VALS[$positionOne];
    my $headerTwo = $VALS[$positionTwo];
    my $headerThree = $VALS[$positionThree];
    my $headerFour = $VALS[$positionFour];

    $_ = $headerFour;
    /^(.*?)-(.*)/;
    if($1 eq "Detection"){
        $headerFour = 'Detection.Pval' . '-' . "$2";
    }

    my $partOne;
    my $partTwo;
    my $partThree;

    $_ = $headerOne;
    /^(.*?)-(.*?)_(.*)/;
    $partOne = $1;
    $partTwo = $2;
    $partThree = $3;
    for($i = 0; $i <= $#IDMap; $i++){
        if($IDMap[$i][0] eq $partThree){
            $partThree = $IDMap[$i][1];
        }
    }
    $headerOne = $partOne . '-' . $partTwo . '_' . $partThree;

    $_ = $headerTwo;
    /^(.*?)-(.*?)_(.*)/;
    $partOne = $1;
    $partTwo = $2;
    $partThree = $3;
    for($i = 0; $i <= $#IDMap; $i++){
        if($IDMap[$i][0] eq $partThree){
            $partThree = $IDMap[$i][1];
        }
    }
    $headerTwo = $partOne . '-' . $partTwo . '_' . $partThree;

    $_ = $headerThree;
    /^(.*?)-(.*?)_(.*)/;
    $partOne = $1;
    $partTwo = $2;
    $partThree = $3;
    for($i = 0; $i <= $#IDMap; $i++){
        if($IDMap[$i][0] eq $partThree){
            $partThree = $IDMap[$i][1];
        }
    }
    $headerThree = $partOne . '-' . $partTwo . '_' . $partThree;

    $_ = $headerFour;
    /^(.*?)-(.*?)_(.*)/;
    $partOne = $1;
    $partTwo = $2;
    $partThree = $3;
    for($i = 0; $i <= $#IDMap; $i++){
        if($IDMap[$i][0] eq $partThree){
            $partThree = $IDMap[$i][1];
        }
    }
    $headerFour = $partOne . '-' . $partTwo . '_' . $partThree;

    print "$headerOne\t$headerTwo\t$headerThree\t$headerFour";
    if($k == ($sampleCount - 1)){
        print "\n";
    }
    else{
        print "\t";
    }
}

#PRINT OUT THE VALUES:

while(<PROFILE>){
    chomp;
    @VALS = split /\t/;

    $IDOne = $VALS[0]; # profile ID #1 = TargetID
    $IDTwo = $VALS[1]; # profile ID #2 = ProbeID

    $found = 0;
    for($j = 0; $j < $annoCount; $j++){
        if(($IDOne eq $arr[$j][0]) && ($IDTwo eq $arr[$j][1])){ #arr[0] = TargetID arr[1] = ProbeID arr[2] = PROBEID/ILMN_..
            print "$arr[$j][2]\t";
            $found = 1;
        }
    }
    if($found == 0){
        print "Couldn't find $IDOne\n";
        exit;
    }

    for($k = 0; $k < $sampleCount; $k++){
        $indexStart = ($k * 8) + 2;
        $positionOne = ($indexStart + 1);
        $positionTwo = ($indexStart + 5);
        $positionThree = ($indexStart + 6);
        $positionFour = ($indexStart + 7);
        print "$VALS[$positionOne]\t$VALS[$positionTwo]\t$VALS[$positionThree]\t$VALS[$positionFour]";
        if($k == ($sampleCount - 1)){
            print "\n";
        }
        else{
            print "\t";
        }
    }
}
close (PROFILE);
