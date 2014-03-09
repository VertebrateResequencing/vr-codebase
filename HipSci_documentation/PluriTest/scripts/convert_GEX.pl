#!/bin/perl -w

use strict;

if($#ARGV != 1){
    print "Usage: $0 <ANNOTATION FILE> <PROFILE FILE>\n";
    exit 0;
}

my(
    $i, $j, $k, 
    @VALS, 
    $ANNOS_FILE, 
    $PROFILE_FILE, 
    $IDOne,
    $IDTwo,
    $IDThree, 
    $annoCount, 
    $found, 
    @arr
);

$ANNOS_FILE = $ARGV[0];
$PROFILE_FILE = $ARGV[1];

$i = 0;

open (ANNOS, "$ANNOS_FILE") || die "cannot open $ANNOS_FILE ($!)";
<ANNOS>;
while(<ANNOS>){
    chomp;
    @VALS = split /\t/;

    $IDOne = $VALS[0];
    $IDTwo = $VALS[1];
    $IDThree = $VALS[4];

    $arr[$i][0] = $IDOne;
    $arr[$i][1] = $IDTwo;
    $arr[$i++][2] = $IDThree;
}
close (ANNOS);

$annoCount = $i;

#for($j = 0; $j < $i; $j++){
#    print "$arr[$j][0]\t$arr[$j][1]\t$arr[$j][2]\n";
#}

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

print "$VALS[1]\t";
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
