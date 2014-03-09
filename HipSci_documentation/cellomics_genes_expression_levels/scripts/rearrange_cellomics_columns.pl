#!/bin/perl -w

use strict;

if($#ARGV != 0){
    print "Usage: $0 <CELLOMICS GENES OUTPUT FILE>\n";
    exit 0;
}

# LINE 40 MUST BE EDITED FOR CORRECT REARRANGEMENT - SEE DOCUMENTATION FOR FURTHER INFO

my(
    $i, 
    $thisNumberOfValuesCount, 
    $expressionValuesCount, 
    @profileValues, 
    $inputFile, 
    $lineNumber
    );

$inputFile = $ARGV[0];
$lineNumber = 0;
$expressionValuesCount = 0;

open (PROFILE_FILE, "$inputFile") || die "cannot open $_ ($!)";
while(<PROFILE_FILE>){
    $lineNumber++;
    chomp;
    @profileValues = split /\t/;

    my $it = 0;
    $thisNumberOfValuesCount = $#profileValues;
    if($lineNumber == 1){
        $expressionValuesCount = $thisNumberOfValuesCount;
    }
    if($thisNumberOfValuesCount != $expressionValuesCount){
        print "$_\n";
    }
    else{
        my @order = (6, 11, 2, 4, 8, 12, 10, 7, 5, 1, 9, 3); # FOR SET 14
        my $listCount = 0;
        print "$profileValues[0]\t";
        foreach(@order){
            $i = $_;
            $listCount++;
            print "$profileValues[$i]\t";
            if($listCount <= $#order){
                print "\t";
            }
        }
        print "\n";        
    }
}
close (PROFILE_FILE);
