#!/usr/bin/env perl

use warnings;
use strict;

if($#ARGV != 0){
    print "Usage: $0 <VCF FILE>\n";
    exit 0;
}

###################################################
# THIS IS STILL IN DEVELOPMENT - VERSION 0.1
#
# PROGRAM FOR IDENTIFYING LOH REGIONS USING VCFs
# GENERATED FROM GENOTYPING DATA FOR HipSci PROJECT
# THIS IS NOT INTENDED AS A GENERIC CALLER
# THIS PROGRAM IS USED TO FIND LOH REGIONS 
# COMPARED TO FIBORBLAST/CONTROL SAMPLE THAT 
# A HipSci iPS SAMPLE HAS BEEN DERIVED FROM
# THROUGH GENETIC REPROGRAMMING
#
###################################################

my(
    $i, $j, $k, $l, 
    $VCFFile, 
    @tmpArr, 
    @sampleNames, 
    @AoA,

    $totalEntryCount, 
    $printValue, 
    $IAIndex, 
    $IBIndex, 
    $thisIndexNumber, 
    $average, 
    $rounded, 

    @allGTMismatchesPerSample, 
    @IAsAtAllGTMismatches, 
    @IBsAtAllGTMismatches, 
    @totalIAs, 
    @totalIBs, 

    $numberOfSamples, 
    $sampleNameLength, 

    $windowSize, 
    $percentageFlexibility, 

    $minRequiredInWindow, 
    $flexibilitySubtractionFromWindow, 
    $flexibilitySubtractionFromWindowRoundedDown, 

    @totalMismatchesPerRow, 
    $totalMismatchesPerRowCount, 
    $lastChr, 
    $lastCoord, 

    @lastInBlockArr, 
    @inBlockArr, 
    @totalBlocks, 

    $thisChr, 
    $thisCoord, 
    $nextChr, 
    $nextCoord, 

    @lastBlockLabel, 

    @blocksAoA
);

##

$VCFFile               = $ARGV[0];
#$windowSize            = $ARGV[1];
#$percentageFlexibility = $ARGV[2];
$windowSize            = 10;
$percentageFlexibility = 80;

print "VCF file                     = $VCFFile\n";
print "Number of SNPs in CNV window = $windowSize\n";
print "Percentage flexibility       = $percentageFlexibility\n";

$minRequiredInWindow = 0;
$flexibilitySubtractionFromWindow = 0;

$flexibilitySubtractionFromWindow            = ($windowSize / 100) * $percentageFlexibility;
$flexibilitySubtractionFromWindowRoundedDown = int $flexibilitySubtractionFromWindow;
$minRequiredInWindow                         = $windowSize - $flexibilitySubtractionFromWindowRoundedDown;

#print "Summary:\n";
#print "flexibilitySubtractionFromWindow:$flexibilitySubtractionFromWindow\n";
#print "roundedDown                     :$roundedDown\n";
#print "minRequiredInWindow:$minRequiredInWindow\n";

if($minRequiredInWindow < 2){
    print "percentageFlexibility is too high, need a minium of 2 SNPs in a window.\n";
    exit;
}

## Skip the header

open(VCF, "$VCFFile") || die "cannot open $VCFFile ($!)";
while (<VCF>){
    chomp;
    if(/^#CHROM/){
        last;
    }
}

# GET THE SAMPLE NAMES FROM THE HEADER LINE OF THE VCF
@tmpArr = split /\t/;
for($i = 9; $i <= $#tmpArr; $i++){
    push(@sampleNames, $tmpArr[$i]);
}
$numberOfSamples = $#sampleNames + 1;

# PRINT THE HEADER OF THE CNV OUTPUT

printCNVHeader();

#### END OF PRINT CNV HEADER

$lastChr = '0';
$lastCoord = '0';

# SET THE STORAGE VALUES TO 0

for($i = 0; $i < $numberOfSamples; $i++){
    $allGTMismatchesPerSample[$i] = 0;
    $IAsAtAllGTMismatches[$i] = 0;
    $IBsAtAllGTMismatches[$i] = 0;
    $totalIAs[$i] = 0;
    $totalIBs[$i] = 0;
    $totalMismatchesPerRow[$i] = 0;
}
$totalEntryCount = 0;

####

readVCFAndShowAllGTDifferences();

#NOW STORED:
#$allGTMismatchesPerSample[$i]  # total number of changes for each sample e.g. 24 2 6 etc.
#$IAsAtAllGTMismatches[$i]   # total  per sample 12.2 3.2 12.8 etc.
#$IBsAtAllGTMismatches[$i]
#push (@AoA, [@tmpArr]); # STORE ALL          LINES
#$totalMismatchesPerRow[$totalMismatchesPerRowCount] = ($totalMismatchesPerRow[$totalMismatchesPerRowCount] + 1); # store the number of mismatches e.g. 0:10 1:0 2:4 etc. = 10 rows have 0 mismatches, 0 rows have 1 mismatch, 4 rows have 2 mismatches.

#### PRINT SUMMARY OF ALL GT DIFFERENCES

printSummaryOfAllGTDifferences();

##########################################################################

# IDENTIFY FLEXIBLE BLOCKS AND STORE

# Create an empty @blocksAoA[][]
for $i (0 .. $#AoA){
    for($j = 0; $j < ($numberOfSamples - 1); $j++){
        $inBlockArr[$j] = 0;
    }
    push (@blocksAoA, [@inBlockArr]);
}
# Create an empty lastBlockLabel
for($j = 0; $j < ($numberOfSamples - 1); $j++){
    $lastBlockLabel[$j] = 0;
}

####

#printAoA();

identifyFlexibleBlocksAndStore();

#printFlexibleBlocksNumbers(); # why is 2nd block labelled 0 for FS03?
printFlexibleBlocks();

#######################

sub printFlexibleBlocksNumbers{
                                                                                                   
    my $blockFound;                                                                                                   
                                                                                                   
    for $i (0 .. ($#blocksAoA)){                                                                   
        $blockFound = 0;                                                                           
        for($j = 0; $j < ($numberOfSamples - 1); $j++){                                            
            if($blocksAoA[$i][$j] > 0){                                                            
                $blockFound++;                                                                     
            }                                                                                      
        }                                                                                          
        if($blockFound > 0){                                                                       
            print "i=$i block#=$blocksAoA[$i][0]; ";                                              
            print "$AoA[$i][0] $AoA[$i][1] $AoA[$i][2]; ";                                        
            for($j = 0; $j < ($numberOfSamples - 1); $j++){                                        
                print "$blocksAoA[$i][$j]; ";                                                      
            }                                                                                      
            print "\n";                                                                            
                                                                                                   
        }                                                                                          
    }                                                  
}

sub printFlexibleBlocks{

    my $blockFound;
    my @lastLocalBlockLabel;
    my @allTrimmedBlocks;
    my $start;
    my $CNVSummaryFile = 'CNVs_summary.txt';
    open(CNVSUMMARY, ">$CNVSummaryFile") || die "cannot open $CNVSummaryFile ($!)";

    for($j = 0; $j < ($numberOfSamples - 1); $j++){
        $lastLocalBlockLabel[$j] = 0; 
    }

    print "CNVs ($windowSize SNPs with $percentageFlexibility% flexibility) in $VCFFile\n";
    print "====================================================\n\n";

    # SCAN blocksAoA, STORE BLOCKS AS TRIMMED BLOCKS FOR SUBSEQUENT PRINTING
    for($j = 0; $j < ($numberOfSamples - 1); $j++){
        my @trimmedBlocks;
        my @countStore;
        my @zeroStore;
        my @mismatchStore;

        #SKIP 0's UP TO THE START OF THE FIRST BLOCK
        $start = 0;
        for $i (0 .. ($#blocksAoA)){
            if($blocksAoA[$i][$j] > 0){
                $start = $i;
                last;
            }
        }
        if($start > 0){
            for $i ($start .. ($#blocksAoA)){
                #0 FOUND - STORE IN ZERO STORE TEMPORARILY UNTIL NEXT NUMBER - IF NEXT NUMBER IS THE 
                #SAME THEN THE ZEROs WILL BE KEPT AND TRANFERRED TO COUNT STORE, IF IT'S DIFFERENT THEY'LL BE FORGOTTEN
                if($blocksAoA[$i][$j] == 0){
                    push(@zeroStore, $i);
                }
                #SAME AS LAST LABEL - TRANSFER STORED 0s TO COUNT STORE AND THEN NUMBER TO COUNTSTORE
                elsif($lastLocalBlockLabel[$j] == $blocksAoA[$i][$j]){
                    if($#zeroStore >= 0){
                        foreach(@zeroStore){
                            push(@countStore, $_);
                            push(@mismatchStore, 0);
                        }
                        $#zeroStore = -1;
                    }
                    push(@countStore, $i);
                    push(@mismatchStore, 1);
                }
                #DIFFERENT TO LAST NUMBER - TRANSFER COUNTSTORE TO TRIMMED BLOCKS FOR FINAL STORAGE
                #RESET COUNTSTORE AND ZEROSTORE, THEN STORE THE FIRST NUMBER IN COUNTSTORE
                elsif($lastLocalBlockLabel[$j] != $blocksAoA[$i][$j]){
                    $lastLocalBlockLabel[$j] = $blocksAoA[$i][$j];
                    if($#countStore >= 0){
                        push (@{$trimmedBlocks[0]}, [@countStore]);
                        push (@{$trimmedBlocks[1]}, [@mismatchStore]);
                    }
                    $#countStore = -1;
                    $#mismatchStore = -1;
                    $#zeroStore = -1;
                    push(@countStore, $i);
                    push(@mismatchStore, 1);
                }
                else{
                    print "UNEXPECTED OPTION\n";
                    exit;
                }
            }
            if($#countStore >= 0){
                push (@{$trimmedBlocks[0]}, [@countStore]);
                push (@{$trimmedBlocks[1]}, [@mismatchStore]);
                #print "Was pushed $#trimmedBlocks\n";
            }
        }
        push (@allTrimmedBlocks, [@trimmedBlocks]);
    }

    #@allTrimmedBlocks structure:
    #                      if 0 length, no CNVs    
    #[sample#][0=SNP overall index pos in AoA;1=mismatches][CNV#?][n for SNP index pos; 0/1 for mismatches]
    for $i (0 .. $#allTrimmedBlocks){
        my $printCNVSummaryString;
        my $aRef = $allTrimmedBlocks[$i];
        my $p = @$aRef - 1;
        print "Sample " . ($i+1) . ":\n";
        $printCNVSummaryString .= "$sampleNames[$i]:";
        if($p == -1){
            print "0 CNVs\n\n";
            $printCNVSummaryString .= "0 CNVs";
            if($i < $#allTrimmedBlocks){
                $printCNVSummaryString .= "; ";
            }
            print CNVSUMMARY $printCNVSummaryString;
        }
        if($p == 1){
            #for $j (0 .. $p){
            $j = 0; # because it's a fixed length i.e. always 2
                my $bRef = $allTrimmedBlocks[$i][$j];
                my $q = @$bRef - 1;
                if($q >= 0){
                    if($q == 0){
                        print "1 CNV\n";
                        $printCNVSummaryString .= "1 CNV (";
                    }
                    else{
                        print '' . ($q+1) . " CNVs\n";
                        $printCNVSummaryString .= ($q+1) . " CNVs (";
                    }
                    for $k (0 .. $q){
                        my $cRef = $allTrimmedBlocks[$i][$j][$k];
                        my $r = @$cRef - 1;
                        print "CNV #". ($k+1) . " = ";
                        if($r >= 0){
                            print '' . ($r+1) . " SNPs ";
                            $printCNVSummaryString .= ($r+1) . " SNPs:";
                            my $mismatchCount = 0;
                            for $l (0 .. $r){
                                if($allTrimmedBlocks[$i][1][$k][$l] == 1){
                                    $mismatchCount++;
                                }
                            }
                            print "$mismatchCount/" . ($r+1) . " mismatches\n";
                            for $l (0 .. $r){
                                my $localIndex = $allTrimmedBlocks[$i][0][$k][$l];
                                print "Chr:$AoA[$localIndex][0] Coord:$AoA[$localIndex][1] SNP:$AoA[$localIndex][2] Mismatch:";
                                if($l == 0){
                                    $printCNVSummaryString .= "Chr$AoA[$localIndex][0] Coords=$AoA[$localIndex][1]";
                                }
                                if($l == $r){
                                    $printCNVSummaryString .= "-$AoA[$localIndex][1]";
                                }
                                if($allTrimmedBlocks[$i][1][$k][$l] == 1){
                                    print "Y\n";
                                }
                                else{
                                    print "N\n";
                                }    
                            }
                            if($k < $q){
                                $printCNVSummaryString .= "; ";
                            }
                            print "\n";
                        }
                        else{
                            print "Shouldn't get here 1\n";
                            exit;
                        }
                    }
                    $printCNVSummaryString .= ")";
                    if($i < $#allTrimmedBlocks){
                        $printCNVSummaryString .= "; ";
                    }
                    print CNVSUMMARY "$printCNVSummaryString";
                }
                else{
                    print "Shouldn't get here 2\n";
                    exit;
                }
            #}
        }
    
    }
    close(CNVSUMMARY);
}

sub identifyFlexibleBlocksAndStore{

    #### AoA lines are modifed VCF format e.g.:
    # 0  1  2              3   4       5     6      7     8      9   10      11    12     13    14
    # Chr Coord Snp        GT  LRR     IA    BAF    IB    GC     GT  LRR     IA    BAF    IB    GC
    # MT 72 2010-08-MT-841 1/1 -0.4763 1.478 0.0981 0.535 0.3296 1/1 -0.2674 1.717 0.0894 0.604 0.3297
    for $i (0 .. ($#AoA - 1)){ # DON'T DO THE LAST LINE BECAUSE NO WINDOW EXISTS

        $thisChr   = $AoA[$i][0];
        $nextChr   = $AoA[($i+1)][0];
        
        # At the end of the file/AoA reduce the window size to only compare up to the last line
        my $thisWindowSize = $windowSize;
        if($i > ($#AoA - ($windowSize - 1))){
            $thisWindowSize = (($#AoA - $i) + 1);
        }
        for($j = 0; $j < ($numberOfSamples - 1); $j++){
            $thisIndexNumber = 3 + (($j + 1) * 6); # 9 15 21 i.e. GT fields
            my $countInWindow = 0;
            my @mismatchPositions;
        
            # count number of mismatches between control and iPS inside the window
            for($k = 1; $k < $thisWindowSize; $k++){
                $mismatchPositions[$k] = 0;
                if((($AoA[$i][3] ne $AoA[$i][$thisIndexNumber]) && ($AoA[$i+$k][3] ne $AoA[$i+$k][$thisIndexNumber])) && ($thisChr eq $nextChr)){
                    if($countInWindow == 0){
                        $countInWindow = 2; # +1 to include the first position
                    }
                    else{
                        $countInWindow++;
                    }
                    $mismatchPositions[$k] = 1; # keep track of where the mismatches were down the window
                }
            }
            #if enough mismatches are found inside this sample's window
            #check if the window already contains a label
            #if it does, label mismatches using that label
            #otherwise label mismatches using last label used for that sample +1
            if($countInWindow >= $minRequiredInWindow){
                my $nonZeroBlockLabel = 0;
        
                # check for any previously set labels in the window
                for($k = 0; $k < $thisWindowSize; $k++){
                    my $AoAWindowIter = $i + $k;
                    if($blocksAoA[$AoAWindowIter][$j] > 0){
                        $nonZeroBlockLabel = $blocksAoA[$AoAWindowIter][$j];
                    }
                }
        
                #if a label previously set was found, use that, otherwise use that last label used for this sample +1
                my $thisBlockLabel = 0;
                if($nonZeroBlockLabel == 0){
                    $thisBlockLabel = $lastBlockLabel[$j] + 1;
                    $lastBlockLabel[$j] = $thisBlockLabel;
                }
                else{
                    $thisBlockLabel = $nonZeroBlockLabel;
                }
        
                #set for the start of the window
                $blocksAoA[$i][$j] = $thisBlockLabel;
        
                #set label for any mismatches that were found in the window
                for($k = 1; $k < $thisWindowSize; $k++){
                    my $thisPos = $i + $k;
                    if($mismatchPositions[$k] == 1){
                        $blocksAoA[$thisPos][$j] = $thisBlockLabel;
                    }
                }
            }
        }

    }
    
}

##

sub printField{
    my $fieldValue = $_[0];
    my $totalSize = $_[1];
    my $valueLength = length($fieldValue);
    print "$fieldValue";
    for(my $j = 0; $j < ($totalSize - $valueLength); $j++){
        print " ";
    }
}

sub printAoA{

    for $i (0 .. $#AoA){
        my $aRef = $AoA[$i];
        my $n = @$aRef - 1;
        for $j (0 .. $n){
            print "$AoA[$i][$j]\t";
        }
        print "\n";
    }

    # OR:

    #print "0:$AoA[$i][0] ";
    #print "$AoA[$i][1] ";
    #print "$AoA[$i][2] ";
    #print "$AoA[$i][3] ";
    #print "$AoA[$i][4] ";
    #print "$AoA[$i][5] ";
    #print "$AoA[$i][6] ";
    #print "$AoA[$i][7] ";
    #print "$AoA[$i][8] ";
    #for($j = 0; $j < ($numberOfSamples - 1); $j++){
    #    $thisIndexNumber = 3 + (($j + 1) * 6); # 9 15 21
    #    print "$thisIndexNumber:$AoA[$i][$thisIndexNumber] ";
    #    print "$AoA[$i][$thisIndexNumber+1] ";
    #    print "$AoA[$i][$thisIndexNumber+2] ";
    #    print "$AoA[$i][$thisIndexNumber+3] ";
    #    print "$AoA[$i][$thisIndexNumber+4] ";
    #    print "$AoA[$i][$thisIndexNumber+5] ";
    #
    #}
    #print "\n";
}

##

sub printCNVHeader{

    print "\nDIFFERENCES TO CONTROL GENOTYPE\n";
    print "===============================\n\n";

    print "Chr  Coord       SNP            ";
    print "Control Sample:          | iPS Samples:\n";
    print "                                ";
    for($i = 0; $i < $numberOfSamples; $i++){
        if($i == 0){
            &printField("$sampleNames[$i]", "25");
            print "| ";
        }
        else{
            &printField($sampleNames[$i], "27");
        }
    }
    print "\n";
    print "                                ";
    for($i = 1; $i <= $numberOfSamples; $i++){
        if($i == 1){
            print "GT       IA       IB     | ";
        }
        else{
            print "GT       IA       IB       ";
        }
    }
    print "\n";
    for($i = 1; $i <= (($numberOfSamples * 28) + 29); $i++){
        print "=";
    }
    print "\n";

}

sub readVCFAndShowAllGTDifferences{

    my $duplicateSNP;
    my $noCallFound;
    my $mismatchFound;

    while(<VCF>){
        chomp;

        #### Skip the header section
        if(/#/){
            next;
        }
        ####

        @tmpArr = split /\t/;
        $totalEntryCount++;

        #### Check for duplicate SNPs on successive lines, skip 2nd instance if found
        $duplicateSNP = 0;
        if(($tmpArr[0] eq $lastChr) && ($tmpArr[1] == $lastCoord)){
            $duplicateSNP = 1;
        }
        $lastChr = $tmpArr[0];
        $lastCoord = $tmpArr[1];
        if($duplicateSNP == 1){
            next;
        }
        ####

        #### Check for no call, skip line if present
        $noCallFound = 0;
        for($i = 1; $i <= $numberOfSamples; $i++){
            $thisIndexNumber = 8 + $i; # TAKE THE 9 15 21 (GTs)
            if($tmpArr[$thisIndexNumber] eq '.'){
                $noCallFound = 1;
            }
        }
        if($noCallFound == 1){
            next;
        }
        ####

        #### Make $_ = CHROM POS ID + the samples' GT info (i.e. don't include REF ALT QUAL FILTER INFO FORMAT)
        $_ = "$tmpArr[0]\t$tmpArr[1]\t$tmpArr[2]";
        for($i = 0; $i < $numberOfSamples; $i++){
            $thisIndexNumber = 9 + $i; # 9 10 11 12 ..
            $_ = $_ . "\t$tmpArr[$thisIndexNumber]";
        }

        #### Put tabs between GT:LRR:IA:BAF:IB:GC and the split everything on tabs
        s/\:/\t/g;
        @tmpArr = split /\t/; # E.g. tmpArr: MT 72 2010-08-MT-841 1/1 -0.4763 1.478 0.0981 0.535 0.3296 1/1 -0.2674 1.717 0.0894 0.604 0.3297 

        #### READ IN THE IAs & IBs FOR EVERYTHING INTO totalIAs and totalIBs
        for($i = 0; $i < $numberOfSamples; $i++){
            $IAIndex = (($i * 6) + 3) + 2;
            $IBIndex = (($i * 6) + 3) + 4;
            $totalIAs[$i] = $totalIAs[$i] + $tmpArr[$IAIndex];
            $totalIBs[$i] = $totalIBs[$i] + $tmpArr[$IBIndex];
        }

        #### Compare all the GTs to find mismatches eg 1/0 1/0 1/1 1/0
        $mismatchFound = 0;
        $totalMismatchesPerRowCount = 0;
        for($i = 1; $i < $numberOfSamples; $i++){
            $thisIndexNumber = 3 + ($i * 6); # TAKE THE 9 15 21 e.g. 1/1
            if($tmpArr[3] ne $tmpArr[$thisIndexNumber]){
                $totalMismatchesPerRowCount++; # count the number of mismatches per row
                $mismatchFound = 1;
            }
        }
        #### Store the sum of mismatches per row  e.g. 0:10 1:0 2:4 etc. = 10 rows have 0 mismatches, 0 rows have 1 mismatch, 4 rows have 2 mismatches.
        $totalMismatchesPerRow[$totalMismatchesPerRowCount] = ($totalMismatchesPerRow[$totalMismatchesPerRowCount] + 1); 

        #### STORE ALL LINES for later processing in the modified tmpArr format e.g. MT 72 2010-08-MT-841 1/1 -0.4763 1.478 0.0981 0.535 0.3296 1/1 -0.2674 1.717 0.0894 0.604 0.3297 
        push (@AoA, [@tmpArr]);

        #### For GT mismatches print out summary e.g. 1    152883711   exm101749      1/1      0.535    0.174  | 1/0*     0.317    0.226    1/1      0.414    0.155    1/0*     0.361    0.217
        if($mismatchFound == 1){

            # Print the basic SNP info e.g. 1    152883711   exm101749
            &printField($tmpArr[0], "5");
            &printField($tmpArr[1], "12");
            &printField($tmpArr[2], "15");

            #Get the GT, IA and IB fields for each sample:
            #0  1  2              3   4       5     6      7     8      9   10      11    12     13    14
            #Chr Coord Snp        GT  LRR     IA    BAF    IB    GC     GT  LRR     IA    BAF    IB    GC
            #MT 72 2010-08-MT-841 1/1 -0.4763 1.478 0.0981 0.535 0.3296 1/1 -0.2674 1.717 0.0894 0.604 0.3297
            for($i = 0; $i < ($numberOfSamples); $i++){                
                for($j = 0; $j < 3; $j++){ # TAKE THE 3 5 7 
                    # Pattern is:
                    # 3,  5,  7  = GT IA IB
                    # 9,  11, 13 = GT IA IB
                    # 15, 17, 19 = GT IA IB
                    # 21, 23, 25 = GT IA IB
                    # etc.
                    if($j == 0){    # For GT
                        $thisIndexNumber = 3 + ($i * 6) + (0);
                    }
                    elsif($j == 1){ # For IA
                        $thisIndexNumber = 3 + ($i * 6) + (2);
                    }
                    elsif($j == 2){ # For IB
                        $thisIndexNumber = 3 + ($i * 6) + (4);
                    }
                    
                    # Count the total number of mismatches per sample e.g. 24 2 6 etc.
                    # and IAs & IBs per sample at mismatches e.g. 12.2 3.2 12.8 etc.
                    if(($i != 0) && ($j == 0)){ # At the iPS GTs
                        if($tmpArr[3] ne $tmpArr[$thisIndexNumber]){ # mismatch found in sample GT compared to the control GT e.g. 1/0 <-> 1/1
                            $printValue = "$tmpArr[$thisIndexNumber]*";
                            $allGTMismatchesPerSample[$i] = ($allGTMismatchesPerSample[$i] + 1); # total number of changes for each sample  e.g. 24 2 6 etc.
                            $IAsAtAllGTMismatches[$i] = $IAsAtAllGTMismatches[$i] + $tmpArr[($thisIndexNumber + 2)]; # total intensity at mismatches per sample 12.2 3.2 12.8 etc.
                            $IBsAtAllGTMismatches[$i] = $IBsAtAllGTMismatches[$i] + $tmpArr[($thisIndexNumber + 4)];
                        }
                        else{ # NO MISMATCH FOUND
                            $printValue = $tmpArr[$thisIndexNumber];
                        }
                    }
                    else{ # NOT GT OR 1st GT
                        $printValue = $tmpArr[$thisIndexNumber];
                    }
                    if(($i == 0) && ($j == 2)){ # MARK THE LAST FIELD OF THE CONTROL
                        &printField("$printValue", "6");
                        print " | ";
                    }
                    else{ # PRINT EVERYTHING THAT'S NOT THE LAST FIELD OF THE CONTROL
                        &printField("$printValue", "9");
                    }
                    
                }
            }
            print "\n";
        }
    }
    close(VCF);

    print "\nGT = Genotype\n";
    print "IA = Intensity of A allele\n";
    print "IB = Intensity of B allele\n";
    print "*  = iPS GT different to control's GT\n";
    print "\n";
}

sub printSummaryOfAllGTDifferences{

    print "SUMMARY OF ALL GT DIFFERENCES\n";
    print "=============================\n\n";
    
    #ALL CHANGES SUMMARY
    print "                                ";
    for($i = 0; $i <= $#sampleNames; $i++){
        if($i == 0){
            &printField("$sampleNames[$i]", "24");
            print " | ";
        }
        else{
            &printField("$sampleNames[$i]", "27");
        }
    }
    print "\n";
    
    print "                                ";
    for($i = 1; $i <= ($numberOfSamples * 27); $i++){
        print "=";
    }
    print "\n";
    
    print "Number of GT changes:           ";
    print "-                        | "; # SELF COMPARISON ALWAYS 0 SO DON'T CALCULATE
    for($i = 1; $i < $numberOfSamples; $i++){
        &printField("$allGTMismatchesPerSample[$i]", 27);
    }
    print "\n";
    
    print "Average IA at GT changes:       ";
    print "-                        | "; # SELF COMPARISON ALWAYS 0 SO DON'T CALCULATE
    for($i = 1; $i < $numberOfSamples; $i++){
        my $thisIAAtAllChanges = $IAsAtAllGTMismatches[$i] / $allGTMismatchesPerSample[$i];
        $rounded = sprintf("%.3f", $thisIAAtAllChanges);
        &printField("$rounded", "27");
    }
    print "\n";
    
    print "Average IA for $totalEntryCount SNPs:     ";
    for($i = 0; $i < $numberOfSamples; $i++){
        $average = $totalIAs[$i] / $totalEntryCount;
        $rounded = sprintf("%.3f", $average);
        if($i == 0){
            &printField("$rounded", "24");
            print " | ";
        }
        else{
            &printField("$rounded", "27");
        }
    }
    print "\n";
    
    print "Average IB at GT changes:       ";
    print "-                        | "; # SELF COMPARISON ALWAYS 0 SO DON'T CALCULATE
    for($i = 1; $i < $numberOfSamples; $i++){
        my $thisIBAtAllChanges = $IBsAtAllGTMismatches[$i] / $allGTMismatchesPerSample[$i];
        $rounded = sprintf("%.3f", $thisIBAtAllChanges);
        &printField("$rounded", "27");
    }
    print "\n";
    
    print "Average IB for $totalEntryCount SNPs:     ";
    for($i = 0; $i < $numberOfSamples; $i++){
        $average = $totalIBs[$i] / $totalEntryCount;
        $rounded = sprintf("%.3f", $average);
        #&printField("$rounded", "27");
        if($i == 0){
            &printField("$rounded", "24");
            print " | ";
        }
        else{
            &printField("$rounded", "27");
        }
    }
    print "\n\n";
    
    #######################
    
    #MARKERS
    
    print "Total SNPs:                                 "; 
    &printField("$totalEntryCount SNPs\n", 10);
    for($i = 0; $i < $numberOfSamples; $i++){
        print "$i / " . ($numberOfSamples - 1) . " samples are different to the control: ";
        #print " SNPs\n";
        &printField("$totalMismatchesPerRow[$i] SNPs", 10);
        print "\n";
    }
    print "\n";

}

###################################################################################

##NOTES

#grep -v '#' EpiFS18.vcf | cut -f 1,2,10-13 | sed 's/:/\t/g' | cut -f 1-3,9,15,21 | grep -v '\.' | awk '($3 != $4 || $3 != $5 || $3 != $6)'

