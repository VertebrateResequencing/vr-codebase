#!/usr/bin/env perl

use warnings;
use strict;

use Getopt::Std;

###################################################
# THIS IS STILL IN DEVELOPMENT - VERSION 0.1
#
# PROGRAM FOR IDENTIFYING LOH REGIONS USING VCFs
# GENERATED FROM GENOTYPING DATA FOR HipSci PROJECT
# THIS IS NOT INTENDED AS A GENERIC CALLER
# THIS PROGRAM IS USED TO FIND LOH REGIONS 
# COMPARED TO FIBROBLAST/CONTROL SAMPLE THAT 
# A HipSci iPS SAMPLE HAS BEEN DERIVED FROM
# THROUGH GENETIC REPROGRAMMING
# 
# THE VCF MUST HAVE THE CONTROL HAS THE FIRST 
# SAMPLE WITH iPS SAMPLES AFTER IT AS YOU GO 
# FROM LEFT TO RIGHT IN THE VCF I.E. 
# CONTROL iPS1 iPS2 iPS3
###################################################

# declare the perl command line flags/options we want to allow
my %options=();
getopts("dsbc:v:l:", \%options);

my $concisePrint = 0;
my $bedPrint = 0;
my $debugPrint = 0;
my $controlSampleName;
my $VCFFile;
my $logfile;

if(defined $options{v}){
    $VCFFile = $options{v};
    if (-e $VCFFile) {
        #print "VCF found OK\n";
    } 
    else{
        print "$VCFFile VCF doesn't exist\n";
        print "Usage: loh.pl -v yourvcf.vcf -c controlSampleName [-s for short output -b BED output -l logfile]\n";
        exit;
    }
}
else{
    print "VCF file required\n";
    print "Usage: loh.pl -v yourvcf.vcf -c controlSampleName [-s for short output -b BED output -l logfile]\n";
    exit;
}
if(defined $options{c}){
    $controlSampleName = $options{c};
}
else{
    print "Control sample name required\n";
    print "Usage: loh.pl -v yourvcf.vcf -c controlSampleName [-s for short output -b BED output -l logfile]\n";
    exit;
}
if(defined $options{l}){
    $logfile = $options{l};
}
else{
    print "logfile required\n";
    print "Usage: loh.pl -v yourvcf.vcf -c controlSampleName -l logfile [-s for short output -b BED output]\n";
    exit;
}
if(defined $options{s}){
    $concisePrint = 1;
}
if(defined $options{b}){
    $bedPrint = 1;
}
if (defined $options{d}){
    $debugPrint = 1;
}
if ($ARGV[0]){
    print "Unexpected arguments:\n";
    foreach (@ARGV){
        print "$_\n";
    }
    print "Usage: loh.pl -v yourvcf.vcf -c controlSampleName [-s for short output -b BED output -l logfile]\n";
    exit;
}

if($concisePrint == 0){
    print "Control sample name = $controlSampleName\n";
}

my(
    $i, $j, $k, $l, $m, $n, 
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
    $nextCoord
);

if($concisePrint == 0){
    print "VCF file = $VCFFile\n";
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

my $controlPosition = 0;
$i = 0;
#print "numberOfSamples = $numberOfSamples\n";
foreach(@sampleNames){
    $i++;
    if($_ eq $controlSampleName){
        $controlPosition = $i;
    }
}
if($controlPosition == 0){
    print "Control sample $controlSampleName not found in $VCFFile\n";
    print "Usage: loh.pl -v yourvcf.vcf -c controlSampleName [-s for short output]\n";
    exit;
}
#$controlPosition = 2; # TESTING
# PRINT THE HEADER OF THE CNV OUTPUT
if($concisePrint == 0){
    printCNVHeader();
}

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
if($concisePrint == 0){
    printSummaryOfAllGTDifferences();
}
#printAoA();

#### AoA lines are modifed VCF format e.g.:
# 0  1  2              3   4       5     6      7     8      9   10      11    12     13    14
# Chr Coord Snp        GT  LRR     IA    BAF    IB    GC     GT  LRR     IA    BAF    IB    GC
# MT 72 2010-08-MT-841 1/1 -0.4763 1.478 0.0981 0.535 0.3296 1/1 -0.2674 1.717 0.0894 0.604 0.3297

#### blocksAoA
#0 0 0
#1 0 0
#1 0 2
#0 0 2

#@allTrimmedBlocks structure:
#                      if 0 length, no CNVs    
#[sample#][0=SNP overall index pos in AoA;1=mismatches][CNV#?][n for SNP index pos; 0/1 for mismatches]

my $windowIter;
my $percentageFlexIter;
my @storedResults;
my @allTrimmedBlocks;
my @blocksAoA;
my @listOfAllTrimmedBlocks;

#my $windowIterStart = 7;
#my $windowIterEnd = 8;
#my $percentageFlexIterStart = 70;
#my $percentageFlexIterEnd = 80;
#my $windowIterStart = 7;
#my $windowIterEnd = 8;
#my $percentageFlexIterStart = 70;
#my $percentageFlexIterEnd = 90;
my $windowIterStart = 7;
my $windowIterEnd = 8;
my $percentageFlexIterStart = 70;
my $percentageFlexIterEnd = 80;
#my $windowIterStart = 2;
#my $windowIterEnd = 15;
#my $percentageFlexIterStart = 10;
#my $percentageFlexIterEnd = 90;

# LOOP THROUGH THE WINDOW SIZE AND PERCENTAGE FLEXIBILITY PARAMETERS AND PRODUCE @listOfAllTrimmedBlocks

my $paramIter = 1;
if($concisePrint == 0){
    print "Performed searches using the following parameter combinations:\n";
}
for($windowIter = $windowIterStart; $windowIter <= $windowIterEnd; $windowIter++){
    for($percentageFlexIter = $percentageFlexIterStart; $percentageFlexIter < $percentageFlexIterEnd; $percentageFlexIter = ($percentageFlexIter + 10)){
        $windowSize            = $windowIter;
        $percentageFlexibility = $percentageFlexIter;
        if($concisePrint == 0){
            print "[$paramIter]: $windowSize SNPs in CNV window & $percentageFlexibility" . "% flexibility\n";
        }
        $paramIter++;

        # STORE EMPTY RESULT AND JUMP TO NEXT ITER IF minRequiredInWindow IS TOO SMALL
        $minRequiredInWindow = 0;
        getMinRequiredInWindow();
        if($minRequiredInWindow < 2){
            if($concisePrint == 0){
                print "percentageFlexibility is too high, need a minium of 2 SNPs in a window.\n";
            }
            my @skippedResult;
            for ($i = 0; $i < ($numberOfSamples - 1); $i++){
                $skippedResult[$i] = -1;
            }
            push (@storedResults, [@skippedResult]);

            $#allTrimmedBlocks = -1;
            push(@listOfAllTrimmedBlocks, [@allTrimmedBlocks]);

            next;
        }
        
        #AoA
        $#blocksAoA = -1;
        # Create an empty @blocksAoA[][]
        for $i (0 .. $#AoA){
            for($j = 0; $j < ($numberOfSamples - 1); $j++){
                $inBlockArr[$j] = 0;
            }
            push (@blocksAoA, [@inBlockArr]);
        }
        identifyFlexibleBlocksAndStore(); 
        if($debugPrint == 1){
            printFlexibleBlocksNumbers();
        }
        #blocksAoa now obtained
        $#allTrimmedBlocks = -1;
        trimBlocks();
        #allTrimmedBlocks now obtained
        if($debugPrint == 1){
            printTrimmedBlocksDebug();
        }
        getCNVSizes();
        #storedResults now obtained
        push(@listOfAllTrimmedBlocks, [@allTrimmedBlocks]); # BUT WHAT ABOUT THE EMPTY storedResult ABOVE?
        #listOfAllTrimmedBlocks now obtained
    }
}

#### AoA lines are modifed VCF format e.g.:
# 0  1  2              3   4       5     6      7     8      9   10      11    12     13    14
# Chr Coord Snp        GT  LRR     IA    BAF    IB    GC     GT  LRR     IA    BAF    IB    GC
# MT 72 2010-08-MT-841 1/1 -0.4763 1.478 0.0981 0.535 0.3296 1/1 -0.2674 1.717 0.0894 0.604 0.3297

#### blocksAoA
#0 0 0
#1 0 0
#1 0 2
#0 0 2

#storedResults
#[0][0] 0 [0][1] 6 [0][2] 4 
#[1][0] 0 [1][1] 9 [1][2] 4 
#[2][0] 0 [2][1] 6 [2][2] 4 
#[3][0] 0 [3][1] 9 [3][2] 4 

#@allTrimmedBlocks structure:
#                      if 0 length, no CNVs    
#[sample#][0=SNP overall index pos in AoA;1=mismatches][CNV#?][n for SNP index pos; 0/1 for mismatches]

#@listOfAllTrimmedBlocks structure:
#                         if 0 length, no CNVs    
#[param combination#][sample#][0=SNP overall index pos in AoA;1=mismatches][CNV#?][n for SNP index pos; 0/1 for mismatches]

if($debugPrint == 1){
    print "The listOfAllTrimmedBlocks is $#listOfAllTrimmedBlocks long\n";
    print "The stored results has a depth of $#storedResults\n";
}

#SCAN storedResults TO FIND THE MAX NUMBER OF SNPs FOR EACH SAMPLE AND STORE IN $maxNumberOfSNPs
my @maxNumberOfSNPs;
my $aRef = $storedResults[0];
my $p = @$aRef - 1;
for $j (0 .. $p){
    $maxNumberOfSNPs[$j][0] = 0;
    $maxNumberOfSNPs[$j][1] = 0;
}
for $i (0 .. $#storedResults){
    my $bRef = $storedResults[$i];
    my $q = @$bRef - 1;
    for $j (0 .. $q){    
        #print "[$i][$j] $storedResults[$i][$j] ";
        if($storedResults[$i][$j] > $maxNumberOfSNPs[$j][0]){
            $maxNumberOfSNPs[$j][0] = $storedResults[$i][$j];
            $maxNumberOfSNPs[$j][1] = $i;
        }
    }
    #print "\n";
}

#@listOfAllTrimmedBlocks structure:
#                         if 0 length, no CNVs    
#[param combination#][sample#][0=SNP overall index pos in AoA;1=mismatches][CNV#?][n for SNP index pos; 0/1 for mismatches]

if($concisePrint == 0){
    print "\n";
    print "LOH REGIONS DETECTED\n";
    print "====================\n\n";
}

my $bedOutputCount = 0;
#if($logfile){
open(LOGFILE, ">$logfile") || die "cannot open $logfile ($!)";
#}

my $cRef = $storedResults[0];
my $r = @$cRef - 1;
for $j (0 .. $r){

    my $sampleNameIndex = $j;
    if(($j+1) >= $controlPosition){
        $sampleNameIndex = $j + 1;
    }
    if($bedPrint == 0){
        print $sampleNames[($sampleNameIndex)] . ":\n";
    #}
    #if($concisePrint == 0){
        for(my $x = 0; $x < length($sampleNames[($j+1)]); $x++){
            print "=";
        }
        print "\n";
    }
    my $thisCombinationPosition = $maxNumberOfSNPs[$j][1];
    my $dRef = $listOfAllTrimmedBlocks[$thisCombinationPosition][$j];
    my $s = @$dRef - 1;
    if($s == -1){
        print "0 CNVs containing 0 SNPs\n\n";
        if($logfile){
            print LOGFILE "$sampleNames[($sampleNameIndex)]" . ": 0 CNVs\n";
        }
    }
    if($s == 1){
        my $eRef = $listOfAllTrimmedBlocks[$thisCombinationPosition][$j][0]; # accessing combo, sample#, 0=SNPs
        my $t = @$eRef - 1;
        if($logfile){
            print LOGFILE "$sampleNames[($sampleNameIndex)]" . ": " . ($t+1) . " CNV";
            if($t > 0){
                print LOGFILE "s";
            }
            print LOGFILE "\n";
        }
        if($concisePrint == 0){
            print "" . ($t+1) . " CNV";
            if($t > 0){
                print "s";
            }
            print " containing ";
            print "$maxNumberOfSNPs[$j][0] SNPs (";

            my $iterCount = 0;
            for($windowIter = $windowIterStart; $windowIter <= $windowIterEnd; $windowIter++){
                for($percentageFlexIter = $percentageFlexIterStart; $percentageFlexIter < $percentageFlexIterEnd; $percentageFlexIter = ($percentageFlexIter + 10)){
                    if($iterCount == $maxNumberOfSNPs[$j][1]){
                        print "window size = $windowIter flexibility = $percentageFlexIter" . "%";
                    }
                    $iterCount++;
                }
            }
            print ")";
            print "\n";
        }
        if($t >= 0){
            for $k (0 .. $t){
                my $fRef = $listOfAllTrimmedBlocks[$thisCombinationPosition][$j][0][$k]; # CNVs
                my $u = @$fRef - 1;
                my $startOfLine = '';
                for $l (0 .. $u){
                    my $localIndex = $listOfAllTrimmedBlocks[$thisCombinationPosition][$j][0][$k][$l]; # actual SNP positions
                    my $localMismatch = $listOfAllTrimmedBlocks[$thisCombinationPosition][$j][1][$k][$l]; # mismatch?
                    if($l == 0){
                        my $endIndex = $listOfAllTrimmedBlocks[$thisCombinationPosition][$j][0][$k][$u]; # actual SNP positions
                        if($bedPrint == 0){
                            $startOfLine = "CNV #" . ($k+1) . ": Summary: " . ($u+1) . "SNPs Chr" . "$AoA[$localIndex][0]:";
                            print "$startOfLine";
                            print "$AoA[$localIndex][1]-";
                            print "$AoA[$endIndex][1]\n";
                        }
                        if($bedPrint == 1){
                            print "$AoA[$localIndex][0]\t";
                            print "$AoA[$localIndex][1]\t";
                            print "$AoA[$endIndex][1]\t";
                            print $sampleNames[($sampleNameIndex)] . "\t";
                            print "" . ($u+1) . " (SNPs)\n";
                            $bedOutputCount++;
                        }
                        if($bedPrint == 0){
                            $startOfLine = "CNV #" . ($k+1) . ": Chr:$AoA[$localIndex][0] ";
                            print "$startOfLine";
                            &printField("Coord:$AoA[$localIndex][1]", "16");
                            &printField("SNP:$AoA[$localIndex][2]", "19");
                            print "Mismatch:$localMismatch\n";
                        }
                    }
                    else{
                        if($bedPrint == 0){
                            for(my $buff = 0; $buff < length($startOfLine); $buff++){
                                print " ";
                            }
                            &printField("Coord:$AoA[$localIndex][1]", "16");
                            &printField("SNP:$AoA[$localIndex][2]", "19");
                            print "Mismatch:$localMismatch\n";
                        }
                    }
                }
            }
        }
    }
    #print "\n";
}

print LOGFILE "$bedOutputCount lines of BED output produced\n";
close(LOGFILE);

#######################

sub getCNVSizes{

    my @storedResult;
    for $i (0 .. $#allTrimmedBlocks){
        $storedResult[$i] = 0;
    }

    for $i (0 .. $#allTrimmedBlocks){
        my $aRef = $allTrimmedBlocks[$i];
        my $p = @$aRef - 1;
        if($p == -1){
            $storedResult[$i] = 0;
        }
        if($p == 1){
            $j = 0; # because it's a fixed length i.e. always 2
            my $bRef = $allTrimmedBlocks[$i][$j];
            my $q = @$bRef - 1;
            if($q >= 0){
                for $k (0 .. $q){
                    my $cRef = $allTrimmedBlocks[$i][$j][$k];
                    my $r = @$cRef - 1;
                    if($r >= 0){
                        $storedResult[$i] += ($r+1);
                    }
                }
            }
        }
    }

#Eg
#[0][0] 0 [0][1] 0 [0][2] 0 
#[1][0] 0 [1][1] 5 [1][2] 0 
#[2][0] 0 [2][1] 0 [2][2] 0 
#[3][0] 0 [3][1] 5 [3][2] 0 
#[4][0] 0 [4][1] 0 [4][2] 0 
#[5][0] 0 [5][1] 5 [5][2] 0 
#[6][0] 0 [6][1] 0 [6][2] 0 
#[7][0] 0 [7][1] 5 [7][2] 0 

    if($debugPrint == 1){
        for $i (0 .. $#allTrimmedBlocks){
            print "storedResult[$i] = $storedResult[$i]\n";
        }
    }
    push (@storedResults, [@storedResult]);
}

sub printTrimmedBlocksDebug{

    print "Printing allTrimmedBlocks debug\n";

    for $i (0 .. $#allTrimmedBlocks){
        my $aRef = $allTrimmedBlocks[$i];
        my $p = @$aRef - 1;
        print "\ni:$i ";
        print "\np:$p ";
        if($p == 1){
            for $j (0 .. $p){
                my $bRef = $allTrimmedBlocks[$i][$j];
                my $q = @$bRef - 1;
                print "j:$j ";
                print "q:$q ";
                if($q >= 0){
                    for $k (0 .. $q){
                        my $cRef = $allTrimmedBlocks[$i][$j][$k];
                        my $r = @$cRef - 1;
                        print "k:$k ";
                        print "r:$r ";
                        if($r >= 0){
                            for $l (0 .. $r){
                                print "[$i][$j][$k][$l]: $allTrimmedBlocks[$i][$j][$k][$l]\t";
                            }
                            print "\n";
                        }
                        else{
                            print "Shouldn't get here 1\n";
                            exit;
                        }
                    }
                }
                else{
                    print "Shouldn't get here 2\n";
                    exit;
                }
            }
        }
        else{
            print "\n";
        }
    }

}

sub trimBlocks{

    #my $blockFound;
    my @lastLocalBlockLabel;

    my $start;

    for($j = 0; $j < ($numberOfSamples - 1); $j++){
        $lastLocalBlockLabel[$j] = 0; 
    }

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

}

sub identifyFlexibleBlocksAndStore{

    #### AoA lines are modifed VCF format e.g.:
    # 0  1  2              3   4       5     6      7     8      9   10      11    12     13    14
    # Chr Coord Snp        GT  LRR     IA    BAF    IB    GC     GT  LRR     IA    BAF    IB    GC
    # MT 72 2010-08-MT-841 1/1 -0.4763 1.478 0.0981 0.535 0.3296 1/1 -0.2674 1.717 0.0894 0.604 0.3297

    my @lastBlockLabel;

    # Create an empty lastBlockLabel
    for($j = 0; $j < ($numberOfSamples - 1); $j++){
        $lastBlockLabel[$j] = 0;
    }

    for $i (0 .. ($#AoA - 1)){ # DON'T DO THE LAST LINE BECAUSE NO WINDOW EXISTS
        $thisChr   = $AoA[$i][0];
        $nextChr   = $AoA[($i+1)][0];
        
        # At the end of the file/AoA reduce the window size to only compare up to the last line
        my $thisWindowSize = $windowSize;
        if($i > ($#AoA - ($windowSize - 1))){
            $thisWindowSize = (($#AoA - $i) + 1);
        }
        for($j = 0; $j < ($numberOfSamples); $j++){
            if(($j+1) == $controlPosition){
                next;
            }
            my $blockStoreIndex = $j;
            if(($j+1) > $controlPosition){
                $blockStoreIndex = ($j - 1);
            }
            my $controlIndexNumber = 3 + (($controlPosition - 1) * 6);
            #$thisIndexNumber = 3 + (($j + 1) * 6); # 9 15 21 i.e. GT fields
            $thisIndexNumber = 3 + (($j) * 6); # 9 15 21 i.e. GT fields
            my $countInWindow = 0;
            my @mismatchPositions;

            # count number of mismatches between control and iPS inside the window
            for($k = 1; $k < $thisWindowSize; $k++){
                $mismatchPositions[$k] = 0;
                #if((($AoA[$i][3] ne $AoA[$i][$thisIndexNumber]) && ($AoA[$i+$k][3] ne $AoA[$i+$k][$thisIndexNumber])) && ($thisChr eq $nextChr)){
                if((($AoA[$i][$controlIndexNumber] ne $AoA[$i][$thisIndexNumber]) && ($AoA[$i+$k][$controlIndexNumber] ne $AoA[$i+$k][$thisIndexNumber])) && ($thisChr eq $nextChr)){
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
                    #if($blocksAoA[$AoAWindowIter][$j] > 0){
                    if($blocksAoA[$AoAWindowIter][$blockStoreIndex] > 0){
                        $nonZeroBlockLabel = $blocksAoA[$AoAWindowIter][$blockStoreIndex];
                        #$nonZeroBlockLabel = $blocksAoA[$AoAWindowIter][$j];
                    }
                }
        
                #if a label previously set was found, use that, otherwise use that last label used for this sample +1
                my $thisBlockLabel = 0;
                if($nonZeroBlockLabel == 0){
                    #$thisBlockLabel = $lastBlockLabel[$j] + 1;
                    $thisBlockLabel = $lastBlockLabel[$blockStoreIndex] + 1;
                    #$lastBlockLabel[$j] = $thisBlockLabel;
                    $lastBlockLabel[$blockStoreIndex] = $thisBlockLabel;
                }
                else{
                    $thisBlockLabel = $nonZeroBlockLabel;
                }
        
                #set for the start of the window
                #$blocksAoA[$i][$j] = $thisBlockLabel;
                $blocksAoA[$i][$blockStoreIndex] = $thisBlockLabel;
        
                #set label for any mismatches that were found in the window
                for($k = 1; $k < $thisWindowSize; $k++){
                    my $thisPos = $i + $k;
                    if($mismatchPositions[$k] == 1){
                        #$blocksAoA[$thisPos][$j] = $thisBlockLabel;
                        $blocksAoA[$thisPos][$blockStoreIndex] = $thisBlockLabel;
                    }
                }
            }
        }

    }
    
}

##

sub getMinRequiredInWindow{

    my (
        $flexibilitySubtractionFromWindow, 
        $flexibilitySubtractionFromWindowRoundedDown);

    $flexibilitySubtractionFromWindow = 0;
    $flexibilitySubtractionFromWindow            = ($windowSize / 100) * $percentageFlexibility;
    $flexibilitySubtractionFromWindowRoundedDown = int $flexibilitySubtractionFromWindow;
    $minRequiredInWindow                         = $windowSize - $flexibilitySubtractionFromWindowRoundedDown;

}

##

sub printFlexibleBlocksNumbers{

    my $blockFound;

    print "Printing flexible block numbers\n";

    for $i (0 .. ($#blocksAoA)){
        $blockFound = 0;
        for($j = 0; $j < ($numberOfSamples - 1); $j++){
            if($blocksAoA[$i][$j] > 0){
                $blockFound++;
            }
        }
        if($blockFound > 0){
            #print "i=$i block#=$blocksAoA[$i][0]; ";
            print "i=$i ";
            print "$AoA[$i][0] $AoA[$i][1] $AoA[$i][2]; ";
            for($j = 0; $j < ($numberOfSamples - 1); $j++){
                print "$blocksAoA[$i][$j]; ";
            }
            print "\n";
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

    #print "Control Sample:          | iPS Samples:\n";
    #print "                               ";
    my $controlIndex = ($controlPosition - 1);

    for($i = 0; $i < $numberOfSamples; $i++){
        if($i == $controlIndex){
            print " | Control Sample:       | ";
        }
        else{
            print "iPS Sample:                ";
        }
    }
    print "\n";
    print "                                ";

    for($i = 0; $i < $numberOfSamples; $i++){
        if($i == $controlIndex){
            print " | ";
            &printField("$sampleNames[$i]", "21");
            print " | ";
        }
        else{
            &printField($sampleNames[$i], "27");
        }
    }

    print "\n";
    print "                                ";
    for($i = 1; $i <= $numberOfSamples; $i++){
        #if($i == 1){
        #    print "GT       IA       IB     | ";
        if($i == ($controlIndex + 1)){
            print " | GT      IA       IB   | ";
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
            if($tmpArr[$thisIndexNumber] eq './.:.:.:.:.:.'){
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

        #### Put tabs between GT:LRR:IA:BAF:IB:GC and then split everything on tabs
        s/\:/\t/g;
        @tmpArr = split /\t/; # E.g. tmpArr: MT 72 2010-08-MT-841 1/1 -0.4763 1.478 0.0981 0.535 0.3296 1/1 -0.2674 1.717 0.0894 0.604 0.3297 

        #### READ IN THE IAs & IBs FOR EVERYTHING INTO totalIAs and totalIBs
        for($i = 0; $i < $numberOfSamples; $i++){
            $IAIndex = (($i * 6) + 3) + 3;
            $IBIndex = (($i * 6) + 3) + 4;

            if(($tmpArr[$IAIndex] eq '.') || ($tmpArr[$IBIndex] eq '.')){
                #next;
            }
            else{
                $totalIAs[$i] = $totalIAs[$i] + $tmpArr[$IAIndex];
                $totalIBs[$i] = $totalIBs[$i] + $tmpArr[$IBIndex];
            }
            # used to be    GT:LRR:IA:BAF:IB:GC
#                            0  1   2   3  4  5
#                            0  1   3   2  4  5
            # format is now GT:LRR:BAF:IA:IB:GC
        }

        #### Compare all the GTs to find mismatches eg 1/0 1/0 1/1 1/0
        $mismatchFound = 0;
        $totalMismatchesPerRowCount = 0;
        for($i = 0; $i < $numberOfSamples; $i++){
        #for($i = 1; $i < $numberOfSamples; $i++){
            if($i+1 == $controlPosition){
                next;
            }
            $thisIndexNumber = 3 + ($i * 6); # TAKE THE 9 15 21 e.g. 1/1
            my $controlGTIndexNumber = 3 + (($controlPosition - 1) * 6); # TAKE THE 9 15 21 e.g. 1/1            
            if($tmpArr[$controlGTIndexNumber] ne $tmpArr[$thisIndexNumber]){
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
            if($concisePrint == 0){
                &printField($tmpArr[0], "5");
                &printField($tmpArr[1], "12");
                &printField($tmpArr[2], "15");
            }
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
                        $thisIndexNumber = 3 + ($i * 6) + (3);
                    }
                    elsif($j == 2){ # For IB
                        $thisIndexNumber = 3 + ($i * 6) + (4);
                    }
                    my $controlIndex = ($controlPosition - 1);
                    my $storeIndex = $i;
                    if($i > $controlIndex){
                        $storeIndex = $i - 1;
                    }
                    # Count the total number of mismatches per sample e.g. 24 2 6 etc.
                    # and IAs & IBs per sample at mismatches e.g. 12.2 3.2 12.8 etc.
                    #if(($i != 0) && ($j == 0)){ # At the iPS GTs
                    if(($i != $controlIndex) && ($j == 0)){ # At the iPS GTs
                        if($tmpArr[(($controlIndex*6)+3)] ne $tmpArr[$thisIndexNumber]){ # mismatch found in sample GT compared to the control GT e.g. 1/0 <-> 1/1
                            $printValue = "$tmpArr[$thisIndexNumber]*";
                            $allGTMismatchesPerSample[$storeIndex] = ($allGTMismatchesPerSample[$storeIndex] + 1); # total number of changes for each sample  e.g. 24 2 6 etc.
                            if(($tmpArr[($thisIndexNumber + 3)] eq '.') || ($tmpArr[($thisIndexNumber + 4)] eq '.')){
                                # do nothing
                            }
                            else{
                                $IAsAtAllGTMismatches[$storeIndex] = $IAsAtAllGTMismatches[$storeIndex] + $tmpArr[($thisIndexNumber + 3)]; # total intensity at mismatches per sample 12.2 3.2 12.8 etc. # CHANGED FOR PETR VCF
                                $IBsAtAllGTMismatches[$storeIndex] = $IBsAtAllGTMismatches[$storeIndex] + $tmpArr[($thisIndexNumber + 4)];
                            }
                        }
                        else{ # NO MISMATCH FOUND
                            $printValue = $tmpArr[$thisIndexNumber];
                        }
                    }
                    else{ # NOT GT OR 1st GT
                        $printValue = $tmpArr[$thisIndexNumber];
                    }
                    #if(($i == 0) && ($j == 2)){ # MARK THE LAST FIELD OF THE CONTROL
                    if(($i == $controlIndex) && ($j == 0)){ # MARK THE LAST FIELD OF THE CONTROL
                        if($concisePrint == 0){
                            print " | ";
                            &printField("$printValue", "6");
                        }
                    }
                    elsif(($i == $controlIndex) && ($j == 2)){ # MARK THE LAST FIELD OF THE CONTROL
                        if($concisePrint == 0){
                            &printField("$printValue", "6");
                            print " | ";
                        }
                    }
                    else{ # PRINT EVERYTHING THAT'S NOT THE LAST FIELD OF THE CONTROL
                        if($concisePrint == 0){
                            &printField("$printValue", "9");
                        }
                    }                    
                }
            }
            if($concisePrint == 0){
                print "\n";
            }
        }
    }
    close(VCF);

    if($concisePrint == 0){
        print "\nGT = Genotype\n";
        print "IA = Intensity of A allele\n";
        print "IB = Intensity of B allele\n";
        print "*  = iPS GT different to control's GT\n";
        print "\n";
    }
}

sub printSummaryOfAllGTDifferences{

    print "SUMMARY OF ALL GT DIFFERENCES\n";
    print "=============================\n\n";

    #ALL CHANGES SUMMARY

    print "                                ";
    my $controlIndex = ($controlPosition - 1);

    for($i = 0; $i < $numberOfSamples; $i++){
        if($i == $controlIndex){
            print " | ";
            &printField("$sampleNames[$i]", "21");
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
    for($i = 0; $i < $numberOfSamples; $i++){
        if($i == $controlIndex){
            print " | -                     | "; # SELF COMPARISON ALWAYS 0 SO DON'T CALCULATE
        }
        else{
            my $thisSampleIndex = $i;
            if($i > $controlIndex){
                $thisSampleIndex = $i - 1;
            }
            &printField("$allGTMismatchesPerSample[$thisSampleIndex]", 27);
        }
    }
    print "\n";
    
    print "Average IA at GT changes:       ";
    for($i = 0; $i < $numberOfSamples; $i++){
        if($i == $controlIndex){
            print " | -                     | "; # SELF COMPARISON ALWAYS 0 SO DON'T CALCULATE
        }
        else{
            my $thisSampleIndex = $i;
            if($i > $controlIndex){
                $thisSampleIndex = $i - 1;
            }
            my $thisIAAtAllChanges;
            if(($IAsAtAllGTMismatches[$thisSampleIndex] == 0) || ($allGTMismatchesPerSample[$thisSampleIndex] == 0)){
                $thisIAAtAllChanges = 0;
            }
            else{
                $thisIAAtAllChanges = $IAsAtAllGTMismatches[$thisSampleIndex] / $allGTMismatchesPerSample[$thisSampleIndex];
            }
            $rounded = sprintf("%.3f", $thisIAAtAllChanges);
            &printField("$rounded", "27");
        }
    }
    print "\n";
    
    print "Average IA for $totalEntryCount SNPs:     ";
    for($i = 0; $i < $numberOfSamples; $i++){
        $average = $totalIAs[$i] / $totalEntryCount;
        $rounded = sprintf("%.3f", $average);

        if($i == $controlIndex){
            print " | ";
            &printField("$rounded", "21");
            print " | ";
        }
        else{
            &printField("$rounded", "27");
        }
    }
    print "\n";

    print "Average IB at GT changes:       ";
    for($i = 0; $i < $numberOfSamples; $i++){
        if($i == $controlIndex){
            print " | -                     | "; # SELF COMPARISON ALWAYS 0 SO DON'T CALCULATE
        }
        else{
            my $thisSampleIndex = $i;
            if($i > $controlIndex){
                $thisSampleIndex = $i - 1;
            }
            my $thisIBAtAllChanges;
            if(($IBsAtAllGTMismatches[$thisSampleIndex] == 0) || ($allGTMismatchesPerSample[$thisSampleIndex] == 0)){
                $thisIBAtAllChanges = 0;
            }
            else{
                $thisIBAtAllChanges = $IBsAtAllGTMismatches[$thisSampleIndex] / $allGTMismatchesPerSample[$thisSampleIndex];
            }
            $rounded = sprintf("%.3f", $thisIBAtAllChanges);
            &printField("$rounded", "27");
        }
    }
    print "\n";

    print "Average IB for $totalEntryCount SNPs:     ";
    for($i = 0; $i < $numberOfSamples; $i++){
        $average = $totalIBs[$i] / $totalEntryCount;
        $rounded = sprintf("%.3f", $average);

        if($i == $controlIndex){
            print " | ";
            &printField("$rounded", "21");
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

