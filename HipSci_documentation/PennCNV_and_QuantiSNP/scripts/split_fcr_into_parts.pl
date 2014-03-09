#!/bin/perl -w

use strict;

if($#ARGV != 1){
    print "Usage: $0 <INPUT GENOTYPING FILE> <OUTPUT DIR>\n";
    exit 0;
}

#SPLIT THE COMBINED GENOTYPING FILE INTO FILES PER SAMPLE
#THIS IS LIKELY TO BE REPLACED BY EITHER SPLITTING PRIOR TO IRODS 
#OR A STEP WITHIN VRPIPE

my(
    $inputFileIter, 
    $filesPrintStartedIter, 
    $inputGenotypingFile, 
    @tmpArr, 
    @filesPrintStarted, 
    $outputDir, 
    $firstLine, 
    $sampleFile, 
    $found
);

$inputGenotypingFile = $ARGV[0];
$outputDir = $ARGV[1];

$inputFileIter = 0;
open(INPUT_FILE, "$inputGenotypingFile") || die "cannot open $inputGenotypingFile ($!)";
$firstLine = <INPUT_FILE>;
while(<INPUT_FILE>){
    chomp;
    @tmpArr = split /\t/;
    $sampleFile = "$outputDir" . '/' . $tmpArr[1];
    $found = 0;
    for($filesPrintStartedIter = 0; $filesPrintStartedIter <= $#filesPrintStarted; $filesPrintStartedIter++){
        if($sampleFile eq $filesPrintStarted[$filesPrintStartedIter]){
            $found++;
        }
    }
    open(SAMPLE_FILE, ">>$sampleFile");
    if($found == 0){
        push(@filesPrintStarted, $sampleFile);
        print SAMPLE_FILE "$firstLine";
    }
    print SAMPLE_FILE "$_\n";
    close(SAMPLE_FILE);
    if((($inputFileIter % 1000000) == 0) && ($inputFileIter > 0)){ 
        print $inputFileIter+1 . " lines of $inputGenotypingFile processed\n"; 
    }
    $inputFileIter++;
}
close(INPUT_FILE);

print "Finished\n";
