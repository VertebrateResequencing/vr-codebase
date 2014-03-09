#!/bin/perl -w

use strict;

if($#ARGV != 1){
    print "Usage: $0 <RAWCNV DIR> <FILTERED CNV DIR>\n";
    exit 0;
}

#USAGE EXAMPLE:
#perl run_filter.pl rawcnv/ filteredcnv/

my(
    $rawCNVDir, 
    $filteredCNVDir, 

    $numberOfSNPs, 
    $filterDistance, 

    $PennCNVDir, 
    @rawCNVFileNames, 

    $sampleName, 
    $filteredFile, 
    $runStr
);

$rawCNVDir       = $ARGV[0];
$filteredCNVDir  = $ARGV[1];

$numberOfSNPs    = 10;
$filterDistance  = 130;

$PennCNVDir      = '/software/vertres/bin-external/PennCNV/';
@rawCNVFileNames = <$rawCNVDir/*.rawcnv>;

foreach(@rawCNVFileNames){
    /$rawCNVDir\/(.*)\.rawcnv/;
    $sampleName = $1;
    $filteredFile = $sampleName . '.filtered';
    $runStr = "perl $PennCNVDir/filter_cnv.pl -numsnp $numberOfSNPs -length " . "$filterDistance" . "k $_ --confidence 10 > $filteredCNVDir/$filteredFile";
    print "$runStr\n";
    system("$runStr");
}
