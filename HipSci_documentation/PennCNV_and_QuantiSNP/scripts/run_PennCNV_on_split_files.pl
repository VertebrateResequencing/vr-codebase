#!/bin/perl -w

use strict;

#USAGE EXAMPLE:
#perl bsub_call_CNV_on_split_Exome1.1_custom.pl split_files/ raw_PennCNV_results/ 

if($#ARGV != 1){
    print "Usage: $0 <GENOTYPE FILES DIR> <RESULTS DIRS>\n";
    exit 0;
}

my(
    $fileCount, 
    $printNumber, 
    @inputFileNames, 
    $resultsDir, 
    $PennCNVDir, 
    $genotypeFilesDir, 
    $hmmFile, 
    $pfbFile
);

$genotypeFilesDir = $ARGV[0]; # e.g. split_files/       
$resultsDir       = $ARGV[1]; # e.g. results_dirs/

$PennCNVDir         = '/software/vertres/bin-external/PennCNV/';
$hmmFile            = $PennCNVDir . 'lib/custom.hmm';
$pfbFile            = $PennCNVDir . 'lib/HumanExome12v1.1.hg19.pfb';

@inputFileNames = <$genotypeFilesDir/*_*>;

if(!(-e "$hmmFile")){
    print "Couldn't find hmm file $hmmFile\n";
    exit;
}
if(!(-e "$pfbFile")){
    print "Couldn't find pfb file $pfbFile\n";
    exit;
}
if(!(-d "$genotypeFilesDir")){
    print "Couldn't find the $genotypeFilesDir\n";
    exit;
}
if(!(-d "$resultsDir")){
    print "Couldn't find the $resultsDir\n";
    exit;
}

$fileCount = 0;
foreach(@inputFileNames){
    /$genotypeFilesDir\/(.*)/;
    if(($fileCount + 1) < 10){
        $printNumber = '0' . ($fileCount+1);
    }
    else{
        $printNumber = ($fileCount+1)
    }

    my $inFile = $_;
    my $shortFileName = $1 . "-$printNumber";
    my $logFile       = "$resultsDir/$shortFileName" . ".log";
    my $rawCNVFile    = "$resultsDir/$shortFileName" . ".rawcnv";
    my $outFile       = "$resultsDir/$shortFileName" . ".out";

    if(-e $logFile)   { system("rm $logFile"); }
    if(-e $rawCNVFile){ system("rm $rawCNVFile"); }
    if(-e $outFile)   { system("rm $outFile"); }

    print "/software/team145/perl-5.18.2/bin/perl $PennCNVDir/detect_cnv.pl -test -hmm $hmmFile -pfb $pfbFile --confidence --region 1-22 --minsnp 1 -log $logFile -out $rawCNVFile $inFile > $outFile\n";
    system "/software/team145/perl-5.18.2/bin/perl $PennCNVDir/detect_cnv.pl -test -hmm $hmmFile -pfb $pfbFile --confidence --region 1-22 --minsnp 1 -log $logFile -out $rawCNVFile $inFile > $outFile";
    $fileCount++;
}

