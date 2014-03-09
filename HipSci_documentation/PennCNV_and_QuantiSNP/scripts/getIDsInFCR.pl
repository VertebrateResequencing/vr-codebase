#!/bin/perl -w

use strict;

if($#ARGV != 0){
    print "Usage: $0 <GENOTYPES FILE>\n";
    exit 0;
}

my $genotypesFile = $ARGV[0];
#my $manifestSummaryFile = $ARGV[1];
#my $thisSampleNumber = $ARGV[2];

$| = 1;

my (
    $i, 
    $j
    );
my $lastSampleName = '';
my @sampleNames;

#my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
#$year += 1900;
#$mon += 1;
#if($mon < 10){
#    $mon = '0' . $mon;
#}
#if($mday < 10){
#    $mday = '0' . $mday;
#}
#my $dateStr = "$year" . "$mon" . "$mday";
#
#print "##fileformat=VCFv4.0\n";
#print "##fileDate=$dateStr\n";
#print "##source=$genotypesFile HipSci genotyping file\n";
#print "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n";
#print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
#print "##FORMAT=<ID=GC,Number=1,Type=Float,Description=\"GenCall score\">\n";
#print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"; # using tabs as it says tab delimited in the VCF spec

open (GENOTYPES, "<", $genotypesFile) || die "cannot open $genotypesFile ($!)";
$_ = <GENOTYPES>;
while(<GENOTYPES>){
    my @columns = split;
    my $sampleName = $columns[1];
    if($lastSampleName ne $sampleName){
        $lastSampleName = $sampleName;
        push(@sampleNames, $sampleName);
    }
}
close(GENOTYPES);

my $thisSampleName;
#my $thisSampleName = '271298_A01_hipscigt5466706';
my $numberOfSamples = ($#sampleNames) + 1;
$i = 1;
foreach(@sampleNames){
    #if($i == $thisSampleNumber){
        print "$_\n";        
    #    $thisSampleName = $_;
    #}
    #$i++;
}
#print "$thisSampleName\n";
##exit;
##Name	Sample ID	Allele1 - Top	Allele2 - Top	GC Score	Theta	R	X	Y	X Raw	Y Raw	Pilot.B Allele Freq	Pilot.Log R Ratio
#my $genotypeFileCount = 0;
#
#my $arrCount = 0;
#$i = 0;
#my @manifestSummaryArr;
##print "$manifestSummaryFile:";
#open(MANIFESTSUMMARY, "<", $manifestSummaryFile) || die "cannot open $manifestSummaryFile ($!)";
#while(<MANIFESTSUMMARY>){
#    $i++;
#    #print "$i : $_";
#    chomp;
#    my @manifestSummaryVals = split (/\t/);
#    if($manifestSummaryVals[5] eq 'UNKNOWN'){
#        next;
#    }
#    else{
#        $manifestSummaryArr[$arrCount][0] = $manifestSummaryVals[0];
#        $manifestSummaryArr[$arrCount][1] = $manifestSummaryVals[1];
#        $manifestSummaryArr[$arrCount][2] = $manifestSummaryVals[2];
#        $manifestSummaryArr[$arrCount][3] = $manifestSummaryVals[3];
#        $manifestSummaryArr[$arrCount][4] = $manifestSummaryVals[4];
#        $manifestSummaryArr[$arrCount][5] = $manifestSummaryVals[5];
#        $arrCount++;
#    }
#}
#close(MANIFESTSUMMARY);
##print "Read in $arrCount from $i\n";
##exit;
#
#my $lineCount = 0;
#open (GENOTYPES, "$genotypesFile") || die "cannot open $genotypesFile ($!)";
#<GENOTYPES>;
#while(<GENOTYPES>){
#    $genotypeFileCount++;
#    chomp;
#    my @vals = split(/\t/);
#    my ($ref, $alt);
#    my $indelType;
#    if($thisSampleName eq $vals[1]){
#        $lineCount++;
#        if($lineCount < 300000){
#            my $genotypeSNPName = $vals[0];
#            my $genotypeSampleID = $vals[1];
#            my $genotypeAlleleOne = $vals[2];
#            my $genotypeAlleleTwo = $vals[3];
#            my $genotypeGCScore = $vals[4];
#            if(($genotypeAlleleOne eq 'D') || ($genotypeAlleleOne eq 'I')){
#                next;
#            }
#            if(($genotypeAlleleOne eq '-') || ($genotypeAlleleOne eq '-')){
#                next;
#            }
#
#            for($i = 0; $i < $arrCount; $i++){
#                my $manifestSummaryGenotypeID    = $manifestSummaryArr[$i][0];
#                if($genotypeSNPName eq $manifestSummaryGenotypeID){
#                    my $manifestSummaryIlmnStrand    = $manifestSummaryArr[$i][1];
#                    my $manifestSummarySNP           = $manifestSummaryArr[$i][2];
#                    my $manifestSummaryChr           = $manifestSummaryArr[$i][3];
#                    my $manifestSummaryMapInfo       = $manifestSummaryArr[$i][4];
#                    my $manifestSummarySourceStrand  = $manifestSummaryArr[$i][5];
#        #}
#        #open(MANIFESTSUMMARY, "$manifestSummaryFile") || die "cannot open $manifestSummaryFile ($!)";
#        #while(<MANIFESTSUMMARY>){
#            #chomp;
#            #my @manifestSummaryVals = split (/\t/);
#            #my $manifestSummaryGenotypeID    = $manifestSummaryVals[0];
#            #if($genotypeSNPName eq $manifestSummaryGenotypeID){
#                #my $manifestSummaryIlmnStrand    = $manifestSummaryVals[1];
#                #my $manifestSummarySNP           = $manifestSummaryVals[2];
#                #my $manifestSummaryChr           = $manifestSummaryVals[3];
#                #my $manifestSummaryMapInfo       = $manifestSummaryVals[4];
#                #my $manifestSummarySourceStrand  = $manifestSummaryVals[5];
#                    if($manifestSummarySourceStrand eq 'UNKNOWN'){
#                        last;
#                    }
#                    print "$manifestSummaryChr\t";
#                    print "$manifestSummaryMapInfo\t";
#                    print "$manifestSummaryGenotypeID\t";
#                    my $firstRefAllele;
#                    my $secondRefAllele;
#                    $_ = $manifestSummarySNP;
#                    /\[([A-Za-z]+)\/([A-Za-z]+)/;
#                    $firstRefAllele = $1;
#                    $secondRefAllele = $2;
################## DO THE ILMN FLIP ON THE VCF SEQ
#                    my $flippedFirstRefAllele = '';
#                    my $flippedSecondRefAllele = '';
#                    if($manifestSummaryIlmnStrand eq 'BOT'){
#                        # FLIP TO COMPLEMENTS
#                        $flippedFirstRefAllele = &flipAllele($firstRefAllele);
#                        $flippedSecondRefAllele = &flipAllele($secondRefAllele);
#                        $firstRefAllele = $flippedFirstRefAllele;
#                        $secondRefAllele = $flippedSecondRefAllele;
#                    }
#                    my $matchOne = 0;
#                    my $matchTwo = 0;
#                    if($genotypeAlleleOne eq $firstRefAllele){
#                        $matchOne = 1;
#                    }
#                    if($genotypeAlleleTwo eq $firstRefAllele){
#                        $matchTwo = 1;
#                    }
#                    my $forwardStrandRefOne;
#                    my $forwardStrandRefTwo;
#                    if($manifestSummarySourceStrand eq '-'){
#                        $forwardStrandRefOne = &flipAllele($firstRefAllele);
#                        $forwardStrandRefTwo = &flipAllele($secondRefAllele);
#                    }
#                    else{
#                        $forwardStrandRefOne = $firstRefAllele;
#                        $forwardStrandRefTwo = $secondRefAllele;
#                    }
#                    print "$forwardStrandRefOne\t$forwardStrandRefTwo\t";
#                    print ".\t.\tNS=1\t";
#                    print "GT:GC\t";
#                    print "$matchOne/$matchTwo:$genotypeGCScore\n";
#                    last;
#                }
#            }
#        }
#        #close (MANIFESTSUMMARY);
#    }
#}
#close (GENOTYPES);
#
#sub flipAllele(){
#    my $localAllele = $_[0];
#    my $localFlippedAllele;
#    for ($j = 0; $j < length($localAllele); $j++) {
#        my $thisBase = substr($localAllele, $j, 1);
#        if($thisBase eq "A"){
#            $localFlippedAllele = $localFlippedAllele . 'T';
#        }
#        elsif($thisBase eq "C"){
#            $localFlippedAllele = $localFlippedAllele . 'G';
#        }
#        elsif($thisBase eq "G"){
#            $localFlippedAllele = $localFlippedAllele . 'C';
#        }
#        elsif($thisBase eq "T"){
#            $localFlippedAllele = $localFlippedAllele . 'A';
#        }
#        elsif($thisBase eq "a"){
#            $localFlippedAllele = $localFlippedAllele . 't';
#        }
#        elsif($thisBase eq "c"){
#            $localFlippedAllele = $localFlippedAllele . 'g';
#        }
#        elsif($thisBase eq "g"){
#            $localFlippedAllele = $localFlippedAllele . 'c';
#        }
#        elsif($thisBase eq "t"){
#            $localFlippedAllele = $localFlippedAllele . 'a';
#        }
#        else{
#            $localFlippedAllele = $localFlippedAllele . $thisBase;
#            print "Flipping unexpected $thisBase\n";
#        }
#    }
#    return $localFlippedAllele;
#}


