#!/usr/bin/env perl

use warnings;
use strict;

if($#ARGV != 2){
    print "Usage: $0 <GENOTYPES FILE> <MANIFEST SUMMARY FILE> <SAMPLE NUMBER>\n";
    exit 0;
}
 
# VCF GENERATION FROM GENOTYPE FILE USING MANIFEST SUMMARY FILE 
# (MANIFEST SUMMARY FILE IS PREVIOUSLY PRODUCED FROM MANIFEST + STRAND FILES)

# USAGE EXAMPLE:
# perl make_VCF-parts.pl pilot.fcr.txt-EDIT VRSEQ-REF_HumanCoreExome-12v1-0_A.csv 3 > VCF_03.txt
# I.e. convert the 3rd sample in the FCR file (samples are in sequential chunks of results in the FCR)

my (
    $sampleIter, 
    $lastSampleName, 
    @sampleNames, 
    $sampleName, 
    $thisSampleName, 
    $genotypesFile, 
    $manifestSummaryFile, 
    $thisSampleNumber, 
    @columns, 
    @manifestSummaryVals, 
    @tmpArr, 
    $genotypeSNPName, 
    $genotypeAlleleOne, 
    $genotypeAlleleTwo, 
    $genotypeGCScore, 
    $manifestSummaryGenotypeID, 
    $manifestSummaryIlmnStrand, 
    $manifestSummarySNP, 
    $manifestSummaryChr, 
    $manifestSummaryMapInfo, 
    $manifestSummarySourceStrand, 
    $firstRefAllele, 
    $secondRefAllele, 
    $flippedFirstRefAllele, 
    $flippedSecondRefAllele, 
    $matchOne, 
    $matchTwo, 
    $forwardStrandRefOne, 
    $forwardStrandRefTwo, 
    $thisSampleNameCount
    );

$genotypesFile = $ARGV[0];       # E.g. pilot.fcr.txt-EDIT
$manifestSummaryFile = $ARGV[1]; # E.g. VRSEQ-REF_HumanCoreExome-12v1-0_A.csv
$thisSampleNumber = $ARGV[2];    # E.g. 3

$| = 1; # AUTOFLUSH ON

# GET THE DATE FOR THE VCF HEADER
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900;
$mon += 1;
if($mon < 10){
    $mon = '0' . $mon;
}
if($mday < 10){
    $mday = '0' . $mday;
}
my $dateStr = "$year" . "$mon" . "$mday";

# VCF HEADER
print "##fileformat=VCFv4.0\n";
print "##fileDate=$dateStr\n";
print "##source=$genotypesFile HipSci genotyping file\n";
print "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n";
print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
print "##FORMAT=<ID=GC,Number=1,Type=Float,Description=\"GenCall score\">\n";
print "##FORMAT=<ID=IA,Number=1,Type=Float,Description=\"Intensity of the A Allele\">\n";
print "##FORMAT=<ID=IB,Number=1,Type=Float,Description=\"Intensity of the B Allele\">\n";
print "##FORMAT=<ID=BAF,Number=1,Type=Float,Description=\"B Allele Frequency\">\n";
print "##FORMAT=<ID=LRR,Number=1,Type=Float,Description=\"Log R Ratio\">\n";
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";

# GET ALL THE SAMPLE NAMES IN THE GT FILE
$lastSampleName = '';
open (GENOTYPES, "<", $genotypesFile) || die "cannot open $genotypesFile ($!)";
$_ = <GENOTYPES>;
while(<GENOTYPES>){
    @columns = split;
    $sampleName = $columns[1];
    if($lastSampleName ne $sampleName){
        $lastSampleName = $sampleName;
        push(@sampleNames, $sampleName);
    }
}
close(GENOTYPES);

# PRINT THE SAMPLE NAME FOR THIS SAMPLE NUMBER
$sampleIter = 1;
foreach(@sampleNames){
    if($sampleIter == $thisSampleNumber){
        $thisSampleName = $_;
    }
    $sampleIter++;
}
#$thisSampleName = '271298_A01_hipscigt5466706';
print "$thisSampleName\n";

# READ IN THE MANIFEST SUMMARY INFO
# In manifest:
# 200610-10	BOT	[T/C]	MT	6753	+
# In GT:
# 1KG_1_100177980	271298_A01_hipscigt5466706	D	D	0.4209	0.034	1.310	1.242	0.067	13158	564	0.0052	-0.0193

my %manifestSummaryHoH = ();

open(MANIFESTSUMMARY, "<", $manifestSummaryFile) || die "cannot open $manifestSummaryFile ($!)";
while(<MANIFESTSUMMARY>){
    chomp;
    @manifestSummaryVals = split (/\t/);
    if($manifestSummaryVals[5] eq 'UNKNOWN'){
        next;
    }
    else{
        $manifestSummaryHoH{ $manifestSummaryVals[0] } = { # TO DO: NEED TO MAKE SURE THERE ARE NO DUPLICATES IN THE INPUT?
            strand  => "$manifestSummaryVals[1]",
            allele  => "$manifestSummaryVals[2]",
            chr     => "$manifestSummaryVals[3]",
            coord   => "$manifestSummaryVals[4]",
            flip    => "$manifestSummaryVals[5]",
        }
    }
}
close(MANIFESTSUMMARY);

open (GENOTYPES, "$genotypesFile") || die "cannot open $genotypesFile ($!)";
<GENOTYPES>;
while(<GENOTYPES>){
    chomp;
    @tmpArr = split(/\t/);
    if($thisSampleName eq $tmpArr[1]){
        $genotypeSNPName = $genotypeAlleleOne = $genotypeAlleleTwo = $genotypeGCScore = '';
        $genotypeSNPName = $tmpArr[0];
        $genotypeAlleleOne = $tmpArr[2];
        $genotypeAlleleTwo = $tmpArr[3];
        $genotypeGCScore = $tmpArr[4];
        if(($genotypeAlleleOne eq 'D') || ($genotypeAlleleOne eq 'I')){
            next;
        }
        if(($genotypeAlleleOne eq '-') || ($genotypeAlleleOne eq '-')){
            next;
        }
        my $intensityA = $tmpArr[7];
        my $intensityB = $tmpArr[8];
        my $BAF = $tmpArr[11];
        my $LRR = $tmpArr[12];
        # FIND THE SNP IN THE MANIFEST INFO
        my $hash_ref = $manifestSummaryHoH{$genotypeSNPName} || die "missing $genotypeSNPName from summary hash?!/n";
        $manifestSummaryIlmnStrand    = $hash_ref->{strand};
        $manifestSummarySNP           = $hash_ref->{allele};
        $manifestSummaryChr           = $hash_ref->{chr};
        $manifestSummaryMapInfo       = $hash_ref->{coord};
        $manifestSummarySourceStrand  = $hash_ref->{flip};
        if($manifestSummarySourceStrand eq 'UNKNOWN'){
            next;
        }

        $manifestSummaryGenotypeID = $genotypeSNPName;        

        print "$manifestSummaryChr\t$manifestSummaryMapInfo\t$manifestSummaryGenotypeID\t";

        $firstRefAllele = $secondRefAllele = '';
        $manifestSummarySNP =~ /\[([A-Za-z]+)\/([A-Za-z]+)/;
        $firstRefAllele = $1;
        $secondRefAllele = $2;

        #FLIPPING : 1) FLIP TO THE ILMN STRAND IF NECESSARY (I.E. 'BOT' IN MANIFEST SUMMARY); 
        #           2) THEN FLIP TO THE HUMAN FORWARD STRAND IF NECESSARY (I.E. '-' IN MANIFEST SUMMARY)

        #ILMN FLIP
        $flippedFirstRefAllele = $flippedSecondRefAllele = '';
        if($manifestSummaryIlmnStrand eq 'BOT'){
            #FLIP TO COMPLEMENTS
            $flippedFirstRefAllele = &flipAllele($firstRefAllele);
            $flippedSecondRefAllele = &flipAllele($secondRefAllele);
            $firstRefAllele = $flippedFirstRefAllele;
            $secondRefAllele = $flippedSecondRefAllele;
        }
        $matchOne = $matchTwo = 0;
        #MATCH THE ALLELES AGAINST THE REF ALLELES
        if($genotypeAlleleOne eq $firstRefAllele){
            $matchOne = 1;
        }
        if($genotypeAlleleTwo eq $firstRefAllele){
            $matchTwo = 1;
        }
        #HUMAN FORWARD STRAND FLIP
        $forwardStrandRefOne = $forwardStrandRefTwo = '';
        if($manifestSummarySourceStrand eq '-'){
            $forwardStrandRefOne = &flipAllele($firstRefAllele);
            $forwardStrandRefTwo = &flipAllele($secondRefAllele);
        }
        else{
            $forwardStrandRefOne = $firstRefAllele;
            $forwardStrandRefTwo = $secondRefAllele;
        }
        print "$forwardStrandRefOne\t$forwardStrandRefTwo\t";
        print ".\t.\tNS=1\t";
        print "GT:GC:IA:IB:BAF:LRR\t";
        print "$matchOne/$matchTwo:$genotypeGCScore:$intensityA:$intensityB:$BAF:$LRR\n";
    }
}
close (GENOTYPES);

sub flipAllele(){
    my $thisBase = shift;
    my $localFlippedAllele;

    if($thisBase eq "A"){
        $localFlippedAllele = $localFlippedAllele . 'T';
    }
    elsif($thisBase eq "C"){
        $localFlippedAllele = $localFlippedAllele . 'G';
    }
    elsif($thisBase eq "G"){
        $localFlippedAllele = $localFlippedAllele . 'C';
    }
    elsif($thisBase eq "T"){
        $localFlippedAllele = $localFlippedAllele . 'A';
    }
    elsif($thisBase eq "a"){
        $localFlippedAllele = $localFlippedAllele . 't';
    }
    elsif($thisBase eq "c"){
        $localFlippedAllele = $localFlippedAllele . 'g';
    }
    elsif($thisBase eq "g"){
        $localFlippedAllele = $localFlippedAllele . 'c';
    }
    elsif($thisBase eq "t"){
        $localFlippedAllele = $localFlippedAllele . 'a';
    }
    else{
        $localFlippedAllele = $localFlippedAllele . $thisBase;
        print "Flipping unexpected $thisBase\n";
    }

    return $localFlippedAllele;
}


