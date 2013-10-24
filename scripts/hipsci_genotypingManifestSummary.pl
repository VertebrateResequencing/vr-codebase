#!/bin/perl -w

use strict;

if($#ARGV != 1){
    print "Usage: $0 <MANIFEST FILE> <STRAND FILE>\n";
    exit 0;
}

# 1). Take SNP info from the Illumina HumanCoreExome manifest file (humancoreexome-12v1-0_a.csv) which 
#     contains info about the SNPs on this genotyping chip: extract 5 fields = SNP ID, Ilumina strand, SNP, chr, coord
# 2). Take strand info from Will Rayner's strand file (humancoreexome-12v1-0_b-b37.strand): 1 field = strand (+ or -)
#     N.B. see http://www.well.ox.ac.uk/~wrayner/strand/ for further info
# 3). Combine fields from 1) and 2) to make a summary of the manifest + strand files.
#
# E.g. 200610-10	bot	[t/c]	mt	6753	+ 
# (where + came from the strand file).
#
# This information is needed to generate the VCF. In addition to the SNP info (SNP ID, SNP, chr, coord) there 
# are two fields that show if the SNP should be flipped during later VCF generation (by hipsci_genotypingfcr2vcf.pl): 
# a). Ilmn strand i.e. if "bot", then flip the SNP to represent the "top" strand;
# b). Take the result from a) and if Will Rayner's strand file indicates "-" then flip to have the final representation 
#     for the human forward strand. 
# I.e. the rule given by Hannah Blackburn from Sanger genotyping is:
# FLIP IF:
# ILMNStrand    RefStrand
# TOP           -
# BOT           + 

#EXAMPLE USAGE:
#perl makeCoreReference_SPLIT.pl HumanCoreExome-12v1-0_A.csv HumanCoreExome-12v1-0_B-b37.strand > VRSEQ-REF_HumanCoreExome-12v1-0_A.csv

my (
    $manifestFileIter, 
    $manifestFile, 
    $strandFile, 
    @manifestVals, 
    $manifestGenotypeID, 
    $manifestIlmnStrand, 
    $manifestSNP, 
    $manifestChr, 
    $manifestMapInfo, 
    $found, 
    @strandTmpArr
    );

$manifestFile = $ARGV[0];    # E.g. HumanCoreExome-12v1-0_A.csv
$strandFile = $ARGV[1];      # E.g. HumanCoreExome-12v1-0_B-b37.strand

$| = 1; # AUTOFLUSH ON

$manifestFileIter = 0;

open(MANIFEST, "$manifestFile") || die "cannot open $manifestFile ($!)";
<MANIFEST>;
<MANIFEST>;
<MANIFEST>;
<MANIFEST>;
<MANIFEST>;
<MANIFEST>;
<MANIFEST>;
while(<MANIFEST>){
    if(/^\[Controls\]/){ # END OF MANIFEST INFO
        last;
    }
    $manifestFileIter++;

    chomp;
    @manifestVals = split (/,/);
    $manifestGenotypeID = $manifestVals[1];
    $manifestIlmnStrand = $manifestVals[2];
    $manifestSNP        = $manifestVals[3];
    $manifestChr        = $manifestVals[9];
    $manifestMapInfo    = $manifestVals[10];

    if(($manifestIlmnStrand eq "TOP") || ($manifestIlmnStrand eq "BOT")){
        print "$manifestGenotypeID\t";
        print "$manifestIlmnStrand\t";
        print "$manifestSNP\t";
        print "$manifestChr\t";
        print "$manifestMapInfo\t";
        # GET THE HUMAN FORWARD STRAND INFO I.E. TO FLIP OR NOT FROM ILMN
        $found = 0;
        open(STRANDFILE, "$strandFile");
        while(<STRANDFILE>){
            @strandTmpArr = split(/\t/);
            if($manifestGenotypeID eq $strandTmpArr[0]){
                print "$strandTmpArr[4]\n";
                $found++;
                last;
            }
        }
        close(STRANDFILE);
        if($found == 0){
            print "UNKNOWN\n";
        }
    }
}
close (MANIFEST);
