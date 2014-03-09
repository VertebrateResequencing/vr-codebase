#!/bin/perl -w

use strict;

if($#ARGV != 0){
    print "Usage: $0 <QUANTILE NORMALISED GENE EXPRESSION FILE>\n";
    exit 0;
}

my(
    $i, 
    $numberOfSamples, 
    @splitValues, 
    $inputFile, 
    $averaged
    );

$inputFile = $ARGV[0];

open (PROFILE_FILE, "$inputFile") || die "cannot open $_ ($!)";
$_ = <PROFILE_FILE>;
chomp;
@splitValues = split /\t/;
$numberOfSamples = ($#splitValues - 1) / 8;
for($i = 0; $i < $numberOfSamples; $i++){
    if($i == 0){
        print "GENE NAME\t";
    }
    my $sampleName = $splitValues[(($i*8)+2)];
    $sampleName =~ s/^MIN_Signal-//;
    print "$sampleName\t";
    if($i == ($numberOfSamples - 1)){
        print "\n";
    }
}
close (PROFILE_FILE);

print "Pluripotent associated:\n";
getGeneInfo("POU5F1");
getGeneInfo("SOX2");
getGeneInfo("NANOG");
print "\n";

print "Endoderm associated:\n";
getGeneInfo("SOX17");
getGeneInfo("FOXA2");
getGeneInfo("GATA4");
print "\n";

print "Mesoderm associated:\n";
getGeneInfo("T");
getGeneInfo("EOMES");
getGeneInfo("MIXL");
print "\n";

print "Neuroectoderm associated:\n";
getGeneInfo("SOX1");
getGeneInfo("NES");
print "\n";

exit;

sub getGeneInfo {

    my $geneName = $_[0];
    my @averages;
    my $geneOccurrencesCount = 0;
    my $lineCount = 0;

    open (PROFILE_FILE, "$inputFile") || die "cannot open $_ ($!)";
    while(<PROFILE_FILE>){
        chomp;
        @splitValues = split /\t/;
        $numberOfSamples = ($#splitValues - 1) / 8;
        $lineCount++;

        if($lineCount == 1){
            for($i = 0; $i < $numberOfSamples; $i++){
                $averages[$i] = 0;
            }
        }

        if($splitValues[0] eq "$geneName"){
            $geneOccurrencesCount++;

            for($i = 0; $i < $numberOfSamples; $i++){
                if($i == 0){
                    if($geneName eq "POU5F1"){
                        print "OCT4/POU5F1\t";
                    }
                    elsif($geneName eq "T"){
                        print "BRACHYURY/T\t";
                    }
                    elsif($geneName eq "NES"){
                        print "NESTIN/NES\t";
                    }
                    else{
                        print "$splitValues[0]\t";
                    }
                }
                print "$splitValues[(($i*8)+3)]\t";
                if($i == ($numberOfSamples - 1)){
                    print "\n";
                }
                $averages[$i] = $averages[$i] + $splitValues[(($i*8)+3)];
            }
        }
    }

    if($geneOccurrencesCount > 1){
        for($i = 0; $i < $numberOfSamples; $i++){
            if($i == 0){
                if($geneName eq "POU5F1"){
                    print "OCT4/POU5F1 Average\t";
                }
                elsif($geneName eq "T"){
                    print "BRACHYURY/T Average\t";
                }
                elsif($geneName eq "NES"){
                    print "NESTIN/NES Average\t";
                }
                else{
                    print "$geneName Average\t";
                }
            }
            my $averagedValue = ($averages[$i] / $geneOccurrencesCount);
            my $rounded = sprintf("%.1f", $averagedValue);
            print "$rounded\t";
            if($i == ($numberOfSamples - 1)){
                print "\n";
            }
        }
    }
    close (PROFILE_FILE);
}
