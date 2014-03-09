#!/bin/perl -w

use strict;

if($#ARGV != 1){
    print "Usage: $0 <CONVERTED TO QSNP DIR> <RUN QSNP DIR>\n";
    exit 0;
}

#USAGE EXAMPLE:
#perl bsub_quantisnp.pl converted_to_QSNP_dirs/converted/ run_QSNP_dirs/

my(
    $convertedToQSNPDir, 
    $runQSNPDir, 
    @convertedFiles, 
    $QSNPDir, 
    $inFile, 
    $shortFileName, 
    $shFile, 
    $lsfFileO, 
    $lsfFileE, 
    $outFile
);

$convertedToQSNPDir = $ARGV[0];
$runQSNPDir = $ARGV[1];

$QSNPDir = '/software/vertres/installs/HipSci/QuantiSNP';

@convertedFiles = <$convertedToQSNPDir/*.converted>;

foreach(@convertedFiles){
    $inFile = $shortFileName = '';
    $inFile = $_;
    /$convertedToQSNPDir\/(.*)\.converted/;
    $shortFileName = $1;
    if(!(-d $convertedToQSNPDir)){
        print "Couldn't find $convertedToQSNPDir\n";
        exit;
    }
    if(!(-d $runQSNPDir)){
        print "Couldn't find $runQSNPDir\n";
        exit;
    }
    if(!(-d "$runQSNPDir/sh_files/")){
        print("mkdir $runQSNPDir/sh_files/");
        system("mkdir $runQSNPDir/sh_files/");
    }
    if(!(-d "$runQSNPDir/lsf_files/")){
        print("mkdir $runQSNPDir/lsf_files/");
        system("mkdir $runQSNPDir/lsf_files/");
    }
    if(!(-d "$runQSNPDir/results/")){
        print("mkdir $runQSNPDir/results/");
        system("mkdir $runQSNPDir/results/");
    }    

    $shFile = "$runQSNPDir" . "/sh_files/" . $shortFileName . ".sh";
    $lsfFileO = "$runQSNPDir" . "/lsf_files/" . $shortFileName . ".o";
    $lsfFileE = "$runQSNPDir" . "/lsf_files/" . $shortFileName . ".e";
    $outFile = "$runQSNPDir" . "/results/" . $shortFileName . ".out";

    if(-e $shFile){ system("rm $shFile"); }
    if(-e $lsfFileO){ system("rm $lsfFileO"); }
    if(-e $lsfFileE){ system("rm $lsfFileE"); }
    if(-e $outFile){ system("rm $outFile"); }

    open (SH_FILE, ">$shFile") || die "cannot open $shFile ($!)";
    print SH_FILE "#!/bin/sh\n\n";
    print SH_FILE "sh $QSNPDir/run_quantisnp2.sh $QSNPDir/v79 --outdir $runQSNPDir/results --input-files $inFile --sampleid $shortFileName --levels $QSNPDir/config/levels.dat --config $QSNPDir/config/params.dat --chr 1:22 > $outFile\n"; 
    close(SH_FILE);

    print("bsub -q long -o $lsfFileO -e $lsfFileE -R \"select[type==X86_64 && mem > 16000] rusage[mem=16000]\" -M16000000 \"sh $shFile\"\n");
    system("bsub -q long -o $lsfFileO -e $lsfFileE -R \"select[type==X86_64 && mem > 16000] rusage[mem=16000]\" -M16000000 \"sh $shFile\"");
}

