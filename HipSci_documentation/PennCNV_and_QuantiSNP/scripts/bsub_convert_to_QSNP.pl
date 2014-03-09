#!/bin/perl -w

use strict;

if($#ARGV != 1){
    print "Usage: $0 <SPLIT FILES DIR> <OUTPUT DIR>\n";
    exit 0;
}

#USAGE EXAMPLE:
#perl bsub_convert_to_QSNP.pl split_files/ bsub_convert_to_QSNP_dirs/ HumanCoreExome-12v1-0_A.csv

my(
    @fileNames, 
    $inFile, 
    $shortFileName, 
    $shFile, 
    $convertedFile, 
    $lsfFileO, 
    $lsfFileE, 
    $splitFilesDir, 
    $outputDir
);

$splitFilesDir = $ARGV[0];
$outputDir = $ARGV[1];

if(!(-d $splitFilesDir)){
    print "Couldn't find $splitFilesDir\n";
    exit;
}
if(!(-d $outputDir)){
    print "Couldn't find $outputDir\n";
    exit;
}

my $shDir        = "$outputDir/sh_files/";
my $convertedDir = "$outputDir/converted/";
my $lsfDir       = "$outputDir/lsf_files/";

if(!(-d $shDir)){
    system("mkdir $shDir");
}
if(!(-d $convertedDir)){
    system("mkdir $convertedDir");
}
if(!(-d $lsfDir)){
    system("mkdir $lsfDir");
}

@fileNames = <$splitFilesDir/*_*>;

foreach(@fileNames){
    /$splitFilesDir\/(.*)/;

    $inFile = $_;
    $shortFileName = $1;
    $shFile = $shDir . $shortFileName . ".sh";
    $convertedFile = $convertedDir . $shortFileName . ".converted";
    $lsfFileO = "$lsfDir" . "$shortFileName" . ".o";
    $lsfFileE = "$lsfDir" . "$shortFileName" . ".e";

    if(-e $shFile){ system("rm $shFile"); }
    if(-e $convertedFile){ system("rm $convertedFile"); }
    if(-e $lsfFileO){ system("rm $lsfFileO"); }
    if(-e $lsfFileE){ system("rm $lsfFileE"); }

    open (SH_FILE, ">$shFile") || die "cannot open $shFile ($!)";
    print SH_FILE "#!/bin/sh\n\n";
    print SH_FILE "perl convert_split_fcr_to_quanti_format.pl $inFile > $convertedFile\n"; 
    close(SH_FILE);

    print("bsub -q long -o $lsfFileO -e $lsfFileE -R \"select[type==X86_64 && mem > 16000] rusage[mem=16000]\" -M16000000 \"sh $shFile\"\n");
    #system("bsub -q long -o $lsfFileO -e $lsfFileE -R \"select[type==X86_64 && mem > 16000] rusage[mem=16000]\" -M16000000 \"sh $shFile\"");
}
