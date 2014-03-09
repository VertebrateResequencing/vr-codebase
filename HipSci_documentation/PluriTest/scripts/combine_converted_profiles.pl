#!/bin/perl -w

use strict;

if($#ARGV != 1){
    print "Usage: $0 <COMBINE FILE 1> <COMBINE FILE 2>\n";
    exit 0;
}

my(
    $i, $j, $k, $l, $m, $n, 
    $count, 
    @splitValues
);

my @files = ($ARGV[0], $ARGV[1]);

my $numberOfFiles = ($#files) + 1;
my @AoA;
my $b = 1;

foreach(@files){
    my $fileName = $_;
    my @arr;
    open (PROFILE_FILE, "$_") || die "cannot open $_ ($!)";
    $_ = <PROFILE_FILE>;
    chomp;
    @splitValues = split /\t/;
    my $numberOfSamples = (($#splitValues) / 4);
    #print "There are $numberOfSamples in $fileName\n";
    close(PROFILE_FILE);

    for($i = 0; $i < $numberOfSamples; $i++){
        open (PROFILE_FILE, "$fileName") || die "cannot open $_ ($!)";
        $_ = <PROFILE_FILE>;
        chomp;
        my $index = ($i * 4);
        $j = 0;
        @splitValues = split /\t/;
        print "ProbeID\t";
        print "" . $splitValues[($index+1)] . "\t";
        print "" . $splitValues[($index+2)] . "\t";
        print "" . $splitValues[($index+3)] . "\t";
        print "" . $splitValues[($index+4)] . "\t";
        $#arr = -1;
        while(<PROFILE_FILE>){
            chomp;
            @splitValues = split /\t/;
            $arr[$j][0] = $splitValues[0];
            $arr[$j][1] = $splitValues[$index+1];
            $arr[$j][2] = $splitValues[$index+2];
            $arr[$j][3] = $splitValues[$index+3];
            $arr[$j][4] = $splitValues[$index+4];
            $j++;
        }
        close (PROFILE_FILE);
        push @AoA, [@arr];
    }
}

my $aRef;
my $currentID;
my $found;

$aRef = $AoA[0];
$n = @$aRef - 1;
for $i(0 .. $n){
    my $line = '';
    $line = "$AoA[0][$i][0]\t$AoA[0][$i][1]\t$AoA[0][$i][2]\t$AoA[0][$i][3]\t$AoA[0][$i][4]\t";
    $currentID = $AoA[0][$i][0];

    for $j(1 .. $#AoA){
	$aRef = $AoA[$j];
	$m = @$aRef - 1;
	$found = 0;
	for $k(0 .. $m){
	    if($currentID eq $AoA[$j][$k][0]){
		$line = $line . "$AoA[$j][$k][1]\t$AoA[$j][$k][2]\t$AoA[$j][$k][3]\t$AoA[$j][$k][4]";
                if($j < $#AoA){
                    $line = $line . "\t";
                }
		$found++;
	    }
	}
	if($found != 1){
	    print "ERROR on $currentID\n";
	    exit;
	}
        if($j == $#AoA){
            print "$line\n";
        }
    }
}

exit;

