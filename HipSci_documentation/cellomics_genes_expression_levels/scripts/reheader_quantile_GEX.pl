#!/bin/perl -w

use strict;

if($#ARGV != 1){
    print "Usage: $0 <ID MAP FILE> <QUANTILE NORMALISED GEX FILE>\n";
    exit 0;
}

my $IDMapFile = $ARGV[0];
my $GEXFile = $ARGV[1];
my (
    $i, $j, 
    @sampleNames
    );
my %IDMapHash;

open (IDs, "<", $IDMapFile) || die "cannot open $IDMapFile ($!)";
while(<IDs>){
    chomp;
    /(\w+)\s+(\w+)/;
    $IDMapHash {"$1"} = "$2";
}
close(IDs);

open (GEX, "<", $GEXFile) || die "cannot open $GEXFile ($!)";
$_ = <GEX>;
chomp;
my @headers = split(/\t/);

foreach(@headers){
    if($_ eq "TargetID"){
        print "TargetID\t";
        next;
    }
    if($_ eq "ProbeID"){
        print "ProbeID\t";
        next;
    }
    my $found = 0;
    while (my ($key, $value) = each %IDMapHash){        
        /(.*)\-(.*)/;
        my $thisHeader = $1;
        my $thisID = $2;
        if(!$thisID){
            print "$thisID not found: '$_'\n";
            exit;
        }
        if($key eq $thisID){
            $found = 1;
            print "$thisHeader-$IDMapHash{$key}\t";
        }
    }
    if($found == 0){
        print "\n\nHeader not found\n";
        exit;
    }
}
print "\n";

while(<GEX>){
    print;
};

close(GEX);

