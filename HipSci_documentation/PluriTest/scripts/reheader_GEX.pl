#!/bin/perl -w

use strict;

if($#ARGV != 1){
    print "Usage: $0 <IDs FILE> <PROFILE FILE>\n";
    exit 0;
}

my $IDsFile = $ARGV[0];
my $GEXFile = $ARGV[1];

$| = 1;

my (
    $i, 
    $j
    );
my @sampleNames;
my %hash;

open (IDs, "<", $IDsFile) || die "cannot open $IDsFile ($!)";
while(<IDs>){
    chomp;
    /(\w+)\s+(\w+)/;
    $hash {"$1"} = "$2";
}
close(IDs);

open (GEX, "<", $GEXFile) || die "cannot open $GEXFile ($!)";
$_ = <GEX>;
my @headers = split(/\t/);

foreach(@headers){
    if($_ eq "ProbeID"){
        print "ProbeID\t";
        next;
    }
    my $found = 0;
    while (my ($key, $value) = each %hash){
        /(.*)\-(.*)/;
        my $thisHeader = $1;
        my $thisID = $2;
        if($key eq $thisID){
            $found = 1;
            print "$thisHeader-$hash{$key}\t";
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
