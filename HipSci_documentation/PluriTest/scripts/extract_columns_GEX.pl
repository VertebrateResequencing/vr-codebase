#!/bin/perl -w

use strict;

if($#ARGV != 0){
    print "Usage: $0 <INPUT FILE>\n";
    exit 0;
}

my(
    $i, 
    $count, 
    @VALS, 
    $inputFile
    );

$inputFile = $ARGV[0];

open (PROFILE_FILE, "$inputFile") || die "cannot open $_ ($!)";
while(<PROFILE_FILE>){
    chomp;
    @VALS = split /\t/;
    $count = $#VALS / 4;
    my ($colOne, $colTwo, $colThree, $colFour);

#my @order = (16, 3, 7, 5);        # cehw      
#my @order = (24, 2, 9, 4);        # giuo      
#my @order = (31, 22, 27, 25);     # hikj      
#my @order = (29, 21, 23);         # iakz      
#my @order = (28, 30, 32, 26);     # leeh      
#my @order = (20, 15, 12, 13);     # ougl      
#my @order = (17, 11, 19, 10, 14); # peop      
#my @order = (18, 6, 8);           # qonr      

#set 14:
#my @order = (6, 11, 2, 4);        # dard      
#my @order = (5, 1, 9, 3);         # vorx      
my @order = (8, 12, 10, 7);        # veqz      

    my $listCount = 0;
    print "$VALS[0]\t";
    foreach(@order){
        $i = $_;
        $i = $i - 1;
        $colOne   = ($i * 4) + 1;
        $colTwo   = $colOne + 1;
        $colThree = $colOne + 2;
        $colFour  = $colOne + 3;
        print "$VALS[$colOne]\t$VALS[$colTwo]\t$VALS[$colThree]\t$VALS[$colFour]";
        $listCount++;
        if($listCount <= $#order){
            print "\t";
        }
    }
    print "\n";
}
close (PROFILE_FILE);

exit;
