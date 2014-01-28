#!/usr/bin/env perl

use Getopt::Long;

use strict;
use warnings;
no warnings 'uninitialized';
use Data::Dumper;

my ( $PROFILE_FILE, $ANNOS_FILE, $MAPPING_FILE, $numberOfSamples, $OUTPUT_FILE, $help );

GetOptions(
    'p|profile=s'               => \$PROFILE_FILE,
    'a|annot=s'                 => \$ANNOS_FILE,
    'm|mapping=s'               => \$MAPPING_FILE,
    's|samples=i'               => \$numberOfSamples,
    'o|out=s'                   => \$OUTPUT_FILE,
    'h|help'                    => \$help,
);

( $PROFILE_FILE && $ANNOS_FILE && $MAPPING_FILE && $numberOfSamples && $OUTPUT_FILE && !$help ) or die <<USAGE;
Usage: $0   
  -p|--profile                <Genome Studio profile file>
  -a|--annot                  <Genome Studio annotation file>
  -m|--mapping                <Sample mapping file>
  -s|--samples                <Number of samples to reformat>
  -o|--out                    <Output file>
  -h|--help                   <this message>
USAGE


my %annothash;
my %idhash;

open OUTF, ">", $OUTPUT_FILE || die "cannot open $OUTPUT_FILE ($!)";

open (ANNOS, "$ANNOS_FILE") || die "cannot open $ANNOS_FILE ($!)";
<ANNOS>;
while ( <ANNOS> ) {
    chomp;
    my @arrayvals = split /\t/;
    $annothash{"$arrayvals[0]_$arrayvals[1]"} = $arrayvals[4];
}
close (ANNOS);

open (IDS, "$MAPPING_FILE") || die "cannot open $MAPPING_FILE ($!)";
while(<IDS>){
    chomp;
    my @arrayvals = split /\t/;
    $arrayvals[1] =~ s/ //g;
    $idhash{$arrayvals[0]} = $arrayvals[1];
}
close (IDS);

open (PROFILE, "$PROFILE_FILE") || die "cannot open $PROFILE_FILE ($!)";
$_ = <PROFILE>;
chomp;
my @arrayvals = split /\t/;

my $valsCount = $#arrayvals;
my $sampleCount = (($valsCount - 2) + 1) / 8;
my $indexStart;


#1
#2
#3#MIN_Signal-9257622034_A	
#4#AVG_Signal-9257622034_A	
#5#MAX_Signal-9257622034_A	
#6#NARRAYS-9257622034_A	
#7#ARRAY_STDEV-9257622034_ABEAD_STDEV-9257622034_A	
#8#Avg_NBEADS-9257622034_A	
#9#Detection-9257622034_A
#3  7  8  9
#11 15 16 17

#PRINT OUT THE HEADER:

$sampleCount = $numberOfSamples;

my $headerline = "$arrayvals[1]\t";
for(my $k = 0; $k < $sampleCount; $k++){
    $indexStart = ($k * 8) + 2;
    my $headerOne = generateHeaderString(\%idhash, $arrayvals[($indexStart + 1)]);
    my $headerTwo = generateHeaderString(\%idhash, $arrayvals[($indexStart + 5)]);
    my $headerThree = generateHeaderString(\%idhash, $arrayvals[($indexStart + 6)]);
    $arrayvals[($indexStart + 7)] =~ s/^Detection/Detection.Pval/;
    my $headerFour = generateHeaderString(\%idhash, $arrayvals[($indexStart + 7)]);
    $headerline = $headerline . "$headerOne\t$headerTwo\t$headerThree\t$headerFour\t";
}
$headerline =~ s/\s+$//;
print OUTF $headerline, "\n";

sub generateHeaderString {
	my ($idhash, $header_string) = @_;
	my %IDhash = %{ $idhash };
	my $searchtag = (split '-', $header_string)[1];
	$header_string =~ s/-(.*?)_(.*)/-$idhash{$searchtag}/ if defined $idhash{$searchtag};
	$header_string =~ s/\s+$//;
	return $header_string;
}
	
#PRINT OUT THE VALUES:

while(<PROFILE>){
    chomp;
    my @arrayvals = split /\t/;
    
    my $searchtag = "$arrayvals[0]_$arrayvals[1]";
    if ( defined $annothash{$searchtag} ){
		print OUTF "$annothash{$searchtag}\t";
    }
    else {
        print OUTF "Couldn't find $arrayvals[0]\n";
        exit;
    }

    my $lineout = '';
    for(my $k = 0; $k < $sampleCount; $k++){
        $indexStart = ($k * 8) + 2;
        $lineout = $lineout . "$arrayvals[($indexStart + 1)]\t$arrayvals[($indexStart + 5)]\t$arrayvals[($indexStart + 6)]\t$arrayvals[($indexStart + 7)]\t";
    }
    $lineout =~ s/\s+$//;
    print OUTF $lineout, "\n";
}
close (PROFILE);
close (OUTF);
