#!/usr/bin/env perl
=pod
Take a list of studies and an existing map of individual <-> sample
and adds new samples as required OR creates new individuals as required

RULE: if an existing individual is a prefix of a new sample - then assign it to the individual
=cut

use strict;
use warnings;
use Getopt::Long;
use Carp;

use Sfind::Sfind;
use Sfind::Study;

my ($studies, $currentMap, $newMap, $newSamples, $species, $taxon, $help);

GetOptions
(
    's|studies=s'       =>  \$studies,
    'm|mapping=s'       =>  \$currentMap,
    'spe|species=s'     =>  \$species,
    't|taxon=s'         =>  \$taxon,
    'nm|new_map=s'      =>  \$newMap,
    'ns|new_sam=s'      =>  \$newSamples,
    'h|help'            =>  \$help,
);

($studies && -f $studies && $currentMap && $newMap && $species && $taxon && $taxon =~ /^\d+$/ && !$help) or die <<USAGE;
    Usage: $0   
                --studies   <file of study names>
                --mapping   <current individual - sample map>
                --taxon     <taxon id>
                --species   <species name - Mus musculus or Homo sapiens>
                --new_map   <new individual - sample map>
                --new_sam   <new samples file>
                --help      <this message>
USAGE

my %currentMap;
open( my $ifh, $currentMap ) or die $!;
while( <$ifh> )
{
    chomp;
    if( $_ =~ /^(.*)\t(.*)$/ )
    {
        if( $currentMap{ $1 } ){push( @{ $currentMap{ $1 } }, $2 );}else{$currentMap{ $1 } = [ $2 ];}
    }
    else{print qq[incorrectly formatted entry: $_;\n];exit;}
}
close( $ifh );

my $sfind = Sfind::Sfind->new();

open( $ifh, $studies ) or die $!;
while( <$ifh> )
{
    chomp;
    my $sname = $_;
    print qq[STUDY: $sname\n];
    my $study = $sfind->get_study_by_name( $sname );
    
    my $samples = $study->samples();

    if( ! $samples ){print qq[ERROR: Please check study exists - no samples found: $sname\n];next}

    foreach( @{ $samples } )
    {
        my $sam_name = $_->name();
        
        my $hasIndividual = 0;
        foreach( keys( %currentMap ) )
        {
            if( $sam_name =~ /^$_/ ) #################RULE HERE - IF SAMPLE IS A PREFIX OF AN INDIVIDUAL - THEN ASSIGN IT TO THE INDIVIDUAL
            {
                my $ind = $_;
                my $found = 0;
                foreach( @{ $currentMap{ $_ } } ){if( $_ eq $sam_name ){$found = 1;}}
                if( ! $found ){print qq[NEW SAMPLE: $sam_name on existing individual $ind\n];push( @{ $currentMap{ $_ } }, $sam_name );}
                $hasIndividual = 1;
                last;
            }
        }
        
        if( ! $hasIndividual )
        {
            print qq[NEW SAMPLE: $sam_name on new individual\n];
            my $ind;
            if( $sam_name =~ /^(.*)(_\d+)$/ ){$ind = $1;}else{$ind = $sam_name;}
            $currentMap{ $ind } = [ $sam_name ];
        }
    }
}

open( my $sfh, qq[>$newSamples] ) or die $!;
foreach( sort keys( %currentMap ) )
{
    print $sfh qq[$_\t\t$_\t$species\t$taxon\tM\t\n];
}
close( $sfh );

open( my $ofh, qq[>$newMap] ) or die $!;
foreach( sort keys( %currentMap ) )
{
    my $k = $_;
    foreach( @{ $currentMap{ $k } } )
    {
        print $ofh qq[$k\t$_\n];
    }
}
close( $ofh );
