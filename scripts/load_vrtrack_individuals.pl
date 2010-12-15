#!/usr/bin/env perl
use strict;
use warnings;
no warnings 'uninitialized';

use Getopt::Long;
use VertRes::Utils::VRTrackFactory;
use VRTrack::Individual;

my ($indivfile, $db, $help);

GetOptions(
    'i|indiv=s'     =>  \$indivfile,
    'd|db=s'        =>  \$db,
    'h|help'	    =>  \$help,
    );


#my $db = 'g1k_track';
#my $db = 'g1k_track_test';

(-f $indivfile && $db && !$help) or die <<USAGE;
    Usage: $0   
                --indiv     <file of individual info>
                --db        <database to update, e.g. g1k_track>
                --help      <this message>

Loads vrtrack database with new individuals.  This is required before running
update_vrtrack to pull in new samples, etc for those individuals.

Individual file is tab-delimited with fields:

1 individual name (e.g. CAST/EiJ)
2 individual alias (for genotype matching.  Can be empty)
3 population name (e.g. CAST_EiJ)
4 species name (e.g. Mus musculus castaneus)
5 taxon id (e.g. 10091)
6 sex (e.g. F)
7 sample accession (e.g. ERS000116.  Can be empty)

USAGE

print "Database: $db\n";
my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $db,
                                                          mode => 'rw');

unless ($vrtrack){
    die "Can't connect to tracking database\n";
}

open (my $INDIV, $indivfile) or die "Can't open $indivfile:$!\n";
while (<$INDIV>){
    chomp;
    my ($name, $alias, $pop, $spp, $taxon, $sex, $acc) = split "\t", $_;
    my $vind = VRTrack::Individual->new_by_name($vrtrack,$name);
   
    next if $vind; # already there, so skip out

    print "Adding individual: $name.\n";
    # need to set sex, alias, population, species.
    $vind = VRTrack::Individual->create($vrtrack, $name);
    $vind->sex($sex);
    $vind->alias($alias) if $alias;
    $vind->acc($acc) if $acc;
    my $vpop = $vind->population($pop);
    unless($vpop){
        print "\tNew population $pop\n";
        $vpop = $vind->add_population($pop);
    }
    my $vspp = $vind->species($spp);
    unless($vspp){
        print "\tNew spp ",$spp,"\n";
        $vspp = $vind->add_species($spp);
        $vspp->taxon_id($taxon) if $taxon;
        $vspp->update;
    }
    $vind->update;
}
