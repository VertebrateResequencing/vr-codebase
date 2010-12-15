#!/usr/bin/env perl
use strict;
use warnings;
no warnings 'uninitialized';
use Getopt::Long;
use VertRes::Utils::VRTrackFactory;
use VRTrack::Library;

my ($spp, $lib, $help);

GetOptions(
    's|spp=s'       =>  \$spp,
    'l|lib=s'       =>  \$lib,
    'h|help'	    =>  \$help,
    );

my %db_for_spp = ( 'mouse'  => 'mouse_reseq_track',
		    'g1k'   => 'g1k_track',
	  );

# For testing
#my %db_for_spp = ( 'mouse'  => 'mouse_cancer_track',
#		    'g1k'   => 'mouse_cancer_track',
#		  );

my $db = $db_for_spp{$spp};

($db && !$help) or die <<USAGE;
    Usage: $0 
                --spp       <species, i.e. g1k or mouse>
                --lib       <name of library to close>
                --help      <this message>

Closes a library in VRTrack database.  Does not delete any information.

USAGE

warn "Species: $spp\n";
warn "Database: $db\n";
my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $db,
                                                          mode => 'rw');

unless ($vrtrack){
    die "Can't connect to tracking database\n";
}

my $library = VRTrack::Library->new_by_name($vrtrack,$lib);
$library = VRTrack::Library->new_by_hierarchy_name($vrtrack,$lib) unless $library;
unless ($library){
    die "Can't get library $lib\n";
}

$library->open(0);
unless ($library->update){
    die "Can't close library $lib\n";
}
