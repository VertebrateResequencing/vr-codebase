#!/usr/bin/env perl
# 
# Author:        tk2
# Maintainer:    tk2
# Created:       
# Updated:       

use Carp;
use strict;
use warnings;

use VertRes::Utils::VRTrackFactory;
use VRTrack::Project;
use VRTrack::Sample;
use VRTrack::Library;
use VRTrack::Lane;

my %db_for_spp = ( 'mouse'  => 'mouse_reseq_track',
                    'g1k'   => 'g1k_track',
                  );

my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => 'g1k_track',
                                                          mode => 'r');

unless ($vrtrack){
    croak 'DB connection failed: '.$DBI::errstr;
}

my $projects = $vrtrack->projects();
foreach( @{$projects} )
{
	my $project = $_;
	my $name = $project->name();
	print "$name\n";
	my $samples = $project->samples();
	foreach( @{$samples} )
	{
		my $sample = $_;
		$name = $sample->name();
		#print "Sample: $name\n";
		my $libraries = $sample->libraries();
		
		foreach( @{$libraries} )
		{
			my $library = $_;
			$name = $library->name();
			#print "Library: $name\n";
			my $lanes = $library->lanes();
			
			if( @$lanes == 0 )
			{
				print "$name has 0 lanes\n";
			}
		}
	}
	

}
