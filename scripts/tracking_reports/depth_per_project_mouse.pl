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

my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => 'mouse_reseq_track',
                                                          mode => 'r');

unless ($vrtrack){
    croak 'DB connection failed: '.$DBI::errstr;
}

print "Project\tTotal Raw Bases\tSampled Bases\tMapped Bases\tRmdup\tDepth Est\n";
my $projects = $vrtrack->projects();
foreach( @{$projects} )
{
	my $project = $_;
	my $name = $project->name();
	print "$name\t";
	
	my $no_qc_lanes = 0;
	my $mapped_bases = 0;
	my $total_passed_raw_bases = 0;
	my $rmdup_bases = 0;
	my $sampled_bases = 0;
	my $samples = $project->samples();
	my $no_qc_raw_bases = 0;
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
			
			foreach( @{$lanes} )
			{
				my $lane = $_;
				$name = $lane->name();
				#print "Lane: $name\n";
				
				if( $lane->qc_status() eq 'passed' && $lane->genotype_status() eq 'confirmed' )
				{
					$total_passed_raw_bases += $lane->raw_bases();
					
					my $mapping = $lane->latest_mapping();
					$mapped_bases += $mapping->bases_mapped();
					$rmdup_bases += $mapping->rmdup_bases_mapped();
					$sampled_bases += $mapping->raw_bases();
				}
				elsif( $lane->qc_status() eq 'no_qc' )
				{
					$no_qc_raw_bases += $lane->raw_bases();
					$no_qc_lanes ++;
				}
			}
		}
	}
	
	if( $total_passed_raw_bases == 0 )
	{
		print "\n";
		next;
	}
	
	print "$total_passed_raw_bases\t";
	print "$sampled_bases\t";
	print "$mapped_bases\t";
	print "$rmdup_bases\t";
	my $depth = ( ( $rmdup_bases / $sampled_bases ) * $total_passed_raw_bases ) / 3000000000;
	print "$depth\t";
	$depth = ( ( ( $rmdup_bases / $sampled_bases ) * $total_passed_raw_bases ) + $no_qc_raw_bases ) / 3000000000;
	print "$depth\n"
}

