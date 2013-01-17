#!/usr/local/bin/perl -T
# 
# Displays QC Lane information with update form
#
# Author:        cj5
# Maintainer:    cj5
# Created:       2012-03-28

BEGIN {
    $ENV{VRTRACK_HOST} = 'mcs10';
    $ENV{VRTRACK_PORT} = 3306;
    $ENV{VRTRACK_RO_USER} = 'vreseq_ro';
    $ENV{VRTRACK_RW_USER} = 'vreseq_rw';
    $ENV{VRTRACK_PASSWORD} = 't3aml3ss';
};

use strict;
use warnings;
use URI;

use SangerPaths qw(core team145);
use SangerWeb;
use VRTrack::Project;
use VertRes::Utils::VRTrackFactory;

#possible states
my $PASSED = 'passed';
my $PENDING = 'pending';
my $NO_QC = 'no_qc';
my $FAILED = 'failed';
my $GT_PENDING = 'gt_pending';
my $INVESTIGATE = 'investigate';
my $CLOSE_LIBRARY = 'close';
my $OPEN_LIBRARY = 'open';

my $qc_grind_script = "./qc_grind.pl";
my $pending_view = "../pending-view/pending_view.pl";
my $LANE_VIEW = 0;
my $LIB_VIEW = 1;
my $SPECIES_VIEW = 4;
my $PROJ_VIEW = 6;

my $sw  = SangerWeb->new({
    'title'   => q(Lane QC ),
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
});

my $cgi = $sw->cgi();
my $SCRIPT_NAME = $cgi->url(-relative=>1);

my $db = $cgi->param('db');
displayError( "No database ID" ) unless $db;
my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $db, mode => 'rw');

my $projectID = $cgi->param('proj_id');
displayError( "No Project ID" ) unless $projectID;
my $project = VRTrack::Project->new( $vrtrack, $projectID );

if ($cgi->param('download')) {
    &downloadLaneData($cgi, $vrtrack, $project);
	exit;
}

print $sw->header();
displayProjectLaneForm($cgi, $vrtrack, $db, $project);
print $sw->footer();
exit;

sub displayProjectLaneForm
{
    my ($cgi, $vrtrack, $database, $project) = @_;

	my ($npg_status_filt, $auto_qc_status_filt, $bases_mapped_filt, $duplication_filt, $rmdup_mapped_filt);

	if( $cgi->param('filter')) {
		$npg_status_filt = $cgi->param('npg_status');
		$auto_qc_status_filt = $cgi->param('auto_qc_status');
		$bases_mapped_filt = $cgi->param('bases_mapped');
		$duplication_filt = $cgi->param('duplication');
		$rmdup_mapped_filt = $cgi->param('rmdup_mapped');
	}
	else {
		$npg_status_filt = 'all';
		$auto_qc_status_filt = 'all';
	}
	#use Data::Dumper;print "<pre>parms:", Dumper($cgi->param()), "FILTER $npg_status_filt  $auto_qc_status_filt  $duplication_filt  $bases_mapped_filt   </pre>"; ## DEBUG
    
    my $samples = $project->samples();
    displayError( "Cant get samples for project: $projectID" ) unless $samples;
    
    my %individuals2Samples;
    foreach( @{ $samples } ) {
        my $sample = $_;
        my $ind = $sample->individual();
        my $indname = $ind->name();
        if( $individuals2Samples{ $indname } ) {
            push( @{ $individuals2Samples{ $indname } }, $sample );
        }
        else{$individuals2Samples{ $indname } = [ $sample ];}
    }
    
    my $pname = $project->name;
    print qq[
    <h2 align="center" style="font: normal 900 1.5em arial"><a href="$qc_grind_script">QC Grind</a></h2>
    <h3 style="font: normal 700 1.5em arial"><a href="$qc_grind_script?mode=$SPECIES_VIEW&amp;db=$database">].
    ucfirst($database).qq{</a> : $pname </h3>
    <h5 style="font: arial"><a href="$pending_view?mode=2&amp;db=$database">Pending Requests</a></h5>};

    print qq[<div align=right><h5><a href="$qc_grind_script?mode=$PROJ_VIEW&amp;proj_id=$projectID&amp;db=$database">Project View</a></h5></div>];
    
    print qq[
        <div class="centerFieldset">
        <fieldset>
        <legend>$pname</legend>
        <table RULES=GROUPS>
        <tr>
		<th style="width: 40px">Pass</th>
		<th style="width: 40px">Fail</th>
		<th style="width: 40px">Invst</th>
		<th style="width: 40px">GTPend</th>
		<th style="width: 40px">Pend</th>
        <th>Individual</th>
        <th>Library</th>
        <th>Run</th>
        <th>NPG Status</th>
        <th>Auto QC Status</th>
        <th>Mapped</th>
        <th>Duplication&#37;</th>
        <th>RmDup Mapped</th>
        <th></th>
        </tr>
    ];

    print $cgi->start_form;
    print $cgi->hidden("db","$database");
    print $cgi->hidden("proj_id","$projectID");

    print qq[ <tr> ];
	print $cgi->th(['','','','','','','','']);
	print "<th>", $cgi->popup_menu( -name => 'npg_status', -values => ['all','pass','fail','pending'], -default => 'all',), "</th>"; 
	print "<th>", $cgi->popup_menu( -name => 'auto_qc_status', -values => ['all','passed','failed','no_qc'], -default => 'all',), "</th>"; 
	print "<th>", $cgi->textfield(-name=>'bases_mapped', -size=>1), "</th>";
	print "<th>", $cgi->textfield(-name=>'duplication', -size=>1), "</th>";
	print "<th>", $cgi->textfield(-name=>'rmdup_mapped', -size=>1), "</th>";
    print "<th>", $cgi->submit(-name => 'filter', -value  => 'Filter'), "</th>";
    print qq[ </tr> ];
    
    foreach( sort( keys( %individuals2Samples ) ) ) {
        my $iname = $_;
        my @samples = @{$individuals2Samples{ $_ }};

		foreach my $sample ( @samples ) {

			my @libraries = sort {$a->name cmp $b->name} @{$sample->libraries()};

# Uncomment if we want to show samples without libs
#			if (@libraries == 0) {		
#				print qq[<tr><td>$iname</td>];
#				print qq[</tr>];
#				next;
#			}

			foreach my $library ( @libraries ) {
				my $lname = $library->name;
				my $lid = $library->id();
				my $lanes = $library->lanes();

				foreach my $lane (sort {$a->name cmp $b->name} @$lanes ) {
					my $mapping_raw_bases = 0;
					my $duplication = 0;

        			my $name = $lane->name;
        			my $npg_qc = $lane->npg_qc_status;
					my $auto_qc_status = $lane->auto_qc_status();

                	my $lane_status_colour = get_colour_for_status($lane->qc_status);
                	my $lane_id = $lane->id;

					next if $npg_status_filt ne 'all' && $npg_status_filt ne $npg_qc;
					next if $auto_qc_status_filt ne 'all' && $auto_qc_status_filt ne $auto_qc_status;

					my $lane_mapstats = getMapStats($lane);
					if (%{$lane_mapstats}) {

						next if $bases_mapped_filt && $lane_mapstats->{mapping_raw_bases} < $bases_mapped_filt;
						next if $duplication_filt && $lane_mapstats->{duplication} < $duplication_filt;
						next if $rmdup_mapped_filt && $lane_mapstats->{rmdup_bases_mapped} < $rmdup_mapped_filt;
	
						print qq[<tr>];

						my $lane_id = $lane->id();
        				my $lane_qc_status = $lane->qc_status;
						my $lane_status_colour = get_colour_for_status($lane_qc_status);

						foreach my $status ($PASSED,$FAILED,$INVESTIGATE,$GT_PENDING,$PENDING) {
							my $state =  $lane_qc_status eq $status ? 'checked' : '';
							print $cgi->td("<input type='radio' name='$lane_id' value='$status' $state DISABLED>");
						}

						print qq[<td>$iname</td> ];
						print qq[
						<td style="background-color:$lane_status_colour;"><a href="$qc_grind_script?mode=$LIB_VIEW&amp;db=$database&amp;lib_id=$lid">$lname</a></td>
						];

						print qq [ <td style="background-color:$lane_status_colour;"><a href="$qc_grind_script?mode=$LANE_VIEW&amp;lane_id=$lane_id&amp;db=$database">$name</a></td> ];

						print $cgi->td([$npg_qc, $auto_qc_status, $lane_mapstats->{mapping_raw_bases}, $lane_mapstats->{duplication}, $lane_mapstats->{rmdup_bases_mapped}]);
						print qq[</tr>];
					}
				} # foreach lane

			} # foreach library
		} # sample
    }

    print qq[ </table> </fieldset> </div> ];

    print "<p align='center'>";
    print $cgi->submit(-name => 'download', -value  => 'Download');
    print "</p>";
    print $cgi->end_form;
}

sub downloadLaneData {
    my ($cgi, $vrtrack, $project) = @_;

	my ($npg_status_filt, $auto_qc_status_filt, $bases_mapped_filt, $duplication_filt, $rmdup_mapped_filt);
	$npg_status_filt = $cgi->param('npg_status');
	$auto_qc_status_filt = $cgi->param('auto_qc_status');
	$bases_mapped_filt = $cgi->param('bases_mapped');
	$duplication_filt = $cgi->param('duplication');
	$rmdup_mapped_filt = $cgi->param('rmdup_mapped');

    my $pname = $project->name;

	print $cgi->header(-type=>'text/tsv',  -attachment=>"$pname.tsv"); 
   # print $cgi->header(), $cgi->start_html(), "<pre>"; # testing
	#use Data::Dumper;print "<pre>parms:", Dumper($cgi->param()), "FILTER $npg_status_filt  $auto_qc_status_filt  $duplication_filt  $bases_mapped_filt   </pre>"; ## DEBUG

    my $samples = $project->samples();
    my %individuals2Samples;
    foreach my $sample ( @{ $samples } ) {
        my $ind = $sample->individual();
        my $indname = $ind->name();
        if( $individuals2Samples{ $indname } )
        {
            push( @{ $individuals2Samples{ $indname } }, $sample );
        }
        else{$individuals2Samples{ $indname } = [ $sample ];}
    }

    print "Sanger ID/Source Sample Name\tLibrary\tRun number\tNPG status\tQC status\tMapped\tDuplication\n";

    foreach( sort( keys( %individuals2Samples ) ) ) {
        my $iname = $_;
        my @samples = @{$individuals2Samples{ $_ }};

		foreach my $sample ( @samples ) {
			my @libraries = sort {$a->name cmp $b->name} @{$sample->libraries()};

# Uncomment if we want to show samples without libs
#			if (@libraries == 0) {
#				print "$iname\t\t\t\t\t\t\t\n";
#				next;
#			}

			foreach my $library ( @libraries ) {

				my $lname = $library->name;
				my $lid = $library->id();
				my $lanes = $library->lanes();

				foreach my $lane (sort {$a->name cmp $b->name} @$lanes ) {
					my $mapping_raw_bases = 0;
					my $duplication = 0;

        			my $name = $lane->name;
        			my $npg_qc = $lane->npg_qc_status;
					my $auto_qc_status = $lane->auto_qc_status();

					next if $npg_status_filt ne 'all' && $npg_status_filt ne $npg_qc;
					next if $auto_qc_status_filt ne 'all' && $auto_qc_status_filt ne $auto_qc_status;

					my $lane_mapstats = getMapStats($lane);
					if (%{$lane_mapstats}) {

						next if $bases_mapped_filt && $lane_mapstats->{mapping_raw_bases} < $bases_mapped_filt;
						next if $duplication_filt && $lane_mapstats->{duplication} < $duplication_filt;
						next if $rmdup_mapped_filt && $lane_mapstats->{rmdup_bases_mapped} < $rmdup_mapped_filt;

						print (join("\t",$iname,$lname,$name,$npg_qc,$auto_qc_status,$lane_mapstats->{mapping_raw_bases},$lane_mapstats->{duplication},$lane_mapstats->{rmdup_bases_mapped}),"\n");
					}
				} # foreach lane

			} # foreach lib
		} # foreach sample
	}
}

sub getMapStats {

	my $lane = shift;
	my %lane_mapstats;

	#work out which mapstats has the QC data (i.e. check for images in the mapstats)
	my @mappings = @{ $lane->mappings() };
	my $mapstats;
	foreach( sort {$a->row_id() <=> $b->row_id()} @mappings ) {
		my $map = $_;
		my $im = $map->images();
		if( @{$im} > 0 ){$mapstats = $map;}
	}
	if ($mapstats) {
		my $bases_mapped = $mapstats->bases_mapped;
		if ($bases_mapped) {   # sometimes the mapping fails
			$lane_mapstats{mapping_raw_bases} = sprintf("%.1f", ($mapstats->raw_bases()/1000000000)); # GB
			$lane_mapstats{duplication} = sprintf("%.1f", ($mapstats->rmdup_bases_mapped()/$mapstats->raw_bases)*100);
			$lane_mapstats{rmdup_bases_mapped} =  sprintf("%.1f", ($mapstats->rmdup_bases_mapped/1000000000)); # GB
		}
	}
	return \%lane_mapstats;
}

sub displayError
{
    my $message = $_[ 0 ];
    
    print qq[<h2>A problem occurred</h2>\n];
    print qq[<p class="error1">$message</p>\n];
    print $sw->footer();
    exit;
}

# returns CSS colour for a QC status
sub get_colour_for_status {
	my $status = shift;
	my $status_colour;

	if( $status eq $CLOSE_LIBRARY ) {
		$status_colour="#FFBF00";
	}
	elsif ($status eq $NO_QC) {
		$status_colour="#FFFFFF";
	}
	elsif ($status eq $PASSED) {
		$status_colour="#C0FFC0";
	}
	elsif ($status eq $FAILED) {
		$status_colour="#FFC0C0";
	}
	else {
		$status_colour="#F5F5F5";
	}

	return $status_colour;
}
