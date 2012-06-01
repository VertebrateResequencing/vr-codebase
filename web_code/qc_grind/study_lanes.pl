#!/usr/local/bin/perl -T
# 
# Displays Study Lane QC status information with update form
#
# Author:        cj5
# Created:       2012-03-28

BEGIN {
    $ENV{VRTRACK_HOST} = 'mcs4a';
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
use VertRes::QCGrind::Util;

my $pending_view = "../pending-view/pending_view.pl";

my $title = 'QC Breeze Lane Status Update';
my $sw  = SangerWeb->new({
    'title'   => $title,
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
});
my $utl = VertRes::QCGrind::Util->new();
my $main_script = $utl->{SCRIPTS}{DATABASES_VIEW};
my $proj_view_script = $utl->{SCRIPTS}{PROJECTS_VIEW};
my $lane_view_script = $utl->{SCRIPTS}{LANE_VIEW};

my $USER = $sw->username();
my $cgi = $sw->cgi();
my $form_submission = ($cgi->param('filter') or $cgi->param('update')) ? 1 : 0;

my $SCRIPT_NAME = $cgi->url(-relative=>1);

my $db = $cgi->param('db');
unless ($db) {
    print $sw->header();
	$utl->displayError( "No database ID",$sw );
}
my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $db, mode => 'rw');

my $projectID = $cgi->param('proj_id');
unless ($projectID) {
    print $sw->header();
	$utl->displayError( "No Project ID",$sw );
}
my $project = VRTrack::Project->new( $vrtrack, $projectID );

if ($cgi->param('download')) {
    &downloadLaneData($cgi, $vrtrack, $project);
	exit;
}

print $sw->header();
if ($cgi->param('update')) {
    &updateLaneData($cgi, $vrtrack, $project);
}

displayProjectLaneForm($cgi, $vrtrack, $db, $project);
print $sw->footer();
exit;

sub displayProjectLaneForm
{
    my ($cgi, $vrtrack, $database, $project) = @_;

	my ($gt_status_filt, $npg_status_filt, $auto_qc_status_filt, $bases_mapped_filt, $duplication_filt, $rmdup_mapped_filt);

	if ($form_submission) {
		$gt_status_filt = $cgi->param('gt_status');
		$npg_status_filt = $cgi->param('npg_status');
		$auto_qc_status_filt = $cgi->param('auto_qc_status');
		$bases_mapped_filt = $cgi->param('bases_mapped');
		$duplication_filt = $cgi->param('duplication');
		$rmdup_mapped_filt = $cgi->param('rmdup_mapped');
	}
	else {
		$gt_status_filt = 'all';
		$npg_status_filt = 'all';
		$auto_qc_status_filt = 'all';
	}
   
    my $samples = $project->samples();
    $utl->displayError( "Cant get samples for project: $projectID",$sw ) unless $samples;
    
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

	print $cgi->h2({-align=>"center", -style=>"font: normal 900 1.5em arial"},"<a href='$main_script'>QC Grind</a> Lane Status Update");
	print $cgi->h3({-align=>"center", -style=>"font: normal 700 1.5em arial"},"Database : <a href='$proj_view_script?db=$database'>$database</a>");
	print $cgi->h3({-align=>"center", -style=>"font: normal 700 1.5em arial"},"Project : $pname");
	print $cgi->h5({-style=>"font: arial"},"<a  target='_blank' href='$pending_view?mode=2&amp;db=$database'>Pending Requests</a>");
    print $cgi->br;

    print qq[
        <div class="centerFieldset">
        <fieldset>
        <legend>Lane status</legend>
        <table RULES=GROUPS>
        <tr>
		<th style="width: 40px">Pass</th>
		<th style="width: 40px">Fail</th>
		<th style="width: 40px">Invst</th>
		<th style="width: 40px">GTPend</th>
		<th style="width: 40px">Pend</th>
        <th>Individual</th>
        <th>Sample</th>
        <th>Library</th>
        <th>Run</th>
    	<th>Genotype</th>
        <th>NPG Status</th>
        <th>Auto QC Status</th>
        <th>Mapped</th>
        <th>Rmdup&#37;</th>
        <th>RmDup Mapped</th>
        <th></th>
        </tr>
    ];

    print $cgi->start_form;
    print $cgi->hidden("db","$database");
    print $cgi->hidden("proj_id","$projectID");

    print qq[ <tr style="background-color:#F5F5F5"> ];

	my @lane_status = qw (passed failed investigate gt_pending pending);
	foreach my $status (@lane_status) {
        print $cgi->th( $cgi->checkbox(-name=>$status, -checked=>1, -label=>''));
	}

	print $cgi->th(['','','','']);
	print "<th>", $cgi->popup_menu( -name => 'gt_status', -values => ['all','confirmed','unconfirmed','unchecked','unknown','candidate','wrong'], -default => 'all',), "</th>"; 
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

			my $sample_name = $sample->name;
			my @libraries = sort {$a->name cmp $b->name} @{$sample->libraries()};

			foreach my $library ( @libraries ) {
				my $libname = $library->name;
				my $lid = $library->id();
				my $lanes = $library->lanes();

				foreach my $lane (sort {$a->name cmp $b->name} @$lanes ) {

					if ($form_submission) { # filter on lane status check button
						next unless $cgi->param($lane->qc_status);
					}

					my $mapping_raw_bases = 0;
					my $duplication = 0;

        			my $lanename = $lane->name;
					my $gt_status = $lane->genotype_status;
        			my $npg_qc = $lane->npg_qc_status;
					my $auto_qc_status = $lane->auto_qc_status();

                	my $lane_status_colour = $utl->get_colour_for_status($lane->qc_status);
                	my $lane_id = $lane->id;

					next if $gt_status_filt ne 'all' && $gt_status_filt ne $gt_status;
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
						my $lane_status_colour = $utl->get_colour_for_status($lane_qc_status);

						my $disabled = $utl->{AUTH_USERS}{$USER} ? '' : 'DISABLED';
						foreach my $status ( qw (passed failed investigate gt_pending pending)) {
							my $state =  $lane_qc_status eq $status ? 'checked' : '';
							print $cgi->td("<input type='radio' name='$lane_id' value='$status' $state $disabled>");
						}

						my $gt_found = $lane_mapstats->{genotype_found};
						my $gt_ratio = $lane_mapstats->{genotype_ratio};
						my $gt_display = $gt_status;
						$gt_display .= " ($gt_found:$gt_ratio)" if $gt_found;
						my $gt_status_colour = get_gt_status_colour($gt_status,$sample_name,$gt_found);

						print $cgi->td($iname);
						print $cgi->td({-style=>"background-color:#F5F5F5"},$sample_name);
						print $cgi->td($libname);

						print qq [ <td style="background-color:$lane_status_colour;"><a href="$lane_view_script?lane_id=$lane_id&amp;db=$database">$lanename</a></td> ];
						print qq [ <td style="background-color:$gt_status_colour;">$gt_display</td> ];
						print $cgi->td([$npg_qc, $auto_qc_status, $lane_mapstats->{mapping_raw_bases}, $lane_mapstats->{duplication}, $lane_mapstats->{rmdup_bases_mapped}]);
						print qq[</tr>];
					}
				} # foreach lane
			} # foreach library
		} # sample
    }

    print qq[ <tr> ];
    if ($utl->{AUTH_USERS}{$USER}) {   
		print qq[ <td colspan=5 align='center'> ];
		print $cgi->submit(-name => 'update', -value  => 'Update');
		print qq[ </td> ];
    }

    print qq[ <td colspan=8 align='right'> ];
    print $cgi->submit(-name => 'download', -value  => 'Download');
    print qq[ </td> ];
    print qq[ </tr> ];
    print qq[ </table> </fieldset> </div> ];
    print $cgi->end_form;
}

sub updateLaneData
{
    my ($cgi, $vrtrack, $project) = @_;

    my @parameters = $cgi->param;
    foreach( @parameters ) {
        my $lane_id = $_;
        if( $lane_id =~ /^\d+$/ ) {
            my $status = $cgi->param($lane_id);
            
            #update the lane qc status
            my $lane= VRTrack::Lane->new($vrtrack, $lane_id);
            
            if( $lane->qc_status() ne $status ) {
                $lane->qc_status($status); #set the new status
                $lane->update;
            }
        }
    }
}

sub downloadLaneData {
    my ($cgi, $vrtrack, $project) = @_;

	my ($gt_status_filt, $npg_status_filt, $auto_qc_status_filt, $bases_mapped_filt, $duplication_filt, $rmdup_mapped_filt);

	$gt_status_filt = $cgi->param('gt_status');
	$npg_status_filt = $cgi->param('npg_status');
	$auto_qc_status_filt = $cgi->param('auto_qc_status');
	$bases_mapped_filt = $cgi->param('bases_mapped');
	$duplication_filt = $cgi->param('duplication');
	$rmdup_mapped_filt = $cgi->param('rmdup_mapped');

    my $pname = $project->name;

	print $cgi->header(-type=>'text/tsv',  -attachment=>"$pname.tsv"); 
   # print $cgi->header(), $cgi->start_html(), "<pre>"; # testing

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

    print join ("\t","Sanger ID","Sample Name","Library","Run number","Genotype","Lane QC status","NPG status","Auto QC status","Mapped","Duplication","RmDup Mapped"), "\n";

    foreach( sort( keys( %individuals2Samples ) ) ) {
        my $iname = $_;
        my @samples = @{$individuals2Samples{ $_ }};

		foreach my $sample ( @samples ) {
			my $sample_name = $sample->name;
			my @libraries = sort {$a->name cmp $b->name} @{$sample->libraries()};

			foreach my $library ( @libraries ) {

				my $libname = $library->name;
				my $lid = $library->id();
				my $lanes = $library->lanes();

				foreach my $lane (sort {$a->name cmp $b->name} @$lanes ) {

					next unless $cgi->param($lane->qc_status); # filter on lane status check button

					my $mapping_raw_bases = 0;
					my $duplication = 0;

        			my $lanename = $lane->name;
					my $gt_status = $lane->genotype_status;
					my $lane_qc_status = $lane->qc_status;
        			my $npg_qc = $lane->npg_qc_status;
					my $auto_qc_status = $lane->auto_qc_status();

					next if $gt_status_filt ne 'all' && $gt_status_filt ne $gt_status;
					next if $npg_status_filt ne 'all' && $npg_status_filt ne $npg_qc;
					next if $auto_qc_status_filt ne 'all' && $auto_qc_status_filt ne $auto_qc_status;

					my $lane_mapstats = getMapStats($lane);
					if (%{$lane_mapstats}) {

						my $gt_found = $lane_mapstats->{genotype_found};
						my $gt_ratio = $lane_mapstats->{genotype_ratio};
						my $gt_display = $gt_status;
						$gt_display .= " ($gt_found:$gt_ratio)" if $gt_found;

						next if $bases_mapped_filt && $lane_mapstats->{mapping_raw_bases} < $bases_mapped_filt;
						next if $duplication_filt && $lane_mapstats->{duplication} < $duplication_filt;
						next if $rmdup_mapped_filt && $lane_mapstats->{rmdup_bases_mapped} < $rmdup_mapped_filt;

						print (join("\t",$iname,$sample_name,$libname,$lanename,$gt_display,$lane_qc_status,$npg_qc,$auto_qc_status,$lane_mapstats->{mapping_raw_bases},$lane_mapstats->{duplication},$lane_mapstats->{rmdup_bases_mapped}),"\n");
					}
				} # foreach lane
			} # foreach lib
		} # foreach sample
	}
}
sub getMapStats {

	my $lane = shift;
	my %lane_mapstats;

	my @mappings = @{ $lane->mappings() };

	foreach my $mapstats ( @mappings ) {
        if ($mapstats->bases_mapped) {
            $lane_mapstats{mapping_raw_bases} = sprintf("%.1f", ($mapstats->bases_mapped()/1000000000)); # GB
            $lane_mapstats{duplication} = sprintf("%.1f", (1.0-($mapstats->rmdup_reads_mapped()/$mapstats->reads_mapped))*100);
            $lane_mapstats{rmdup_bases_mapped} =  sprintf("%.1f", ($mapstats->rmdup_bases_mapped/1000000000)); # GB
            $lane_mapstats{genotype_found} = $mapstats->genotype_found;
            $lane_mapstats{genotype_ratio} = sprintf("%.3f",$mapstats->genotype_ratio) if $mapstats->genotype_ratio;
            last;
        }
    }
    return \%lane_mapstats;
}

sub get_gt_status_colour
{
	my (($status,$iname,$gt_found)) = @_;

	my $green = "#CCFFCC";
	my $red = "#FFC0C0";
	my $status_colour;

	if ($status eq 'confirmed') {
		$status_colour=$green;
	}
	elsif ($status eq 'wrong') {
		$status_colour=$red;
	}
	elsif ($status eq 'unknown') {
		$status_colour=$red;
	}
	elsif ($gt_found) {
		if ($gt_found eq $iname) {
			$status_colour=$green;
		}
		else {
			$status_colour=$red;
		}
	}
	else {
		$status_colour="#F5F5F5";
	}
	return $status_colour;
}

