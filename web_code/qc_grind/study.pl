#!/usr/local/bin/perl -T
# 
# Displays Study QC status information at 'full' lane level with update form
#
# Author:        cj5
# Created:       2012-03-28

use strict;
use warnings;
use lib '/var/www/lib';

use CGI::Carp qw(fatalsToBrowser);
use URI;

use GD::Graph::bars;
use SangerPaths qw(core team145);
use SangerWeb;
use VRTrack::Project;
use VertRes::Utils::VRTrackFactory;
use VertRes::QCGrind::Util;
use Data::Dumper;

my $pending_view = "../pending-view/pending_view.pl";

my $title = 'QC Grind Study';
my $sw  = SangerWeb->new({
    'title'   => $title,
    'banner'  => q(),
    'jsfile'  => ['http://code.jquery.com/jquery-latest.js','/js/qc.js','/js/jquery.tablesorter.min.js'],
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
});
my $utl = VertRes::QCGrind::Util->new();
my $main_script = $utl->{SCRIPTS}{DATABASES_VIEW};
my $proj_view_script = $utl->{SCRIPTS}{PROJECTS_VIEW};
my $lane_view_script = $utl->{SCRIPTS}{LANE_VIEW};
my $study_lanes_script = $utl->{SCRIPTS}{STUDY_LANES};
my ($lane_name_filt, $sample_filt);

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
 
print $sw->header();

my $projectID = $cgi->param('proj_id');
unless ($projectID) {
	$utl->displayError( "No Project ID",$sw );
}
my $project = VRTrack::Project->new( $vrtrack, $projectID );

displayProjectLaneForm($cgi, $vrtrack, $db, $project);
print $sw->footer();
exit;

sub displayProjectLaneForm
{
    my ($cgi, $vrtrack, $database, $project) = @_;

	if ($form_submission) {
		$lane_name_filt = $cgi->param('lanename');
        $sample_filt = $cgi->param('sample');
    }
    else {
        $sample_filt = 'all';
    }
   
    my $samples = $project->samples();
    $utl->displayError( "Cant get samples for project: $projectID",$sw ) unless $samples;

    # Get lanelet and sample info
    my %lanes; # HoA, array of lanelets per lane hash
    my %lanes_samples; # Sample counts per lane
    my %sample_names;   # for sample filter menu

    foreach my $sample ( @{ $samples } ) {
        my $libraries = $sample->libraries();

        foreach my $library ( @$libraries ) {
            my $lanelets = $library->lanes();
            foreach my $lanelet ( @$lanelets ) {

                next if $lanelet->is_withdrawn;

                $sample_names{$sample->name}++;
                next if $sample_filt ne 'all' && $sample_filt ne $sample->name;

                my $lane_nm = $lanelet->name;
                $lane_nm =~ s/#.*//;
                push @{$lanes{$lane_nm}}, $lanelet;
                $lanes_samples{$lane_nm}++;
            }
        }
    }

    if ($cgi->param('update')) {
        &update_lanelet_status($cgi, $vrtrack, $project, %lanes);
    }

    my $pname = $project->name;

	print $cgi->h2({-align=>"center", -style=>"font: normal 900 1.5em arial"},"<a href='$main_script'>QC Grind</a> Project Lanes Update");
	print $cgi->h3({-align=>"center", -style=>"font: normal 700 1.5em arial"},"Database : <a href='$proj_view_script?db=$database'>$database</a>");
	print $cgi->h3({-align=>"center", -style=>"font: normal 700 1.5em arial"},"Project : $pname");
	print $cgi->h5({-style=>"font: arial"},"<a  target='_blank' href='$pending_view?mode=2&amp;db=$database'>Pending Requests</a>");
    print $cgi->br;

    print $cgi->start_form({name=>'qc_status'});
    print $cgi->hidden("db","$database");
    print $cgi->hidden("proj_id","$projectID");

    print $cgi->br;

    print qq[
        <div class="centerFieldset">
        <fieldset>
        <legend>Lane status</legend>
        <table RULES=GROUPS cellpadding="4" class="sortable">
        <thead> 
        <tr>
		<th style="width: 40px">Pass</th>
		<th style="width: 40px">Fail</th>
		<th style="width: 40px">Invst</th>
		<th style="width: 40px">GTPend</th>
		<th style="width: 40px">Pend</th>
        <th>Lane</th>
        <th colspan=2 title = "Click thumbnails to enlarge graphs">Lanelets</th>
        <th>Samples</th>
    	<th>GT confirmed</th>
        <th>NPG Passed</th>
        <th>Auto QC Passed</th>
        <th title="Lane raw bases in Gb">Raw Tot</th>
        <th title="Mean raw bases in Gb">Mean</th>
        <th title="Max raw bases in Gb">Max</th>
        <th title="Min raw bases in Gb">Min</th>
        <th title="Mapstats mapped bases in Gb">Mapped Tot</th>
        <th></th>
        </tr>
    ];

    # filters 
    print qq[ <tr style="background-color:#F5F5F5"> ];

	my @lane_status = qw (passed failed investigate gt_pending pending);
	foreach my $status (@lane_status) {
        print $cgi->td( $cgi->checkbox(-name=>$status, -checked=>1, -label=>'', -title=>"Filter on all lanelets $status"));
	}
	print "<td>", $cgi->textfield(-name=>'lanename', -size=>6), "</td>";
	print $cgi->td(['','']);
    my @sn = sort keys(%sample_names);
	print "<td>", $cgi->popup_menu( -name => 'sample', -values => ['all',@sn], -default => 'all', ), "</td>"; 
	print $cgi->td(['','','','','','','','']);
    print "<td>", $cgi->submit(-name => 'filter', -value  => 'Filter'), "</td>";
    print qq[ </tr> ];
    print qq[ </thead> ];
    print qq[<tbody>];

    foreach( sort( keys( %lanes) ) ) {

        my $lane_name = $_;
        next if $lane_name_filt && $lane_name !~ /^$lane_name_filt/;

        my $lanelets = $lanes{$lane_name};

        my %lane_data;

        foreach my $lanelet (sort {$a->name cmp $b->name} @$lanelets ) {

            $lane_data{num_lanelets}++;
            $lane_data{raw_bases} += sprintf("%.2f", ($lanelet->raw_bases/1000000000));
            
            $lane_data{raw_max} = $lanelet->raw_bases if (!exists $lane_data{raw_max} || $lanelet->raw_bases > $lane_data{raw_max});
            $lane_data{raw_min} = $lanelet->raw_bases if (!exists $lane_data{raw_min} || $lanelet->raw_bases < $lane_data{raw_max});

            my $lqc_status = $lanelet->qc_status;
            $lane_data{$lqc_status}++;

            my $gt_status = $lanelet->genotype_status;
            $lane_data{gt_confirmed}++ if $gt_status eq 'confirmed';

            my $npg_qc = $lanelet->npg_qc_status;
            $lane_data{npg_passed}++ if $npg_qc eq 'pass';

            my $auto_qc_status = $lanelet->auto_qc_status();
            $lane_data{auto_qc_passed}++ if $auto_qc_status eq 'passed';

            my @mappings = @{ $lanelet->mappings() };
            foreach my $mapstats ( @mappings ) {
                if ($mapstats->bases_mapped) {
                    $lane_data{bases_mapped} += sprintf("%.2f", ($mapstats->bases_mapped()/1000000000)); # GB
                    last;
                }
            }
        }
        ## filter out lane from results if Lanes status unchecked and all lanelets are of this status?
        if ($form_submission) { # filter on lane status checkboxes
            my $do_filter=0;
            foreach my $status (@lane_status) {
                next if $cgi->param($status);
                $do_filter++ if $lane_data{$status} == $lane_data{num_lanelets};
            }
            next if $do_filter;
        }

        $lane_data{raw_mean} = sprintf("%.2f", ($lane_data{raw_bases} / $lane_data{num_lanelets}));
        $lane_data{raw_min} = sprintf("%.2f", ($lane_data{raw_min}/1000000000));
        $lane_data{raw_max} = sprintf("%.2f", ($lane_data{raw_max}/1000000000));

        # data rows
        print qq[<tr>];
        my $lane_status_colour = '';

        my $disabled = $utl->{AUTH_USERS}{$USER} ? '' : 'DISABLED';
        foreach my $status (@lane_status) {
            print $cgi->td("<input type='radio' id='$lane_name' name='$lane_name' value='$status' $disabled>$lane_data{$status}");
            $lane_status_colour = $utl->get_colour_for_status($status) if $lane_data{$status} == $lane_data{num_lanelets};
        }

        print qq[ <td style="background-color:$lane_status_colour;"><a href="$study_lanes_script?;proj_id=$projectID&amp;db=$database&amp;runname=$lane_name" target="_blank">$lane_name</a></td> ];

        print $cgi->td({-align => "center"},$lane_data{num_lanelets});

        # generate a sizeable lanelet histogram
        my $img = lanelet_histogram($lanelets);
        print qq[<td style="padding:0px">]; 

        my $uri = URI->new("data:");
        $uri->media_type("image/png");
        $uri->data($img);

        print qq[
            <img src="$uri" width="30" height="30" title = "Click to enlarge/minimise" style="border:0px dotted #83A4C3;" onclick="if (this.width==30) {this.width=500;this.height=500;} else {this.width=30;this.height=30;}">
        ];
        print qq[</td>];
        print $cgi->td({-align => "center"},$lanes_samples{$lane_name});
        
        foreach my $col (qw(gt_confirmed npg_passed auto_qc_passed )) {
            print $cgi->td({-align => "center", -style=>"background-color:#F5F5F5"},$lane_data{$col});
        }
        foreach my $col (qw(raw_bases raw_mean raw_max raw_min bases_mapped )) {
            print $cgi->td({-align => "right"},$lane_data{$col});
        }

        print qq[</tr>];

    } # lane
    print qq[</tbody>];

    if ($utl->{AUTH_USERS}{$USER}) {   

        print qq[ <tfoot><tr style="background-color:#F5F5F5"> ];

        # toggle lane status radios
        foreach my $status (@lane_status) {
            print $cgi->td("<input type='radio' name='statusAll' title='Toggle $status' class='togglestatus' id='$status'/>");
        }

		print qq[ <td colspan=12 > ];
		print $cgi->submit(-name => 'update', -value  => 'Update');
        print $cgi->button(-id =>'resetRadios', -value=>'Reset');
		print qq[ </td> ];
        print qq[ </tr> </tfoot> ];
    }

    print qq[ </table> </fieldset> </div> ];
    print $cgi->end_form;
}

sub update_lanelet_status
{
    my ($cgi, $vrtrack, $project, %lanes) = @_;

    my @parameters = $cgi->param;
    foreach my $param ( @parameters ) {

        if( $param =~ /^\d+_\d$/ ) { # must be a lane_id radio button like 7260_1
            my $status = $cgi->param($param);
            
            foreach my $lanelet ( @{$lanes{$param}} ) {

                if( $lanelet->qc_status() ne $status ) {
                    $lanelet->qc_status($status);
                    $lanelet->update;
                }
                
                $utl->set_graph_lane_node_qc_status($lanelet->name, $status);
            }
        }
    }
}

sub lanelet_histogram 
{
    my ($lanelets) = @_;

    my (@names,@raw_bases);
    foreach my $lanelet (sort {$a->name cmp $b->name} @$lanelets ) {
            my ($nm) = $lanelet->name =~ /(#\S+)$/;
            push (@names,$nm);
            push (@raw_bases,sprintf("%.2f", ($lanelet->raw_bases/1000000000)));
    }
    my @data = (\@names,\@raw_bases);

    my $graph = new GD::Graph::bars(450,450);

    $graph->set( 
        x_label           => 'Lanelet',
        y_label           => 'Gb',
        title             => 'Lane raw bases',
    ) or die $graph->error;

    my $gd = $graph->plot(\@data) or die $graph->error;

    return $gd->png;
}

