#!/usr/local/bin/perl -T
# qc_grind.pl
# 
# Displays QC information for Sanger short-read sequencing to allow lanes
# to be passed/failed
#
# Author:        jws
# Maintainer:    jws
# Created:       2010-05-17

use strict;
use warnings;
use URI;

use lib '.';

#use SangerPaths qw(core team145);
use SangerPaths qw(core);
use lib '..';
use VRTrack::VRTrack;
use VRTrack::Project;
use VRTrack::Sample;
use VRTrack::Library;
use VRTrack::Lane;
use VertRes::Utils::VRTrackFactory;

use SangerWeb;

$ENV{PATH}= '/usr/local/bin';   # solely to stop taint from barfing

$|++;


#different modes/views possible
my $LANE_VIEW = 0;
my $LIB_VIEW = 1;
my $ERROR_DISPLAY = 2;
my $LIBS_VIEW = 3;
my $SPECIES_VIEW = 4;
my $LIB_UPDATE = 5;
my $PROJ_VIEW = 6;
my $SAMP_VIEW = 7;
my $LANE_UPDATE = 10;
my $LANES_UPDATE = 9;
my $FILTER_LANES = 11;

#possible states
my $PASSED = 'passed';
my $PENDING = 'pending';
my $NO_QC = 'no_qc';
my $FAILED = 'failed';
my $GT_PENDING = 'gt_pending';
my $INVESTIGATE = 'investigate';
my $CLOSE_LIBRARY = 'close';
my $OPEN_LIBRARY = 'open';

# Use SSO to authenticate, and then authorise against a list of people allowed
# to update the database.
my %AUTH_USERS = (  'jws' => 1, # Jim Stalker
                    'tk2' => 1, # Thomas Keane
                    'pd3' => 1, # Petr Danecek
                    'sb10'=> 1, # Sendu Bala
                    'rd'  => 1, # Richard Durbin
                    'kb1' => 1, # Karen McLaren
                    'ylx' => 1, #Yali Xue from Chris's group
                 );

###############################CSS Stuff#############################

my $css = <<CSS ;

.centerFieldset {
text-align:center;
}

.centerFieldset fieldset {
margin-left:auto;
margin-right:auto;
/* INHERITED ALIGNMENT IS CENTER. ONLY INCLUDE THIS IF YOU WANT */
/* TO CHANGE THE ALIGNMENT OF THE CONTENTS OF THE FIELDSET */
text-align:left;
}

.centerFieldset table {
    margin-left:auto;
    margin-right:auto;
}

table.summary {
    border-collapse: collapse;
    font: 0.9em Verdana, Arial, Helvetica, sans-serif;
}

table.summary td {
    white-space: nowrap;
    text-align: right;
    padding-right: 1em;
    padding-left: 1em;
    padding-top: 2px;
    padding-bottom: 2px;
    border-bottom: #83a4c3 dotted 1px;
}

table.summary tr.header th, table.summary tr.header td {
    font-weight:bold;
    border-bottom: #83a4c3 solid 1px;
    text-align: left;
    vertical-align:middle;
}

table.summary tr.level th, table.summary tr.level td {
    font-weight:bold;
    border-bottom: #83a4c3 dotted 1px;
    text-align: left;
    vertical-align:middle;
}

table.summary tr.total th, table.summary tr.total td {
    font-weight:bold;
    background-color: #CBDCED;
    border-bottom: #83a4c3 dotted 1px;
    text-align: right;
    vertical-align:middle;
}



input.btn {
  font: bold 150% 'trebuchet ms',helvetica,sans-serif;
  border:1px solid;
  border-color: #707070 #000 #000 #707070;
}

input.btnhov { 
  cursor:pointer;
  border-color: #c63 #930 #930 #c63; 
}

.clear
{
    clear: both;
    display: block;
    overflow: hidden;
    visibility: hidden;
    width: 0;
    height: 0;
}

img.preview:hover {
    width: 480px;
    height: 480px;
}

.thumbnail{
position: relative;
z-index: 0;
}

.thumbnail:hover{
background-color: transparent;
z-index: 50;
}

.thumbnail span{ /*CSS for enlarged image*/
position: absolute;
padding: 5px;
left: -1000px;
visibility: hidden;
text-decoration: none;
}

.thumbnail span img{ /*CSS for enlarged image*/
border-width: 0;
padding: 2px;
}

.thumbnail:hover span{ /*CSS for enlarged image on hover*/
visibility: visible;
top: 0;
left: -480px; /*position where enlarged image should offset horizontally */

}

CSS


###############################Entry Points#############################

my $sw  = SangerWeb->new({
    'title'   => q(QC Grind v3),
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
    'style'   => $css,
});

my $USER = $sw->username();
my $cgi = $sw->cgi();
#script name for self links
my $SCRIPT_NAME = $cgi->url(-relative=>1);


#decide on the entry point
my $mode = $cgi->param('mode');

if( defined( $mode ) && $mode == $ERROR_DISPLAY ) 
{
    my $message = $cgi->param('error_msg');
    
    print $sw->header();
    displayError( $message );
    print $sw->footer();
    
    exit;
}

# List available databases.
if( !defined( $cgi->param('db')) ) {
    print $sw->header();
    displayDatabasesPage();
    print $sw->footer();
    exit;
}

# All other entry points require a database
my $db = $cgi->param('db');
if( ! defined $db ) {
    redirectErrorScreen( $cgi, "Database must be defined!" );
    exit;
}

if( ! isDatabase( $db ) ) {
    redirectErrorScreen( $cgi, "Invalid database name!" );
    exit;
}

my $vrtrack = connectToDatabase( $db ); 
redirectErrorScreen( $cgi, "Failed to connect to database: $db" ) unless defined( $vrtrack );

if( ! defined( $mode ) && defined( $cgi->param('lane') ) ) {
    #when you just want to look at a lane by URL fiddling
    my $lane_name = $cgi->param('lane');
    if( $lane_name =~ /^\d+_\d+$/ ) {
        
        #get the row_id for the lane
        my $lane = VRTrack::Lane->new_by_hierarchy_name($vrtrack,$lane_name);
        
        print $sw->header();
        displayLane( $cgi, $vrtrack,$db, $lane->id(), undef );
        print $sw->footer();
        exit;
    }
    else {
        #error
        redirectErrorScreen( $cgi, "Invaid parameters for direct lane view (required are database and lane name)" );
        exit;
    }
}
elsif( $mode == $SPECIES_VIEW ) {
    print $sw->header();
    displayProjectsPage( $cgi, $vrtrack, $db );
    print $sw->footer();
    exit;
}
elsif( $mode == $PROJ_VIEW ) {
    my $pid = $cgi->param('proj_id');
    if( ! defined( $pid ) )
    {
        redirectErrorScreen( $cgi, "Must provide a project ID" );
        exit;
    }
    
    print $sw->header();
    displayProjectPage($cgi, $vrtrack, $db, $pid);
    print $sw->footer();
    exit;
}
elsif( $mode == $LIB_VIEW || $mode == $LIB_UPDATE ) {
    
    my $libID = $cgi->param('lib_id');
    
    if( $mode == $LIB_UPDATE )
    {
        my $state = $cgi->param("lib_update");
        updateLibrary($cgi, $vrtrack, $libID,$state);
    }
    
    print $sw->header();
    displayLibrary( $cgi, $vrtrack, $db,$libID);
    print $sw->footer();
    exit;
}
elsif( $mode == $LANE_VIEW || $mode == $LANE_UPDATE ) {
    
    my $laneID = $cgi->param('lane_id');
    
    if( $mode == $LANE_UPDATE )
    {
        my $state = $cgi->param("lane_update");
        updateLane($cgi,$vrtrack, $laneID,$state);
    }
    
    print $sw->header();
    displayLane( $cgi, $vrtrack, $db, $laneID);
    print $sw->footer();
    exit;
}
elsif( $mode == $LANES_UPDATE ) {
    
    updateLanes($cgi, $vrtrack, $db);

    # updateLanes can be called from a library view or a filter project view
    # return to the appropriate one.
    my $redir_url = $cgi->param('return_to');
    print $cgi->redirect( -URL => $redir_url, -method   => 'GET', -status   => 302 );
    exit;
}
elsif( $mode == $FILTER_LANES ) {
    my $projID = $cgi->param('proj_id');
    my $filter = $cgi->param('filter');
    unless ( defined $projID  && defined $filter )
    {
        redirectErrorScreen( $cgi, "Must provide a project ID & filter" );
        exit;
    }
    print $sw->header();
    displayFilteredLanes( $cgi, $vrtrack, $db, $projID,$filter);
    print $sw->footer();
    exit;
}
else
{
    redirectErrorScreen( $cgi, "Invalid mode!" );
}

########################################################################

sub updateLibrary {
    my ($cgi, $vrtrack, $libID,$state) = @_;

    #update the lane in the db
    my $library = VRTrack::Library->new( $vrtrack, $libID );
    if( $state ne $PASSED && $state ne $PENDING && $state ne $FAILED && $state ne $CLOSE_LIBRARY && $state ne $OPEN_LIBRARY )
    {
        redirectErrorScreen( $cgi, "Invalid library state found: $state" );
    }
    else
    {
        if( $state eq $CLOSE_LIBRARY || $state eq $OPEN_LIBRARY )
        {
            eval
            {
                if( $state eq $CLOSE_LIBRARY )
                {
                    $library->open( 0 );
                }
                else
                {
                    $library->open( 1 );
                }
                $library->update;
            };
        }
        else
        {
            eval
            {
                $library->qc_status( $state );
                $library->update;
            };
        }
        redirectErrorScreen( $cgi, "Failed to update library: $libID" ) unless ! $@;
    }
}


sub updateLane {
    my ($cgi, $vrtrack, $laneID, $state) = @_;
    
    #update the lane in the db
    my $lane = VRTrack::Lane->new( $vrtrack, $laneID );
    if( $state ne $PASSED && $state ne $PENDING && $state ne $FAILED && $state ne $GT_PENDING && $state ne $INVESTIGATE )
    {
        redirectErrorScreen( $cgi, "Invalid lane state found: $state" );
    }
    else
    {
        eval
        {
            $lane->qc_status( $state );
            $lane->update;
        };
        redirectErrorScreen( $cgi, "Failed to update lane: $laneID" ) unless ! $@;
    }

}


# multiple lane update from a library
sub updateLanes {
    my ($cgi, $vrtrack) = @_;
    my @parameters = $cgi->param;
    foreach( @parameters )
    {
        my $lane_id = $_;
        if( $lane_id =~ /^\d+$/ ) #assume its a lane ID
        {
            my $status = $cgi->param($lane_id);
            
            redirectErrorScreen( $cgi, "Invalid state for lane: $lane_id" ) unless $status eq $PASSED || $status eq $FAILED || $status eq $PENDING || $status eq $GT_PENDING || $status eq $INVESTIGATE;
            
            #update the lane qc status
            my $lane= VRTrack::Lane->new($vrtrack, $lane_id);
            
            if( $lane->qc_status() ne $status )
            {
                $lane->qc_status($status); #set the new status
                $lane->update;
            }
        }
    }
}


sub displayDatabasesPage {

    print qq[
        <h2 align="center" style="font: normal 900 1.5em arial">QC Grind</h2>
        <div class="centerFieldset">
        <fieldset style="width: 500px">
        <legend>Select dataset to QC</legend>
    ];
    
    my @dbs = VertRes::Utils::VRTrackFactory->databases();
    foreach( @dbs )
    {
        print qq[
            <p><a href="$SCRIPT_NAME?mode=$SPECIES_VIEW&amp;db=$_">].ucfirst( $_ ).qq[</a></p>
        ];
    }
    
    print qq[
        </fieldset>
        </div>
    ];


}


sub displayProjectsPage
{
    my ($cgi, $vrtrack,$database ) = @_;
    
    my @projects = sort {$a->name cmp $b->name} @{$vrtrack->projects()};
    
    print qq[
    <h2 align="center" style="font: normal 900 1.5em arial"><a href="$SCRIPT_NAME">QC Grind</a></h2>
    <h3 style="font: normal 700 1.5em arial">].ucfirst($database).qq[</h3>
    ];
    
    my $t = ucfirst( $database );
    print qq[
        <div class="centerFieldset">
        <fieldset style="width: 700px">
        <legend>Projects</legend>
        <table width="90%">
        <tr>
        <th>Name</th>
        <th>Accession</th>
        </tr>
    ];
    
    foreach( @projects )
    {
        my $project = $_;
        my $study = $project->study();
        my $acc = '-';
        if( defined( $study ) )
        {
            $acc = $study->acc();
        }
        my $pid = $project->id();
        
        my $name = $project->name;
        print qq[
            <tr>
                <td><a href="$SCRIPT_NAME?db=$database&amp;mode=$PROJ_VIEW&amp;proj_id=$pid">$name</a></td>
                <td>$acc</td>
            </tr>
        ];
    }
    print qq[
        </tr>
        </table>
        </fieldset>
        </div>
    ];
}

sub displayFilteredLanes {
    my ($cgi, $vrtrack, $database, $projectID, $filter) = @_;
    
    my $project = VRTrack::Project->new( $vrtrack, $projectID );
    displayError( "Cant get project: $projectID" ) unless $project;
    
    my $samples = $project->samples();
    displayError( "Cant get samples for project: $projectID" ) unless $samples;
    
    my $pname = $project->name;
    print qq[
    <h2 align="center" style="font: normal 900 1.5em arial"><a href="$SCRIPT_NAME">QC Grind</a></h2>
    <h3 style="font: normal 700 1.5em arial">
    <a href="$SCRIPT_NAME?mode=$SPECIES_VIEW&amp;db=$database">].ucfirst($database).qq[</a> : 
    <a href="$SCRIPT_NAME?mode=$PROJ_VIEW&amp;proj_id=$projectID&amp;db=$database">$pname</a> :
    Filter $filter lanes</h3><br/>
        ];
    
    my @filtlanes;

    foreach( sort { $a->ssid() <=> $b->ssid() } @$samples){
        my @libraries = sort {$a->name cmp $b->name} @{$_->libraries()};
        foreach( @libraries ) {
            my $library = $_;
            my $lanes = $library->lanes();
            foreach( @$lanes ) {
                my $lane = $_;
                my $lanename = $lane->name;
                my $id = $lane->id;
                if ($filter && $filter ne 'all'){ 
                    next unless $lane->qc_status eq $filter;
                }
                push @filtlanes, $lane;
            }
        }
    }
    my $current_url = $cgi->url(-query=>1,-relative=>1);
    displayQCLaneList(\@filtlanes,$vrtrack,$database,$current_url);
}

sub displayProjectPage
{
    my ($cgi, $vrtrack, $database, $projectID) = @_;
    
    my $project = VRTrack::Project->new( $vrtrack, $projectID );
    displayError( "Cant get project: $projectID" ) unless $project;
    
    my $samples = $project->samples();
    displayError( "Cant get samples for project: $projectID" ) unless $samples;
    
    my $pname = $project->name;
    print qq[
    <h2 align="center" style="font: normal 900 1.5em arial"><a href="$SCRIPT_NAME">QC Grind</a></h2>
    <h3 style="font: normal 700 1.5em arial"><a href="$SCRIPT_NAME?mode=$SPECIES_VIEW&amp;db=$database">].
    ucfirst($database).qq{</a> : $pname </h3><br/>};
    
    print qq[
        <div class="centerFieldset">
        <fieldset style="width: 700px">
        <legend>$pname</legend>
        <table RULES=GROUPS width="100%">
        <tr>
        <th>Sample</th>
        <th>Sanger</th>
        <th>Library</th>
        <th>Lanes</th>
        <th>Passed</th>
        <th>Pass seq</th>
        <th>Depth</th>
        <th><a href="$SCRIPT_NAME?mode=$FILTER_LANES&amp;filter=$PENDING&amp;proj_id=$projectID&amp;db=$database">Pending</a></th>
        </tr>
    ];
    
    my $grandTotalLanes = 0;
    my $grandTotalPassed = 0;
    my $grandTotalDepth = 0;
    my $grandTotalNoQC = 0;
    my $grandTotalSamples = 0;
    foreach( sort { $a->ssid() <=> $b->ssid() } @$samples){
        my $sample = $_;
        my $sname = $sample->name;
        my $is_ours = 1;
        if( defined( $project->study() ) )
        {
                $is_ours = $sample->is_sanger_sample();
        }
        
        print qq[
            <tr>
        ];
        
        if( ! $is_ours && $database =~ '.*g1k.*' )
        {
            print qq[
                <td bgcolor="red">$sname</td>
                <td bgcolor="red">NO</td>
            ];
        }
        else
        {
            print qq[
                <td>$sname</td>
                <td></td>
            ];
        }
        
        my @libraries = sort {$a->name cmp $b->name} @{$sample->libraries()};
        my $firstL = 1;
        my $sampleLanes = 0;
        my $sampleDepth = 0;
        my $samplePassed = 0;
        my $samplePassseq = 0;
        my $sample_no_qcLanes = 0;
        foreach( @libraries )
        {
            my $library = $_;
            my $lname = $library->name;
            my $lid = $library->id();
            my $lanes = $library->lanes();
            my $passedLanes = 0;
            my $passedBases = 0;
            my $sampledBases = 0;
            my $total_passed_bases = 0;
            my $lib_no_qcLanes = 0;
            foreach( @$lanes )
            {
                my $lane = $_;
                $sampleLanes ++;
                if( $lane->qc_status() eq $PASSED )
                {
                    $passedLanes++;
                    my $mapping = $lane->latest_mapping();
                    if( $mapping )
                    {
                        $passedBases += $mapping->rmdup_bases_mapped;
                        $sampledBases += $mapping->raw_bases();
                        $total_passed_bases += $lane->raw_bases();
                        
                        $samplePassed ++;
                    }
                }
                elsif( $lane->qc_status() eq $PENDING )
                {
                    $lib_no_qcLanes ++;
                }
            }
            
            my $depth = 0;
            my $pass_seq = 0;
            if( $total_passed_bases > 0 )
            {
                $depth = ( ( $passedBases / $sampledBases ) * $total_passed_bases ) / 3000000000;
                $depth = sprintf("%.2f", $depth);
                $sampleDepth += $depth;
                $samplePassseq += $total_passed_bases;
                $pass_seq = bp_to_nearest_unit($total_passed_bases, 1);
            }
            
            my $colour = get_colour_for_status( $library->open() ? $library->qc_status() : $CLOSE_LIBRARY );
            my $numLanes = @$lanes;
            if( ! $firstL ){print qq[<tr><td></td><td></td>];$firstL=0;}
            print qq[
                <td style="background-color:$colour;"><a href="$SCRIPT_NAME?mode=$LIB_VIEW&amp;db=$database&amp;lib_id=$lid">$lname</a></td>
            ];
            
            print $numLanes > 0 ? qq[<td>$numLanes</td><td>$passedLanes</td>] : qq[<td></td><td></td>];
            print $total_passed_bases > 0 ? qq[<td>$pass_seq</td>] : qq[<td></td>];
            print $depth > 0 ? qq[<td>$depth].qq[x</td>] : qq[<td></td>];
            
            if( $lib_no_qcLanes > 0 )
            {
                print qq[
                    <td>$lib_no_qcLanes</td>
                    </tr>
                ];
            }
            else{print qq[<td></td></tr>];}
            $firstL = 0;
            $sample_no_qcLanes += $lib_no_qcLanes;
        }
        $samplePassseq = bp_to_nearest_unit($samplePassseq, 1);
        print qq[<tr><th></th><th></th><th></th><th>$sampleLanes</th><th>$samplePassed</th><th>$samplePassseq</th><th>$sampleDepth].qq[x</th>];
        
        if( $sample_no_qcLanes > 0 )
        {
            print qq[<th>$sample_no_qcLanes</th>];
        }
        else
        {
            print qq[<th></th>];
        }
        
        print qq[</tr>];
        $grandTotalLanes += $sampleLanes;
        $grandTotalPassed += $samplePassed;
        $grandTotalDepth += $sampleDepth;
        $grandTotalNoQC += $sample_no_qcLanes;
                $grandTotalSamples ++;
    }
    
    print qq[<tr><tfoot><th>$grandTotalSamples</th><th></th><th></th><th>$grandTotalLanes</th><th>$grandTotalPassed</th><th>$grandTotalDepth].qq[x];
    
    if( $grandTotalNoQC > 0 )
    {
        print qq[<th>$grandTotalNoQC</th>];
    }
    else
    {
        print qq[<th></th>];
    }
    
    print qq[
        </tfoot></tr>
        </table>
        </fieldset>
        </div>
    ];
}

sub displayLane
{
    my ( $cgi, $vrtrack, $database, $laneID) = @_;
    
    my $lane;
    if( $laneID == -1 )
    {
        #pick a random library
        my $lanes = $vrtrack->qc_filtered_lane_names();
        $lane = VRTrack::Lane->new_by_hierarchy_name($vrtrack,$lanes->[0]);
        
        displayError("Null lane returned: $lane->[0]" ) unless defined( $lane );
    }
    else
    {
        $lane = VRTrack::Lane->new( $vrtrack, $laneID );
    }

    unless ($lane){
        displayError("Can't retrieve lane $laneID\n");
    }
    
    my $name = $lane->name;
    my ($project,$sample,$library) = get_hierarchy_for_lane($vrtrack, $lane);
    
    my $pname = $project->name;
    my $pid = $project->id();
    my $sname = $sample->name;
    my $sid = $sample->id();
    my $libname = $library->name;
    my $lid = $library->id();
    $name =~ /(\d+)_(\d+)/;
    my $run = $1;
    print qq[
        <h2 align="center" style="font: normal 900 1.5em arial">QC Grind</h2>
        <h3 style="font: normal 700 1.5em arial">
        <a href="$SCRIPT_NAME?mode=$SPECIES_VIEW&amp;db=$database">].ucfirst($database).qq[</a> : 
        <a href="$SCRIPT_NAME?mode=$PROJ_VIEW&amp;proj_id=$pid&amp;db=$database">$pname</a> :
        <a href="$SCRIPT_NAME?mode=$PROJ_VIEW&amp;proj_id=$pid&amp;db=$database">$sname</a> :
        <a href="$SCRIPT_NAME?mode=$LIB_VIEW&amp;lib_id=$lid&amp;db=$database">$libname</a> :
        $name 
        ];
                
        if( $run !~ /SRR|ERR/ )
        {
                print qq{[<a href="http://intweb.sanger.ac.uk/perl/prodsoft/npg/npg/run/$run">npg</a>]};
        }
        
        print "</h3>";
    
    my $status = $lane->qc_status;
    my $auto_qc_status = $lane->auto_qc_status();
    
    #work out which mapstats has the QC data (i.e. check for images in the mapstats)
    my @mappings = @{ $lane->mappings() };
    my $mapstats;
    foreach( sort {$a->row_id() <=> $b->row_id()} @mappings )
    {
        my $map = $_;
        my $im = $map->images();
        if( @{$im} > 0 ){$mapstats = $map;}
    }
    
    my $npg_qc = $lane->npg_qc_status;
    my ($error_rate, $adapter_perc, $reads_mapped, $bases_mapped, $reads_paired, $rmdup_reads_mapped, $rmdup_bases_mapped, $clip_bases,$adapter_reads);
    my $cycles = $lane->read_len;
    my $raw_reads = $lane->raw_reads;
    my $raw_bases = $lane->raw_bases;
    my $reads = $raw_reads; # if no mapping
    my $bases = $raw_bases; # if no mapping
    my $gt_status = $lane->genotype_status;
    my $total_bases = $bases;
    
    my $status_colour = get_colour_for_status($status);
    my $printed = 0;
    if ($mapstats)
    {
        my $images = $mapstats->images;
        
        # over-ride raw_reads and raw_bases if we have a mapping
        # This handles sampling for mapping
        $reads = $mapstats->raw_reads;
        $bases = $mapstats->raw_bases;
        $clip_bases = $mapstats->clip_bases;    # bases after clipping
        $total_bases = defined $clip_bases ? $clip_bases : $bases;
        $reads_mapped = $mapstats->reads_mapped;
        $bases_mapped = $mapstats->bases_mapped;
        $reads_paired = $mapstats->reads_paired;
        $rmdup_reads_mapped = $mapstats->rmdup_reads_mapped;
        $rmdup_bases_mapped = $mapstats->rmdup_bases_mapped;
        $error_rate = sprintf("%.3f",$mapstats->error_rate);
        if (defined $mapstats->adapter_reads)
        {
            $adapter_reads = $mapstats->adapter_reads;
            $adapter_perc = sprintf("%.1f",($adapter_reads/$reads)*100);
        }
        
        # genotypes
        my $gt_found = $mapstats->genotype_found;
        my $gt_ratio = sprintf("%.3f",$mapstats->genotype_ratio);
        my $gt_display = "$gt_status ($gt_found:$gt_ratio)";

        if ($bases_mapped)
        {   # sometimes the mapping fails
            my $dupe_rate_perc       = sprintf("%.2f", (1-($rmdup_reads_mapped/$reads_mapped))*100); #($bases_mapped - $rmdup_bases_mapped)/$bases_mapped);
            my $reads_mapped_perc   = sprintf("%.1f", ($reads_mapped/$reads)*100);
            my $reads_paired_perc   = sprintf("%.1f", ($reads_paired/$reads)*100);
            my $bases_mapped_perc   = sprintf("%.1f", ($bases_mapped/$total_bases)*100);
            my $rmdup_reads_mapped_perc   = sprintf("%.1f", ($rmdup_reads_mapped/$reads)*100);
            my $rmdup_bases_mapped_perc   = sprintf("%.1f", ($rmdup_bases_mapped/$total_bases)*100);
            my $clip_bases_perc    = sprintf("%.1f", ($clip_bases/$bases)*100);
            my $error_rate_perc    = sprintf("%.2f", ($error_rate*100));
            
            $reads_mapped = commify($reads_mapped);
            $bases_mapped = commify($bases_mapped);
            $reads_paired = commify($reads_paired);
            $rmdup_reads_mapped = commify($rmdup_reads_mapped);
            $rmdup_bases_mapped = commify($rmdup_bases_mapped);
            $reads = commify($reads);
            $bases = commify($bases);
            $raw_reads = commify($raw_reads);
            $raw_bases = bp_to_nearest_unit($raw_bases);
            $clip_bases = commify($clip_bases) if $clip_bases;
            
            print qq[
                <div class="centerFieldset">
                <fieldset style="width: 800px">
                <legend>QC Plots</legend>
                <table width="80%">
                <tr>
            ];
            
            #map the image names to the object
            my $imglist = $mapstats->images;
            my %images = map {$_->name => $_} @$imglist;
            my $img = $images{'gc-content.png'};
            if( $img )
            {
                print qq[<td>];printFullImage($img);print qq[</td>];
                delete( $images{'gc-content.png'} );
            }
            
            $img = $images{ "insert-size.png" };
            if( $img )
            {
                print qq[<td>];printFullImage($img);print qq[</td>];
                delete( $images{ "insert-size.png" } );
            }
            
            print qq[
                </tr>
                </table>
                <table>
                <tr>
            ];
            
            my @i = sort( { $a->name() cmp $b->name() } @{$images} );
            $images = \@i;
            
            foreach( keys( %images ) )
            {
                print qq[<td>];printPreviewImage($images{ $_ });print qq[</td>];
            }
            
            print qq[
                </tr>
                </table>
                </fieldset>
                </div>
            ];
            
            
            # Mapping stats
            my ($adapter_str, $rmapped_str, $rpaired_str, $rmdup_rmapped_str);
            my ($clip_str, $bmapped_str, $rmdup_bmapped_str);
            
            $adapter_str = "$adapter_reads ($adapter_perc\%)" if $adapter_reads;
            $rmapped_str = "$reads_mapped ($reads_mapped_perc\%)" if $reads_mapped;
            $rpaired_str = "$reads_paired ($reads_paired_perc\%)" if $reads_paired;
            $rmdup_rmapped_str = "$rmdup_reads_mapped ($rmdup_reads_mapped_perc\%)" if $rmdup_reads_mapped;

            if (defined $clip_bases)
            {
                $clip_str = "$clip_bases ($clip_bases_perc\% raw)";
            }
            else 
            {  
                $clip_str = "unclipped";
            }
            $bmapped_str = "$bases_mapped ($bases_mapped_perc\%)" if $bases_mapped;
            $rmdup_bmapped_str = "$rmdup_bases_mapped ($rmdup_bases_mapped_perc\%)" if $rmdup_bases_mapped;
            
            print qq[
                <br/>
                <div class="centerFieldset">
                <fieldset style="width: 900px">
                <legend>Mapping data</legend>
                <table width="100%">
                <tr>
                <td><table>
                <tr><td>Total Reads: </td><td>$raw_reads</td></tr>
                <tr><td>QC Reads: </td><td>$reads</td></tr>
                <tr><td>Reads w/adapter</td><td>$adapter_str</td></tr>
                <tr><td>Reads mapped: </td><td>$rmapped_str</td></tr>
                <tr><td>Reads paired: </td><td>$rpaired_str</td></tr>
                <tr><td>Reads mapped (rmdup): </td><td>$rmdup_rmapped_str</td></tr>
                </table>
                </td>
                <td><table>
                <tr><td>Total Bases: </td><td>$raw_bases</td></tr>
                <tr><td>QC Bases: </td><td>$bases</td></tr>
                <tr><td>Bases postclip: </td><td>$clip_str</td></tr>
                <tr><td>Bases mapped: </td><td>$bmapped_str</td></tr>
                <tr><td>&nbsp;</td><td>&nbsp;</td></tr>
                <tr><td>Bases mapped (rmdup): </td><td>$rmdup_bmapped_str</td></tr>
                </table>
                </td>
                <td><table>
                <tr><td>&nbsp;</td><td></td></tr>
                <tr><td>Cycles: </td><td>$cycles</td></tr>
                <tr><td>NPG QC: </td><td>$npg_qc</td></tr>
                <tr><td>Error rate: </td><td>$error_rate_perc%</td></tr>
                <tr><td>Genotype check: </td><td>$gt_display</td></tr>
                <tr><td>Duplication rate: </td><td>$dupe_rate_perc%</td></tr>
                </table>
                </td>
                </tr>
                </table>
                </fieldset>
                </div>
            ];
        }
            print qq[
        <div class="centerFieldset">
        <fieldset style="background-color: $status_colour;width: 500px">
        <legend>Lane QC : $status &nbsp;&nbsp;&nbsp;(Auto QC : $auto_qc_status)</legend>
    ];
    
    if ($AUTH_USERS{$USER})
    {
        my $laneID = $lane->id();
        print qq[
            <form action="$SCRIPT_NAME">
            <input type="hidden" name="mode" value="$LANE_UPDATE">
            <input type="hidden" name="db" value="$database">
            <input type="hidden" name="lane_id" value="$laneID">
        ];
        
        print qq[
            <p>
            <table width="400" align="center">
            <tr>
            <td align="center"><input type="submit" name="lane_update" value="$PENDING" class="btn" style="background-color: #F5F5F5;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
            <td align="center"><input type="submit" name="lane_update" value="$GT_PENDING" class="btn" style="background-color: #F5F5F5;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
            <td align="center"><input type="submit" name="lane_update" value="$INVESTIGATE" class="btn" style="background-color: #F5F5F5;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
            <td align="center"><input type="submit" name="lane_update" value="$FAILED" class="btn" style="background-color: #FFC0C0;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
            <td align="center"><input type="submit" name="lane_update" value="$PASSED" class="btn" style="background-color: #C0FFC0;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
            </tr>
            </table>
            </p>
            </form>
        ];
    }
    print qq[
        </fieldset>
        </div>
    ];
    }
}

sub displayLibrary
{
    my ( $cgi, $vrtrack, $database, $libID ) = @_;
    
    my $library;
    if( $libID == -1 )
    {
        #pick a random library
        my $libs = $vrtrack->qc_filtered_lib_hnames();
        $library = VRTrack::Library->new_by_hierarchy_name($vrtrack,$libs->[0]);
        
        displayError( "Can't get library ".$libs->[0] ) unless defined( $library );
    }
    else #jump into a particular library
    {
        $library = VRTrack::Library->new( $vrtrack, $libID );
        displayError( "Cant get library: $libID" ) unless $library;
    }
    
    my ($project,$sample) = get_hierarchy_for_library($vrtrack, $library);
    
    my $name = $library->name;
    my $project_name = $project->name;
    my $pid = $project->id();
    my $sample_name = $sample->name;
    my $sid = $sample->id();
    my $lib_ssid = $library->ssid();
    my $sample_ssid = $sample->ssid();
    print qq[
    <h2 align="center" style="font: normal 900 1.5em arial">QC Grind</h2>
    <h3 style="font: normal 700 1.5em arial">
    <a href="$SCRIPT_NAME?mode=$SPECIES_VIEW&amp;db=$database">].ucfirst($database).qq[</a> :
    <a href="$SCRIPT_NAME?mode=$PROJ_VIEW&amp;db=$database&amp;proj_id=$pid">$project_name</a> : 
    <a href="$SCRIPT_NAME?mode=$PROJ_VIEW&amp;db=$database&amp;proj_id=$pid">$sample_name</a> : 
    $name
    ];
    
    print qq[
    <a href="http://intweb.sanger.ac.uk/perl/prodsoft/npg/npg/search?query=$name">[NPG]</a>
        ];
    print "</h3>";
    
    #print the pass fail buttons on top
    my $lib_status = $library->qc_status;
    my $lib_status_colour = get_colour_for_status($lib_status);
    my $libcount = scalar @{$library->lane_ids};
    
    print qq[<br/>
    <div class="centerFieldset">
    <fieldset style="background-color: $lib_status_colour;width: 500px">
    <legend>Library QC : $lib_status</legend>
    ];
    
    if ($AUTH_USERS{$USER})
    {
        my $libID = $library->id();
        my $opentoggle =  $library->open() == 1 ? $CLOSE_LIBRARY : $OPEN_LIBRARY;
        my $opentoggle_bg = $opentoggle == $CLOSE_LIBRARY ? "#FFC0C0" : "#C0FFC0";

        print qq[
            <form action="$SCRIPT_NAME">
            <input type="hidden" name="mode" value="$LIB_UPDATE">
            <input type="hidden" name="db" value="$database">
            <input type="hidden" name="lib_id" value="$libID">
            <p>
            <table width="400" align="center">
            <tr>
            <td align="center"><input type="submit" name="lib_update" value="$PENDING" class="btn" style="background-color: #F5F5F5;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
            <td align="center"><input type="submit" name="lib_update" value="$FAILED" class="btn" style="background-color: #FFC0C0;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
            <td align="center"><input type="submit" name="lib_update" value="$PASSED" class="btn" style="background-color: #C0FFC0;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
            <td align="center"><input type="submit" name="lib_update" value="$opentoggle" class="btn" style="background-color: $opentoggle_bg;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
            ];
            
            print qq[
            </tr>
            <tr>
                <td align="left"></td>
                <td colspan="2" align="center"><b>$libcount</b> lanes in library</td>
                <td align="right"></td>
            </tr>
            </table>
            </p>
            </form>
        ];
    }
    
    print qq[
    </fieldset>
    </div>
    ];
    
    #print some info on the library
    my $insert = $library->insert_size();
    print qq[
    <div class="centerFieldset">
        <table>
        <tr>
        <td>Insert</td>
        <td>$insert</td>
    ];

    my $typeID = $library->library_type_id();
    if( defined( $typeID ) && $typeID > 0 )
    {
        my $type = VRTrack::Library_type->new($vrtrack, $typeID);
        my $type_ = $type->name();
        
        print qq[
                <td>Type</td>
                <td>$type_</td>
        ];
    }
    
    print qq[
        </tr>
        </table>
        </div>
    ];
    
    #get all the lanes for the library
    my $lanes = $library->lanes();
    my $current_url = $cgi->url(-query=>1,-relative=>1);
    displayQCLaneList($lanes,$vrtrack,$database,$current_url);
}

    
sub displayQCLaneList {
    my ($lanes,$vrtrack,$database,$return_url) = @_;
    
    my $count = 0;
    #style="width: 1200px">
    print qq[
    <div class="centerFieldset">
    <fieldset > 
    <legend>Lane data</legend>
    <form action="$SCRIPT_NAME">
    <input type="hidden" name="mode" value="$LANES_UPDATE">
    <input type="hidden" name="db" value="$database">
    <input type="hidden" name="return_to" value="$return_url">
    <table width="100%">
    <tr>
    <th style="width: 40px">Pass</th>
    <th style="width: 40px">Fail</th>
    <th style="width: 40px">Invst</th>
    <th style="width: 40px">GTPend</th>
    <th style="width: 40px">Pend</th>
    <th>Name</th>
    <th>Auto QC</th>
    <th>Cycles</th>
    <th>Bases</th>
    <th>Post-clip</th>
    <th>Adapter</th>
    <th>Mapped</th>
    <th>Paired</th>
    <th>Rmdup</th>
    <th>Error</th>
    <th>Genotype</th>
    </tr>
    ];
    
    my $noQC = 0;
    #sort by the lane name
    foreach( sort( { $a->name() cmp $b->name() } @{$lanes} ))
    {
        my $lane = $_;
        my ($project,$sample,$library) = get_hierarchy_for_lane($vrtrack, $lane);
        
        my $status = $lane->qc_status;
        my $name = $lane->name;
        
        #work out which mapstats has the QC data (i.e. check for images in the mapstats)
        my @mappings = @{ $lane->mappings() };
        my $mapstats;
        foreach( sort {$a->row_id() <=> $b->row_id()} @mappings )
        {
            my $map = $_;
            my $im = $map->images();
            if( @{$im} > 0 ){$mapstats = $map;}
        }
        
        my $npg_qc = $lane->npg_qc_status;
        my $imglist = [];
        my ($error_rate, $adapter_perc, $reads_mapped, $bases_mapped, $reads_paired, $rmdup_reads_mapped, $rmdup_bases_mapped, $clip_bases,$adapter_reads);
        my $cycles = $lane->read_len;
        my $reads = $lane->raw_reads;
        my $bases = $lane->raw_bases;
        my $gt_status = $lane->genotype_status;
        my $auto_qc_status = $lane->auto_qc_status();
        my $total_bases = $bases;
        
        my $printed = 0;
        if ($mapstats)
        {
            my $images = $mapstats->images;
            
            # over-ride raw_reads and raw_bases if we have a mapping
            # This handles sampling for mapping
            $reads = $mapstats->raw_reads;
            $bases = $mapstats->raw_bases;
            $clip_bases = $mapstats->clip_bases;    # bases after clipping
            $total_bases = defined $clip_bases ? $clip_bases : $bases;
            $reads_mapped = $mapstats->reads_mapped;
            $bases_mapped = $mapstats->bases_mapped;
            $reads_paired = $mapstats->reads_paired;
            $rmdup_reads_mapped = $mapstats->rmdup_reads_mapped;
            $rmdup_bases_mapped = $mapstats->rmdup_bases_mapped;
            $error_rate = sprintf("%.3f",$mapstats->error_rate);
            if (defined $mapstats->adapter_reads)
            {
                $adapter_reads = $mapstats->adapter_reads;
                $adapter_perc = sprintf("%.1f",($adapter_reads/$reads)*100);
            }
            
            # genotypes
            my $gt_found = $mapstats->genotype_found;
            my $gt_ratio = sprintf("%.3f",$mapstats->genotype_ratio);
            my $gt_display = "$gt_status ($gt_found:$gt_ratio)";

            if ($bases_mapped)
            {   # sometimes the mapping fails
                my $dupe_rate       = sprintf("%.4f", (1-$rmdup_reads_mapped/$reads_mapped)); #($bases_mapped - $rmdup_bases_mapped)/$bases_mapped);
                my $reads_mapped_perc   = sprintf("%.1f", ($reads_mapped/$reads)*100);
                my $reads_paired_perc   = sprintf("%.1f", ($reads_paired/$reads)*100);
                my $bases_mapped_perc   = sprintf("%.1f", ($bases_mapped/$total_bases)*100);
                my $rmdup_reads_mapped_perc   = sprintf("%.1f", ($rmdup_reads_mapped/$reads)*100);
                my $rmdup_bases_mapped_perc   = sprintf("%.1f", ($rmdup_bases_mapped/$total_bases)*100);
                my $clip_bases_perc    = sprintf("%.1f", ($clip_bases/$bases)*100);
                my $error_rate_perc     = sprintf("%.2f",($error_rate*100));
                
                $reads_mapped = commify($reads_mapped);
                $bases_mapped = commify($bases_mapped);
                $reads_paired = commify($reads_paired);
                $rmdup_reads_mapped = commify($rmdup_reads_mapped);
                $rmdup_bases_mapped = commify($rmdup_bases_mapped);
                $reads = commify($reads);
                $bases = commify($bases);
                $clip_bases = commify($clip_bases) if $clip_bases;
                
                my $id = $lane->id();
                my $lane_status_colour = get_colour_for_status($status);
                my $lane_auto_status_colour = defined($auto_qc_status) ? get_colour_for_status($auto_qc_status) : '';
                #print the table row
                print qq[
                <tr>
                <td><input type="radio" name="$id" value="$PASSED" ];
                print $status eq $PASSED ? 'checked' : '';
                print qq[></td>

                <td><input type="radio" name="$id" value="$FAILED" ];
                print $status eq $FAILED ? 'checked' : '';
                print qq[></td>


                <td><input type="radio" name="$id" value="$INVESTIGATE" ];
                print $status eq $INVESTIGATE ? 'checked' : '';
                print qq[></td>

                <td><input type="radio" name="$id" value="$GT_PENDING" ];
                print $status eq $GT_PENDING ? 'checked' : '';
                print qq[></td>

                <td><input type="radio" name="$id" value="$PENDING" ];
                print $status eq $PENDING ? 'checked' : '';
                print qq[></td>

                <td style="background-color:$lane_status_colour;"><a href="$SCRIPT_NAME?mode=$LANE_VIEW&amp;lane_id=$id&amp;db=$database">$name</a></td>
                ];
                
                print defined($auto_qc_status) ? qq[<td style="background-color:$lane_auto_status_colour;"></td>] : qq[<td>undef</td>];

                print qq[
                <td>$cycles</td>
                <td>$bases</td>
                <td>$clip_bases_perc%</td>
                <td>$adapter_perc%</td>
                <td>$reads_mapped_perc%</td>
                <td>$reads_paired_perc%</td>
                <td>$rmdup_reads_mapped_perc%</td>
                <td>$error_rate_perc%</td>
                <td>$gt_status ($gt_ratio)</td>
                ];
                
                #print the images
                foreach( @$images )
                {
                    my $image = $_;
                    my $iname = $image->name;
                    if( $iname eq "gc-content.png" || $iname eq "insert-size.png" || $iname eq "gc-depth.png" )
                    {
                        print qq[<td>];
                        printPreviewImage($image);
                        print qq[</td>];
                    }
                }
                print "</tr>";
                $printed = 1;
            }
        }
        else #no mapstats so just print the lane name etc.
        {
            print qq[
                <tr>
                <td></td><td></td><td></td>
                <td>$name</td>
                <td>$cycles</td>
                </tr>
            ];
        }
        
        if( ! $printed )
        {
            $noQC ++;
        }
    }
    print qq[
        </table>
    ];
    
    if ($AUTH_USERS{$USER})
    {   
        print qq[<input align="center" type="submit" name="lanes_update" value="Submit" class="btn" style="background-color: #F5F5F5;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" />];
    }
    
    print qq[
        </form>
        </fieldset>
        </div>
    ];

}


sub displayError
{
    my $message = $_[ 0 ];
    
    print qq[<h2>A problem occurred</h2>\n];
    print qq[<p class="error1">$message</p>\n];
    print $sw->footer();
    exit;
}

sub connectToDatabase
{
    my $database = $_[ 0 ];
    my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $db,
                                                          mode => 'rw');
    return $vrtrack;
}

sub redirectErrorScreen
{
    my $cgi = $_[ 0 ];
    my $error = 'An unknown error occurred';
    if( @_ == 2 )
    {
        $error = $_[ 1 ];
    }
    
    my $location = "$SCRIPT_NAME?mode=$ERROR_DISPLAY&amp;error_msg=$error";
    #print $cgi->header(); 
    #print "window.location=\"$location\";\n\n";
    print $cgi->redirect( -URL => $location, -method   => 'GET', -status   => 302 );
    #print "Location: $location";
    exit;
}

# returns CSS colour for a QC status
sub get_colour_for_status {
    my $status = shift;
    my $status_colour;

        if( $status eq $CLOSE_LIBRARY )
        {
                $status_colour="#FFBF00";
        }
        else
        {
                if ($status eq $NO_QC)
                {
                        $status_colour="#FFFFFF";
                }
                if ($status eq $PASSED)
                {
                        $status_colour="#C0FFC0";
                }
                elsif ($status eq $FAILED)
                {
                        $status_colour="#FFC0C0";
                }
                else 
                {
                        $status_colour="#F5F5F5";
                }
    }

    return $status_colour;
}

# add commas in to big numbers
sub commify { 
    my $number = shift;
    $number =  reverse $number;
    $number =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g; 
    return scalar reverse $number; 
}

# returns the project, sample, library objects for a given lane
sub get_hierarchy_for_lane 
{
    my $vrtrack = shift;
    my $lane = shift;
    my ($proj, $samp, $lib);
    eval {
        $lib = VRTrack::Library->new($vrtrack,$lane->library_id);
        if ($lib){
            $samp = VRTrack::Sample->new($vrtrack,$lib->sample_id);
            if ($samp){
                $proj = VRTrack::Project->new($vrtrack,$samp->project_id);
            }
        }
    };

    unless ($proj && $samp && $lib)
    {
        print_and_exit(sprintf('Cannot retrieve hierarchy for lane %s', $lane->hierarchy_name));
    }
    return ($proj, $samp, $lib);
}

# returns the project, sample, library objects for a given lane
sub get_hierarchy_for_library
{
    my $vrtrack = shift;
    my $library = shift;
    my ($proj, $samp);
    eval {
        $samp = VRTrack::Sample->new($vrtrack,$library->sample_id);
        if ($samp)
        {
            $proj = VRTrack::Project->new($vrtrack,$samp->project_id);
        }
    };

    unless ($proj && $samp)
    {
        print_and_exit(sprintf('Cannot retrieve hierarchy for lane %s', $library->hierarchy_name));
    }
    return ($proj, $samp);
}

sub printFullImage
{
    my $imageObj = shift;
    
    $imageObj->name =~ /\.(\w+)$/;
    my $type = $1;
    my $uri = URI->new("data:");
    $uri->media_type("image/$type");
    $uri->data($imageObj->image);
    my $caption = $imageObj->caption;
    
    print qq[
        <img src="$uri" alt="$caption" width="400" height="400" style="border:1px dotted #83A4C3;">
    ];
}

sub printPreviewImage
{
    my $imageObj = shift;
    
    $imageObj->name =~ /\.(\w+)$/;
    my $type = $1;
    my $uri = URI->new("data:");
    $uri->media_type("image/$type");
    $uri->data($imageObj->image);
    my $caption = $imageObj->caption;
    
    print qq[
        <a class="thumbnail" href="#">
        <img src="$uri" width="100" height="100" alt="$caption" style="border:1px dotted #83A4C3;">
        <span><img src="$uri" alt="$caption" style="border:1px dotted #83A4C3;"/></span></a>
    ];
}

sub isDatabase
{
        my $db = shift;
        
        my @dbs = VertRes::Utils::VRTrackFactory->databases();
        foreach( @dbs ){if( $db eq $_ ){return 1;}}
        return 0;
}


sub bp_to_nearest_unit {
    my ($bp,$dp) = @_;
    $dp = 2 unless defined $dp;

    my @units = qw( bp Kb Mb Gb Tb );

    my $power_ranger = int( ( length( abs($bp) ) - 1 ) / 3 );
    my $unit = $units[$power_ranger];
    my $unit_str;

    my $value = int( $bp / ( 10 ** ( $power_ranger * 3 ) ) );

    if ( $unit ne "bp" ){
        $unit_str = sprintf( "%.${dp}f%s", $bp / ( 10 ** ( $power_ranger * 3 ) ), " $unit" );
    }
    else{
        $unit_str = "$value $unit";
    }
    return $unit_str;
}

