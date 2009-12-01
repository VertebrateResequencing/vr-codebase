#!/usr/local/bin/perl -T
# qcgrind_v2.pl
# 
# Displays QC information for Sanger short-read sequencing to allow lanes
# to be passed/failed
#
# Author:        tk2
# Maintainer:    tk2
# Created:       Tue Sep 29 14:22:41 BST 2009 @599 /Internet Time/

use strict;
use warnings;
no warnings 'uninitialized';
use File::Basename;
use URI;

use SangerPaths qw(core team145);
use VRTrack::VRTrack;
use VRTrack::Project;
use VRTrack::Sample;
use VRTrack::Library;
use VRTrack::Lane;
use SangerWeb;

$ENV{PATH}= '/usr/local/bin';   # solely to stop taint from barfing

$|++;

#script name for self links
my $SCRIPT_NAME = basename( $0 );

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
my $SELECT_SPECIES_VIEW = 8;
my $LANES_UPDATE = 9;

my %DB_FOR_SPECIES = ( 'mouse'  => 'mouse_reseq_track', 'g1k'   => 'g1k_track', 'g1kmeta'	=> 'g1k_meta');

#possible filters for libaries/lanes
my $PASSED_FILTER = 'passed';
my $PENDING_FILTER = 'pending';
my $NO_QC = 'no_qc';
my $FAILED_FILTER = 'failed';
my $CLOSE_LIBRARY = 'close';
my $OPEN_LIBRARY = 'open';

# Use SSO to authenticate, and then authorise against a list of people allowed
# to update the database.
my %AUTH_USERS = (  'jws' => 1, # Jim Stalker
                    'tk2' => 1, # Thomas Keane
                    'pd3' => 1, # Petr Danecek
                    'sb10'=> 1, # Sendu Bala
                    'rd'  => 1, # Richard Durbin
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
    'title'   => q(QC Grind v2),
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
    'style'   => $css,
});

my $USER = $sw->username();
my $cgi = $sw->cgi();

#decide on the entry point
my $mode = $cgi->param('mode');

if( $mode == $SELECT_SPECIES_VIEW || ! defined( $mode ) )
{
    print $sw->header();
    print qq[
        <h2 align="center" style="font: normal 900 1.5em arial">QC Grind</h2>
        <div class="centerFieldset">
        <fieldset style="width: 900px">
        <legend>Select dataset to QC</legend>
    ];
    
    foreach( keys( %DB_FOR_SPECIES ) )
    {
        print qq[
            <p><a href="$SCRIPT_NAME?mode=$SPECIES_VIEW&amp;sp=$_">].ucfirst($_).qq[</a></p>
        ];
    }
    
    print qq[
        </fieldset>
        </div>
    ];
    print $sw->footer();
    exit;
}
elsif( $mode == $SPECIES_VIEW )
{
    my $species = $cgi->param('sp');
    if( ! defined $species ) 
    {
        redirectErrorScreen( $cgi, "Species must be defined!" );
        exit;
    }
    
    my $database = $DB_FOR_SPECIES{ lc( $species ) };
    if( ! defined $database ) 
    {
        redirectErrorScreen( $cgi, "Cant find the database name for species: $species" );
        exit;
    }
    
    print $sw->header();
    displayProjectsPage( $cgi, $database, $species );
    print $sw->footer();
    exit;
}
elsif( $mode == $PROJ_VIEW )
{
    my $species = $cgi->param('sp');
    
    if( ! defined $species ) 
    {
        redirectErrorScreen( $cgi, "Species must be defined!" );
        exit;
    }
    
    my $pid = $cgi->param('proj_id');
    if( ! defined( $pid ) )
    {
        redirectErrorScreen( $cgi, "Must provide a project ID" );
        exit;
    }
    
    my $database = $DB_FOR_SPECIES{ lc( $species ) };
    
    if( ! defined $database ) 
    {
        redirectErrorScreen( $cgi, "Cant find the database name for species: $species" );
        exit;
    }
    
    print $sw->header();
    displayProjectPage($cgi, $database, $species, $pid);
    print $sw->footer();
    exit;
}
elsif( $mode == $LIB_VIEW || $mode == $LIB_UPDATE )
{
    my $species = $cgi->param('sp');
    
    if( ! defined $species ) 
    {
        redirectErrorScreen( $cgi, "Species must be defined!" );
        exit;
    }
    
    my $libID = $cgi->param('lib_id');
    my $filter = $cgi->param('filter');
    my $database = $DB_FOR_SPECIES{ lc( $species ) };
    
    if( ! defined $database ) 
    {
        redirectErrorScreen( $cgi, "Cant find the database name for species: $species" );
        exit;
    }
    
    if( $mode == $LIB_UPDATE )
    {
        my $state = $cgi->param("lib_update");
        
        #connect to the db
        my $vrtrack = connectToDatabase( $database );
        
        redirectErrorScreen( $cgi, "Failed to connect to database: $database" ) unless defined( $vrtrack );
        
        #update the lane in the db
        my $library = VRTrack::Library->new( $vrtrack, $libID );
        if( $state ne $PASSED_FILTER && $state ne $PENDING_FILTER && $state ne $FAILED_FILTER && $state ne $CLOSE_LIBRARY && $state ne $OPEN_LIBRARY )
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
    
    print $sw->header();
    displayLibrary( $cgi, $database, defined( $libID ) ? $libID : -1, defined( $filter ) ? $filter : undef, $species );
    print $sw->footer();
    exit;
}
elsif( $mode == $LANE_VIEW || $mode == $LANE_UPDATE )
{
    my $species = $cgi->param('sp');
    
    if( ! defined $species )
    {
        redirectErrorScreen( $cgi, "Species must be defined!" );
        exit;
    }
    
    my $laneID = $cgi->param('lane_id');
    my $filter = $cgi->param('filter');
    my $database = $DB_FOR_SPECIES{ lc( $species ) };
    
    if( ! defined $species ) 
    {
        redirectErrorScreen( $cgi, "Cant find the database name for species: $species" );
        exit;
    }
    
    if( $mode == $LANE_UPDATE )
    {
        my $state = $cgi->param("lane_update");
        
        #connect to the db
        my $vrtrack = connectToDatabase( $database );
        
        redirectErrorScreen( $cgi, "Failed to connect to database: $database" ) unless defined( $vrtrack );
        
        #update the lane in the db
        my $lane = VRTrack::Lane->new( $vrtrack, $laneID );
        if( $state ne $PASSED_FILTER && $state ne $PENDING_FILTER && $state ne $FAILED_FILTER )
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
    
    print $sw->header();
    displayLane( $cgi, $database, defined( $laneID ) ? $laneID : -1, defined( $filter ) ? $filter : undef, $species );
    print $sw->footer();
    exit;
}
elsif( $mode == $LANES_UPDATE )
{
    my $species = $cgi->param('sp');
    my @parameters = $cgi->param;
    my $libID = $cgi->param('lib_id');
    my $filter = $cgi->param('filter');
    
    my $database = $DB_FOR_SPECIES{ lc( $species ) };
    #connect to the db
    my $vrtrack = connectToDatabase( $database );
    
    redirectErrorScreen( $cgi, "Failed to connect to database: $database" ) unless defined( $vrtrack );
    
    foreach( @parameters )
    {
        my $lane_id = $_;
        if( $lane_id =~ /^\d+$/ ) #assume its a lane ID
        {
            my $status = $cgi->param($lane_id);
            
            redirectErrorScreen( $cgi, "Invalid state for lane: $lane_id" ) unless $status eq $PASSED_FILTER || $status eq $FAILED_FILTER || $status eq $PENDING_FILTER;
            
            #update the lane qc status
            my $lane= VRTrack::Lane->new($vrtrack, $lane_id);
            
            if( $lane->qc_status() ne $status )
            {
                $lane->qc_status($status); #set the new status
                $lane->update;
            }
        }
    }

    print $sw->header();
    displayLibrary( $cgi, $database, defined( $libID ) ? $libID : -1, defined( $filter ) ? $filter : undef, $species );
    print $sw->footer();
    exit;
}
elsif( $mode == $ERROR_DISPLAY )
{
    my $message = $cgi->param('error_msg');
    
    print $sw->header();
    displayError( $message );
    print $sw->footer();
    
    exit;
}
else
{
    redirectErrorScreen( $cgi, "Invalid mode!" );
}

###############################Common Functions#########################

sub displayProjectsPage
{
    my ($cgi, $database, $species) = @_;
    
    #connect to the db
    my $vrtrack = connectToDatabase( $database );
    
    redirectErrorScreen( $cgi, "Failed to connect to database: $database" ) unless defined( $vrtrack );
    
    my @projects = sort {$a->name cmp $b->name} @{$vrtrack->projects()};
    
    print qq[
        <h2 align="center" style="font: normal 900 1.5em arial">QC Grind</h2>
        <h3 style="font: normal 700 1.5em arial">].ucfirst($species).qq[</h3>
    ];
    
    my $t = ucfirst( $species );
    print qq[
        <div class="centerFieldset">
        <fieldset style="width: 800px">
        <legend>$species</legend>
        <table width="100%">
        <tr>
        <th>Project</th>
        <th>Accession</th>
        <th>Samples</th>
        <th>Libraries</th>
        <th>Passed</th>
        <th>Depth</th>
        <th>Pending</th>
        </tr>
    ];
    
    foreach( @projects )
    {
        my $project = $_;
        my $acc = $project->study()->acc();
        my $pid = $project->id();
        
        my $name = $project->name;
        print qq[
            <tr>
            <td><a href="$SCRIPT_NAME?sp=$species&amp;mode=$PROJ_VIEW&amp;proj_id=$pid">$name</a></td>
            <td>$acc</td>
        ];
        
        my $first = 1;
        my @samples = sort {$a->ssid cmp $b->ssid} @{$project->samples()};
        
        my $totalLibraries = 0;
        my $totalPassedLibraries = 0;
        my $totalDepth = 0;
        foreach( @samples )
        {
            my $sample = $_;
            my $sname = $sample->name;
            my $sid = $sample->id();
            my $libraries = $sample->libraries();
            my $libs = @$libraries;
            if( !$first ){print qq[<tr><td></td><td></td>];}
            print qq[
                <td>$sname</td>
                <td>$libs</td>
            ];
            
            my $passed = 0;
            my $depth = 0;
            my $pending = 0;
            foreach( @$libraries )
            {
                my $library = $_;
                if( $library->qc_status() eq $PASSED_FILTER )
                {
                    $passed ++;
                }
                
                my $lanes = $library->lanes();
                foreach( @$lanes )
                {
                    if( $_->qc_status() eq $PENDING_FILTER )
                    {
                        $pending ++;
                    }
                }
                
                $depth += $library->projected_passed_depth(3000000000);
                my $lids = $library->lane_ids();
            }
            
            $totalDepth += $depth;
            print qq[
                <td>$passed</td>
                <td>$depth].qq[x</td>
            ];
            
            print $pending > 0 ? qq[<td>$pending</td>] : qq[<td></td>];
            print qq[</tr>];
            
            $totalLibraries += $libs;
            $totalPassedLibraries += $passed;
            $first = 0;
        }
        print qq[
            <tr>
            <th></th><th></th><th>Total</th><th>$totalLibraries</th><th>$totalPassedLibraries</th>
        ];
        
        if( $totalDepth > 0 )
        {
            print qq[
                <th>$totalDepth].qq[x</th>
            ];
        }
        else{print qq[<th></th>];}
    }
    print qq[
        </tr>
        </table>
        </fieldset>
        </div>
    ];
}

sub displayProjectPage
{
    my ($cgi, $database, $species, $projectID) = @_;
    
    #connect to the db
    my $vrtrack = connectToDatabase( $database );
    
    redirectErrorScreen( $cgi, "Failed to connect to database: $database" ) unless defined( $vrtrack );
    
    my $project = VRTrack::Project->new( $vrtrack, $projectID );
    redirectErrorScreen( $cgi, "Cant get project: $projectID" ) unless $project;
    
    my $samples = $project->samples();
    redirectErrorScreen( $cgi, "Cant get samples for project: $projectID" ) unless $samples;
    
    my $pname = $project->name;
    my $ssid = $project->ssid();
	if( defined( $ssid ) && $ssid > 0 )
	{
		print qq[
    	    <h2 align="center" style="font: normal 900 1.5em arial">QC Grind</h2>
			<h3 style="font: normal 700 1.5em arial"><a href="$SCRIPT_NAME?mode=$SPECIES_VIEW&amp;sp=$species">].ucfirst($species).qq{</a> :
			$pname [<a href="http://psd-production.internal.sanger.ac.uk:6600/projects/$ssid/workflows/1">SS</a>] 
			</h3><br/>
		};
	}
    
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
        <th>Depth</th>
        <th>Pending</th>
        </tr>
    ];
    
    my @samples_ = @$samples;
    @samples_ = sort { $a->ssid() <=> $b->ssid() } @samples_;
    my $grandTotalLanes = 0;
    my $grandTotalPassed = 0;
    my $grandTotalDepth = 0;
    my $grandTotalNoQC = 0;
    foreach( @samples_ )
    {
        my $sample = $_;
        my $sname = $sample->name;
        my $is_ours = $sample->is_sanger_sample();
        
        print qq[
            <tr>
        ];
        
        if( ! $is_ours && $species eq 'g1k' )
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
                if( $lane->qc_status() eq $PASSED_FILTER )
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
                elsif( $lane->qc_status() eq $PENDING_FILTER )
                {
                    $lib_no_qcLanes ++;
                }
            }
            
            my $depth = 0;
            if( $total_passed_bases > 0 )
            {
                $depth = ( ( $passedBases / $sampledBases ) * $total_passed_bases ) / 3000000000;
                $depth = sprintf("%.2f", $depth);
                $sampleDepth += $depth;
            }
            
			my $colour = get_colour_for_status( $library->open() ? $library->qc_status() : $CLOSE_LIBRARY );
            my $numLanes = @$lanes;
            if( ! $firstL ){print qq[<tr><td></td><td></td>];$firstL=0;}
            print qq[
                <td style="background-color:$colour;"><a href="$SCRIPT_NAME?mode=$LIB_VIEW&amp;sp=$species&amp;lib_id=$lid">$lname</a></td>
            ];
            
            print $numLanes > 0 ? qq[<td>$numLanes</td><td>$passedLanes</td>] : qq[<td></td><td></td>];
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
        print qq[<tr><th></th><th></th><th></th><th>$sampleLanes</th><th>$samplePassed<th>$sampleDepth].qq[x</th>];
        
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
    }
    
    print qq[<tr><tfoot><th>Grand Total</th><th></th><th></th><th>$grandTotalLanes</th><th>$grandTotalPassed</th><th>$grandTotalDepth].qq[x];
    
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
    my ( $cgi, $database, $laneID, $filter, $species ) = @_;
    
    #connect to the db
    my $vrtrack = connectToDatabase( $database );
    
    redirectErrorScreen( $cgi, "Failed to connect to database: $database" ) unless defined( $vrtrack );
    
    my @filter = defined( $filter ) ? ( $filter ) : ();
    my $lane;
    if( $laneID == -1 )
    {
        #pick a random library
        my $lanes = $vrtrack->qc_filtered_lane_names( @filter );
        $lane = VRTrack::Lane->new_by_hierarchy_name($vrtrack,$lanes->[0]);
        
        redirectErrorScreen( $cgi, "Null lane returned: $lane->[0]" ) unless defined( $lane );
    }
    else
    {
        $lane = VRTrack::Lane->new( $vrtrack, $laneID );
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
        <a href="$SCRIPT_NAME?mode=$SPECIES_VIEW&amp;sp=$species">].ucfirst($species).qq[</a> : 
        <a href="$SCRIPT_NAME?mode=$PROJ_VIEW&amp;proj_id=$pid&amp;sp=$species">$pname</a> :
        <a href="$SCRIPT_NAME?mode=$PROJ_VIEW&amp;proj_id=$pid&amp;sp=$species">$sname</a> :
        <a href="$SCRIPT_NAME?mode=$LIB_VIEW&amp;lib_id=$lid&amp;sp=$species">$libname</a> :
        $name [<a href="http://intweb.sanger.ac.uk/perl/prodsoft/npg/npg/run/$run">npg</a>]
        </h3>
    ];
    
    my $status = $lane->qc_status;
    my $mapstats = $lane->latest_mapping;
    my $npg_qc = $lane->npg_qc_status;
    my ($error_rate, $adapter_perc, $reads_mapped, $bases_mapped, $reads_paired, $rmdup_reads_mapped, $rmdup_bases_mapped, $clip_bases,$adapter_reads);
    my $cycles = $lane->read_len;
    my $reads = $lane->raw_reads;
    my $bases = $lane->raw_bases;
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
            my $dupe_rate       = sprintf("%.4f", (1-$rmdup_reads_mapped/$reads_mapped)); #($bases_mapped - $rmdup_bases_mapped)/$bases_mapped);
            my $reads_mapped_perc   = sprintf("%.1f", ($reads_mapped/$reads)*100);
            my $reads_paired_perc   = sprintf("%.1f", ($reads_paired/$reads)*100);
            my $bases_mapped_perc   = sprintf("%.1f", ($bases_mapped/$total_bases)*100);
            my $rmdup_reads_mapped_perc   = sprintf("%.1f", ($rmdup_reads_mapped/$reads)*100);
            my $rmdup_bases_mapped_perc   = sprintf("%.1f", ($rmdup_bases_mapped/$total_bases)*100);
            my $clip_bases_perc    = sprintf("%.1f", ($clip_bases/$bases)*100);
            
            $reads_mapped = commify($reads_mapped);
            $bases_mapped = commify($bases_mapped);
            $reads_paired = commify($reads_paired);
            $rmdup_reads_mapped = commify($rmdup_reads_mapped);
            $rmdup_bases_mapped = commify($rmdup_bases_mapped);
            $reads = commify($reads);
            $bases = commify($bases);
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
                <tr><td>Reads: </td><td>$reads</td></tr>
                <tr><td>Reads w/adapter</td><td>$adapter_str</td></tr>
                <tr><td>Reads mapped: </td><td>$rmapped_str</td></tr>
                <tr><td>Reads paired: </td><td>$rpaired_str</td></tr>
                <tr><td>Reads mapped (rmdup): </td><td>$rmdup_rmapped_str</td></tr>
                </table>
                </td>
                <td><table>
                <tr><td>Bases: </td><td>$bases</td></tr>
                <tr><td>Bases postclip: </td><td>$clip_str</td></tr>
                <tr><td>Bases mapped: </td><td>$bmapped_str</td></tr>
                <tr><td>&nbsp;</td><td>&nbsp;</td></tr>
                <tr><td>Bases mapped (rmdup): </td><td>$rmdup_bmapped_str</td></tr>
                </table>
                </td>
                <td><table>
                <tr><td>Cycles: </td><td>$cycles</td></tr>
                <tr><td>NPG QC: </td><td>$npg_qc</td></tr>
                <tr><td>Error rate: </td><td>$error_rate</td></tr>
                <tr><td>Genotype check: </td><td>$gt_display</td></tr>
                <tr><td>Duplication rate: </td><td>$dupe_rate</td></tr>
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
        <legend>Lane QC : $status</legend>
    ];
    
    if ($AUTH_USERS{$USER})
    {
        my $laneID = $lane->id();
        print qq[
            <form action="$SCRIPT_NAME">
            <input type="hidden" name="mode" value="$LANE_UPDATE">
            <input type="hidden" name="sp" value="$species">
            <input type="hidden" name="lane_id" value="$laneID">
        ];
        if( defined( $filter ) )
        {
            print qq[<input type="hidden" name="filter" value="$filter">];
        }
        
        print qq[
            <p>
            <table width="400" align="center">
            <tr>
            <td align="center"><input type="submit" name="lane_update" value="$PENDING_FILTER" class="btn" style="background-color: #F5F5F5;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
            <td align="center"><input type="submit" name="lane_update" value="$FAILED_FILTER" class="btn" style="background-color: #FFC0C0;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
            <td align="center"><input type="submit" name="lane_update" value="$PASSED_FILTER" class="btn" style="background-color: #C0FFC0;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
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
    my ( $cgi, $database, $libID, $filter, $species ) = @_;
    
    #connect to the db
    my $vrtrack = connectToDatabase( $database );
    
    redirectErrorScreen( $cgi, "Failed to connect to database: $database" ) unless defined( $vrtrack );
    
    my @filter = defined( $filter ) ? ( $filter ) : ();
    my $library;
    if( $libID == -1 )
    {
        #pick a random library
        my $libs = $vrtrack->qc_filtered_lib_names( @filter );
        $library = VRTrack::Library->new_by_hierarchy_name($vrtrack,$libs->[0]);
        
        redirectErrorScreen( $cgi, "Null library returned: $libs->[0]" ) unless defined( $library );
    }
    else #jump into a particular library
    {
        $library = VRTrack::Library->new( $vrtrack, $libID );
        redirectErrorScreen( $cgi, "Cant get library: $libID" ) unless $library;
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
        <a href="$SCRIPT_NAME?mode=$SPECIES_VIEW&amp;sp=$species">].ucfirst($species).qq[</a> :
        <a href="$SCRIPT_NAME?mode=$PROJ_VIEW&amp;sp=$species&amp;proj_id=$pid">$project_name</a> : 
        <a href="$SCRIPT_NAME?mode=$PROJ_VIEW&amp;sp=$species&amp;proj_id=$pid">$sample_name</a> : 
        $name
	];
	
	if( defined( $lib_ssid ) && $lib_ssid > 0 )
	{
		print qq[
        <a href="http://psd-production.internal.sanger.ac.uk:6600/workflow_samples/$sample_ssid/items/$lib_ssid">[SS]</a>
        <a href="http://intweb.sanger.ac.uk/perl/prodsoft/npg/npg/search?query=$name">[NPG]</a>
        </h3>
		];
	}
    
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
        print qq[
            <form action="qcgrind_v2.pl">
            <input type="hidden" name="mode" value="$LIB_UPDATE">
            <input type="hidden" name="sp" value="$species">
            <input type="hidden" name="lib_id" value="$libID">
            <p>
            <table width="400" align="center">
            <tr>
            <td align="center"><input type="submit" name="lib_update" value="$PENDING_FILTER" class="btn" style="background-color: #F5F5F5;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
            <td align="center"><input type="submit" name="lib_update" value="$FAILED_FILTER" class="btn" style="background-color: #FFC0C0;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
            <td align="center"><input type="submit" name="lib_update" value="$PASSED_FILTER" class="btn" style="background-color: #C0FFC0;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
            ];
			
			if( $library->open() == 1 )
			{
				print qq[
					<td align="center"><input type="submit" name="lib_update" value="$CLOSE_LIBRARY" class="btn" style="background-color: #FFC0C0;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
				];
			}
			else
			{
				print qq[
					<td align="center"><input type="submit" name="lib_update" value="$OPEN_LIBRARY" class="btn" style="background-color: #C0FFC0;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
				];
			}
			
			print qq[
			</tr>
            <tr>
            <td align="left"></td>
            <td align="center"><b>$libcount</b> lanes in library</td>
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
    my $typeID = $library ->library_type_id();
    my $type = VRTrack::Library_type->new($vrtrack, $typeID);
    my $type_ = $type->name();
    print qq[
    <div class="centerFieldset">
    <table>
    <tr>
    <td>Insert</td>
    <td>$insert</td>
    <td>Type</td>
    <td>$type_</td>
    </tr>
    </table>
    </div>
    ];
    
    #get all the lanes for the library (max out at 20)
    my $lanes = $library->lanes();
    
    #sort by the lane name
    my @l = sort( { $a->name() cmp $b->name() } @{$lanes} );
    $lanes = \@l;
    
    my $count = 0;
    
    print qq[
    <div class="centerFieldset">
    <fieldset style="width: 1200px">
    <legend>Lane data</legend>
    <form action="$SCRIPT_NAME">
    <input type="hidden" name="mode" value="$LANES_UPDATE">
    <input type="hidden" name="sp" value="$species">
    <input type="hidden" name="lib_id" value="$libID">
    <input type="hidden" name="filter" value="$filter">
    <table width="100%">
    <tr>
    <th>Pass</th>
    <th>Fail</th>
    <th>Pend</th>
    <th>Name</th>
    <th>Auto QC</th>
    <th>Cycles</th>
    <th>Bases</th>
    <th>Post-clip</th>
    <th>Adapter</th>
    <th>Mapped</th>
    <th>Paired</th>
    <th>Rmdup</th>
    <th>Genotype</th>
    </tr>
    ];
    
    my $noQC = 0;
    foreach( @$lanes )
    {
        my $lane = $_;
        my ($project,$sample,$library) = get_hierarchy_for_lane($vrtrack, $lane);
        
        my $status = $lane->qc_status;
        my $name = $lane->name;
        my $mapstats = $lane->latest_mapping;
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
                <td><input type="radio" name="$id" value="$PASSED_FILTER"
                ];
                print $status eq $PASSED_FILTER ? 'checked' : '';
                print qq[></td>
                <td><input type="radio" name="$id" value="$FAILED_FILTER"];
                print $status eq $FAILED_FILTER ? 'checked' : '';
                print qq[></td>
                <td><input type="radio" name="$id" value="$PENDING_FILTER"];
                print $status eq $PENDING_FILTER ? 'checked' : '';
                print qq[></td>
                <td style="background-color:$lane_status_colour;"><a href="$SCRIPT_NAME?mode=$LANE_VIEW&amp;lane_id=$id&amp;sp=$species">$name</a></td>
                ];
                
                print defined($auto_qc_status) ? qq[<td style="background-color:$lane_auto_status_colour;">$auto_qc_status</td>] : qq[<td>undef</td>];

                print qq[
                <td>$cycles</td>
                <td>$bases</td>
                <td>$clip_bases_perc%</td>
                <td>$adapter_perc%</td>
                <td>$reads_mapped_perc%</td>
                <td>$reads_paired_perc%</td>
                <td>$rmdup_reads_mapped_perc%</td>
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
}

sub connectToDatabase
{
    my $database = $_[ 0 ];
    my $vrtrack = VRTrack::VRTrack->new({ host => 'mcs4a',
                                    port => 3306,
                                    user => 'vreseq_rw',
                                    password => 't3aml3ss',
                                    database => $database,
                                   });
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
    
    my $location = "qcgrind_v2.pl?mode=$ERROR_DISPLAY&amp;error_msg=$error";
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
		if ($status eq 'no_qc')
		{
			$status_colour="#FFFFFF";
		}
		if ($status eq 'passed')
		{
			$status_colour="#C0FFC0";
		}
		elsif ($status eq 'failed')
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
