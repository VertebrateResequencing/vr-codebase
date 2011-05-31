#!/usr/local/bin/perl -T
# Displays pending sequencing requests
#
# Author:        tk2/jm23
# Maintainer:    ??
# Created:       2011-04-08

use strict;
use warnings;
use URI;

use lib '.';
#use SangerPaths qw(core team145);
use SangerPaths qw(core);
use lib './modules';
use VRTrack::VRTrack;
use VRTrack::Project;
use VRTrack::Sample;
use VRTrack::Library;
use VRTrack::Lane;
use VertRes::Utils::VRTrackFactory;
use Data::Dumper;

use SangerWeb;
$ENV{PATH}= '/usr/local/bin';   # solely to stop taint from barfing

$|++;

#different modes/views possible
my $ERROR_DISPLAY = 0;
my $PROJ_VIEW = 1;
my $DB_VIEW = 2;

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

tr:nth-child(2n+1) {
	background-color: #ecf1ef;
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

my $sw  = SangerWeb->new({
    'title'   => q(Pending View v1),
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
    'style'   => $css,
});

#script name for self links
my $cgi = $sw->cgi();
my $SCRIPT_NAME = $cgi->url(-relative=>1);
my $qc_grind = "../qc_grind/qc_grind.pl";

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

if( $mode == $DB_VIEW )
{
    print $sw->header();
    displayPendingRequestsPage( $cgi, $vrtrack, $db ); 
    print $sw->footer();
    exit;
}
# elsif( $mode == $PROJ_VIEW ) 
# {
#     my $pid = $cgi->param('proj_id');
#     if( ! defined( $pid ) )
#     {
#         redirectErrorScreen( $cgi, "Must provide a project ID" );
#         exit;
#     }
#     
#     print $sw->header();
# 	displayPendingRequestsPage($cgi, $vrtrack, $db, $pid);
#     print $sw->footer();
#     exit;
# }
else
{
    redirectErrorScreen( $cgi, "Invalid mode!" );
}

sub displayDatabasesPage 
{
    print qq[
        <h2 align="center" style="font: normal 900 1.5em arial">Pending View</h2>
        <div class="centerFieldset">
        <fieldset style="width: 500px">
        <legend>Select dataset to View</legend>
    ];
    
    my @dbs = VertRes::Utils::VRTrackFactory->databases(1);
    foreach( @dbs )
    {
        print qq[
            <p><a href="$SCRIPT_NAME?mode=$DB_VIEW&amp;db=$_">].$_.qq[</a></p>
        ];
    }
    
    print qq[
        </fieldset>
        </div>
    ];
}

sub displayPendingRequestsPage 
{

    my $cgi = shift;
    my $vrtrack = shift;
    my $db = shift;    

	my @projects = @{ $vrtrack->projects() };

	my @libpendings;
	my @seqpendings;
	my $lib_type = "library";
	my $seq_type = "sequence";

	print qq[
    	<h2 align="center" style="font: normal 900 1.5em arial"><a href="$SCRIPT_NAME">Pending Requests</a></h2>
		<div class="centerFieldset">
    ];

    foreach my $project (@projects)
    {
        my $projectID = $project->id();
    
    my $samples = $project->samples();
    displayError( "Cant get samples for project: $projectID" ) unless $samples;
     
    my $pname = $project->name;

    my $pendingflag = "pending";

    foreach( @{$samples} ) 
    {	
		my $sample = $_;
		my $sname = $sample->name();
		my $library_requests = $sample->library_requests();
		foreach ( @{$library_requests}) {
			my $librequest = $_;
			my $lib_status = $librequest->prep_status();
			if ($lib_status eq $pendingflag) {
				my $lib_date = $librequest->changed();
				my $ssid = $librequest->ssid();
				push @libpendings, [$projectID, $pname, $sname, $lib_status, $lib_date, $ssid];
			}
		}
		my $libraries = $sample->libraries();
		foreach ( @{$libraries}) {
			my $library = $_;
			my $seq_requests = $library->seq_requests();
			foreach ( @{$seq_requests} ) {
				my $seqrequest = $_;
				my $seq_status = $seqrequest->seq_status();
				if ($seq_status eq $pendingflag) {
					my $seq_date = $seqrequest->changed();
					my $ssid = $seqrequest->ssid();
					push @seqpendings, [$projectID, $pname, $sname, $seq_status, $seq_date, $ssid];
				}
			}
		}
    }
	}
	printPendingRows(\@libpendings, $lib_type, $db, $cgi);
	printPendingRows(\@seqpendings, $seq_type, $db, $cgi);
    print qq[
    	</div>
    ];
}

sub printPendingRows
{
	my ($pendinglist, $req_type, $db, $cgi) = @_;
    print qq[
    <fieldset > 
    <legend>Pending $req_type requests</legend>
    <table width="80%">
    <tr>
    <th>QC Link</th>
    <th>Project</th>
    <th>Sample Name</th>
	<th>Status</th>
	<th>Date changed</th>
	<th>SSID</th>
	</tr>
    ];
	my @pendings = @$pendinglist;
	my $current_project = '';
  	if (@pendings) {
		for my $i ( 0 .. $#pendings ) {
			my @pending = @{$pendings[$i]};
			print qq[<tr>];
			if ($current_project ne $pending[1]) {
				print qq[<td><a href="$qc_grind?db=$db&mode=6&proj_id=$pending[0]">QC_grind</a></td>];
				$current_project = $pending[1];
			}
			else {
				print qq[<td></td>];
			}			
			for my $j ( 1 .. ($#pending-1) ) {
				print qq[<td>$pending[$j]</td>];
			}
			print qq[<td><a href="http://psd-production.internal.sanger.ac.uk:6600/requests/$pending[$#pending]">$pending[$#pending]</a></td></tr>];
		}
	}
  	else {
		print qq[<tr><td colspan=5>No pending $req_type requests were found in the $db database.</td></tr>];
	}	
    print qq[
        </table>
        </fieldset><br/>
    ];
}

sub isDatabase
{
        my $db = shift;
        
        my @dbs = VertRes::Utils::VRTrackFactory->databases(1);
        foreach( @dbs ){if( $db eq $_ ){return 1;}}
        return 0;
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
   # print $cgi->header(); 
   # print "window.location=\"$location\";\n\n";
    print $cgi->redirect( -URL => $location, -method   => 'GET', -status   => 302 );
    #print "Location: $location";
    exit;
}
