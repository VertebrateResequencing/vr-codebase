#!/usr/local/bin/perl -T
# Displays QC information for Sanger short-read sequencing to allow lanes
# to be passed/failed
#
# Author:        tk2    
# Maintainer:    tk2
# Created:       2011-03-12

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
    'title'   => q(Map View v1),
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
    'style'   => $css,
});

#script name for self links
my $cgi = $sw->cgi();
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

if( $mode == $DB_VIEW )
{
    print $sw->header();
    displayDatabasePage( $cgi, $vrtrack, $db ); 
    print $sw->footer();
    exit;
}
elsif( $mode == $PROJ_VIEW ) 
{
    my $pid = $cgi->param('proj_id');
    if( ! defined( $pid ) )
    {
        redirectErrorScreen( $cgi, "Must provide a project ID" );
        exit;
    }
    
    print $sw->header();
    displayProjectLanesPage($cgi, $vrtrack, $db, $pid);
    print $sw->footer();
    exit;
}
else
{
    redirectErrorScreen( $cgi, "Invalid mode!" );
}

sub displayDatabasesPage 
{
    print qq[
        <h2 align="center" style="font: normal 900 1.5em arial">Map View</h2>
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

sub displayDatabasePage
{
    my $cgi = shift;
    my $vrtrack = shift;
    my $db = shift;
    print qq[
         <h2 align="center" style="font: normal 900 1.5em arial"><a href="$SCRIPT_NAME">Map View</a></h2>
        <div class="centerFieldset">
        <fieldset style="width: 500px">
        <legend>Select study to View</legend>
    ];
    
    my @projects = @{ $vrtrack->projects() };
    foreach my $project (@projects)
    {
        my $id = $project->id();
        print qq[
            <p><a href="$SCRIPT_NAME?mode=$PROJ_VIEW&amp;proj_id=$id&amp;db=$db">].$project->name().qq[</a></p>
        ];
    }
    
    print qq[
        </fieldset>
        </div>
    ];
}

sub displayProjectLanesPage 
{
    my ($cgi, $vrtrack, $database, $projectID) = @_;
    
    my $project = VRTrack::Project->new( $vrtrack, $projectID );
    displayError( "Cant get project: $projectID" ) unless $project;
    
    my $samples = $project->samples();
    displayError( "Cant get samples for project: $projectID" ) unless $samples;
    
    my $pname = $project->name;
    print qq[
    <h2 align="center" style="font: normal 900 1.5em arial"><a href="$SCRIPT_NAME">Map View</a></h2>
        ];
    
    print qq[
    <div class="centerFieldset">
    <fieldset > 
    <legend>Lane data</legend>
    <table width="60%">
    <tr>
    <th>Library</th>
    <th>Name</th>
    ];
    
    my @lanes;
    my %mappers;
    foreach( sort { $a->ssid() <=> $b->ssid() } @$samples)
    {
        my @libraries = sort {$a->name cmp $b->name} @{$_->libraries()};
        foreach( @libraries ) 
        {
            my $library = $_;
            my  @lanes_ = @{ $library->lanes() };
            foreach my $lane ( @lanes_ ) 
            {
                push( @lanes, $lane );
                my @mapstats = @{ $lane->mappings() };
                foreach my $mapstat (@mapstats)
                {
#                    print $mapstat->mapper()->name();exit;
                    if( $mapstat->mapper() )
                    {
                        my $mapper = $mapstat->mapper();
                        $mappers{ $mapper->name().qq[_].$mapper->version().qq[_].$mapper->id() } = 1;
                    }
                }
            }
        }
    }
    
    foreach my $mname ( sort( keys( %mappers ) ) )
    {
        print qq[<th>$mname</th>];
    }
    print qq[</tr>];
    
    foreach my $lane ( @lanes )
    {
        print qq[<tr>];
        my $id = $lane->id();
        my $library = VRTrack::Library->new($vrtrack, $lane->library_id());my $libName = $library->name();
        print qq[<td>$libName</td><td><a href="http://intwebdev.sanger.ac.uk/cgi-bin/teams/team145/qc_grind/qc_grind.pl?mode=0&lane_id=$id&db=$database">].$lane->name().qq[</a></td>];
        my @mappings = @{ $lane->mappings() };
        my %lane_mappers;
        foreach my $mapstat ( @mappings )
        {
            if( $mapstat->mapper() && $mapstat->raw_bases() && $mapstat->raw_bases() > 0 )
            {
                my $mapper = $mapstat->mapper();
                $lane_mappers{ $mapper->name().qq[_].$mapper->version().qq[_].$mapper->id() } = 1;
            }
        }
        foreach my $mapper ( sort( keys( %mappers ) ) )
        {
            if( $lane_mappers{ $mapper } )
            {
                print qq[<td>yes</td>];
            }
            else{print qq[<td>no</td>];}
        }
        print qq[</tr>];
    }
    print qq[
        </table>
        </fieldset>
        </div>
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
    #print $cgi->header(); 
    #print "window.location=\"$location\";\n\n";
    print $cgi->redirect( -URL => $location, -method   => 'GET', -status   => 302 );
    #print "Location: $location";
    exit;
}
