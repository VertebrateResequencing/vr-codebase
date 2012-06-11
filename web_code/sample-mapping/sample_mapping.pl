#!/usr/local/bin/perl -T
# Displays sample id mappings (Sanger sample name, Supplier name, Accession number) for all major vrtrack databases
# Data can be downloaded in tsv or csv (for excel) formats
#
# Author:        jm23 
# Maintainer:    jm23
# Created:       2011-05-08

use strict;
use warnings;
use URI;

BEGIN {
   $ENV{VRTRACK_HOST} = 'mcs4a';
   $ENV{VRTRACK_PORT} = 3306;
   $ENV{VRTRACK_RO_USER} = 'vreseq_ro';
   $ENV{VRTRACK_RW_USER} = 'vreseq_rw';
   $ENV{VRTRACK_PASSWORD} = 't3aml3ss';
};

use lib '.';
use SangerPaths qw(core team145);
#use SangerPaths qw(core);
#use lib './modules';
use VRTrack::VRTrack;
use VRTrack::Project;
use VRTrack::Sample;
use VRTrack::Library;
use VRTrack::Lane;
use VRTrack::Multiplex_pool;
use VertRes::Utils::VRTrackFactory;
use Data::Dumper;
use DBI;

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

/* Sortable tables */
table.sortable thead {
    background-color:#eee;
    color:#666666;
    font-weight: bold;
    cursor: default;
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
    'title'   => q(Sample Mapping v1),
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
    #'stylesheet' => '/Teams/Team145/view.css'
    'jsfile'  => 'http://js.sanger.ac.uk/sorttable_v2.js',
    'style'   => $css,
});

#script name for self links
my $cgi = $sw->cgi();
my $SCRIPT_NAME = $cgi->url(-relative=>1);
my $web_db = 'vrtrack_web_index';

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

my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $web_db, mode => 'r');
redirectErrorScreen( $cgi, "Failed to connect to database: $web_db" ) unless defined( $vrtrack );

# List available databases.
if( !defined( $cgi->param('db')) ) {
    print $sw->header();
    displayDatabasesPage($vrtrack);
    print $sw->footer();
    exit;
}

# All other entry points require a database
my $db = $cgi->param('db');

if ($cgi->param('download')) {
	my $pid = $cgi->param('proj_id');
    if( defined( $pid ) ) {
        if ($cgi->param('download') eq 'TSV') {
    		&downloadMappingsTSV($cgi, $vrtrack, $db, $pid);
    		exit;
    	}
    	elsif ($cgi->param('download') eq 'CSV (Excel)') {
    		&downloadMappingsCSV($cgi, $vrtrack, $db, $pid);
    		exit;
    	}
    }	
}

# if ($cgi->param('csvfile')) {
# 	my $pid = $cgi->param('proj_id');
#     if( defined( $pid ) ) {
#     	&downloadMappingsCSV($cgi, $vrtrack, $db, $pid);
#     	exit;
#     }	
# }

if( ! defined $db ) {
    redirectErrorScreen( $cgi, "Database must be defined!" );
    exit;
}

if( ! isDatabase( $db, $vrtrack ) ) {
    redirectErrorScreen( $cgi, "Invalid database name!" );
    exit;
}

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
    displayProjectMappingsPage($cgi, $vrtrack, $db, $pid);
    print $sw->footer();
    exit;
}

else
{
    redirectErrorScreen( $cgi, "Invalid mode!" );
}



sub displayDatabasesPage 
{
    my $vrtrack = shift;

   print qq[
        <h2 align="center" style="font: normal 900 1.5em arial">Sample ID Mapping</h2>
        <div class="centerFieldset">
        <fieldset style="width: 500px">
        <legend>VRTrack databases</legend>
    ];
    my @vrtrack_dbs = qw(vrtrack_human_wgs vrtrack_human_wes vrtrack_mouse_wgs vrtrack_mouse_wes vrtrack_cerevisiae_wgs);
    foreach( @vrtrack_dbs )
    {
        print qq[
            <p><a href="$SCRIPT_NAME?mode=$DB_VIEW&amp;db=$_">$_</a></p>
        ];
    }
    print qq[
        </fieldset>
        </div>
        <div class="centerFieldset">
        <fieldset style="width: 500px">
        <legend>UK10K Databases</legend>
    ];
    my @uk10k_dbs = qw(vrtrack_uk10k_cohort vrtrack_uk10k_neuro vrtrack_uk10k_obesity vrtrack_uk10k_rare);
    foreach( @uk10k_dbs )
    {
        print qq[
            <p><a href="$SCRIPT_NAME?mode=$DB_VIEW&amp;db=$_">$_</a></p>
        ];
    }    
    print qq[ </fieldset>
              </div> ];
}

sub displayDatabasePage
{
    my $cgi = shift;
    my $vrtrack = shift;
    my $db = shift;
    
    my %projects = fetchProjects($vrtrack, $db);
    my @proj_names = keys %projects;
    #if scalar @proj_names == 1 then go straight to map view
    if (scalar @proj_names == 1 ) {
    	$mode =$PROJ_VIEW;
    	my $proj_name = $proj_names[0];
        displayProjectMappingsPage($cgi, $vrtrack, $db, $projects{$proj_name});
    }
    else {
    	print qq[
         <h2 align="center" style="font: normal 900 1.5em arial"><a href="$SCRIPT_NAME">Sample ID Mapping</a></h2>
        <div class="centerFieldset">
        <fieldset style="width: 500px">
        <legend>Select study to View</legend>
    	];
    
    	foreach (@proj_names)
    	{
        	print qq[
            	<p><a href="$SCRIPT_NAME?mode=$PROJ_VIEW&amp;proj_id=$projects{$_}&amp;db=$db">$_</a></p>
        	];
    	}
    
    	print qq[
        	</fieldset>
        	</div>
    	];
	}
}

sub displayProjectMappingsPage 
{
    my ($cgi, $vrtrack, $database, $projectID) = @_;
    
    my $db_id = getDatabaseID ($vrtrack, $database);
        
    my %sampleMappings = getSampleMappings ($vrtrack, $db_id, $projectID);

    my $pname = fetchProjectName($vrtrack, $projectID, $db_id);
    print qq[
    <h2 align="center" style="font: normal 900 1.5em arial"><a href="$SCRIPT_NAME">Sample ID Mappings</a></h2>
    <h5 style="font: arial"><p><a href="$SCRIPT_NAME?mode=$DB_VIEW&amp;db=$database">$database</a>: $pname</p></h5>];
    
    print qq[<div class="centerFieldset">
    <fieldset style = "border:none">
    <table width="60%">];
    
    print $cgi->start_form;
    print $cgi->hidden("db","$database");
    print $cgi->hidden("proj_id","$projectID");
    print qq[
    <tr></tr>
    <tr>
    <td colspan = 3 align='right'>Download:];
    print $cgi->submit(-name => 'download', -value  => 'TSV');
    print $cgi->submit(-name => 'download', -value  => 'CSV (Excel)');
    print qq[
    	</td>
    	</tr>
    	</table>
    	</fieldset>
        </div>];
    print qq[
    <div class="centerFieldset">
    <fieldset > 
    <legend>Sample mappings</legend>
    <table class='sortable' width="60%">
    <tr>
    <th>Sanger sample name</th>
    <th>Supplier name</th>
    <th>Accession number</th>
    </tr>
    ];
    
    foreach my $sang ( sort keys %sampleMappings )
    {
    	my @sampData = @{ $sampleMappings{$sang} };
        print qq[<tr><td align="center">$sang</td><td align="center">$sampData[0]</td><td align="center">$sampData[1]</td></tr>];
    }
    print qq[
        </table>
        </fieldset>
        </div>
    ];
    print $cgi->end_form;
}


sub isDatabase
{
        my $db = shift;
        my @dbs = fetchTrackingDatabases($vrtrack);
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

sub fetchTrackingDatabases
{
	my $vrtrack = shift;
    my @dbs;
	my $sql = qq[SELECT db_name FROM tracking_database];
	my $sth = $vrtrack->{_dbh}->prepare($sql);
	if ($sth->execute()) {
		my ($col1);
		$sth->bind_col(1, \$col1);
		while ($sth->fetch) {
			push @dbs, $col1;
		}
	}
	return @dbs;
}

sub fetchProjects
{
	my $vrtrack = shift;
	my $db = shift;
    my %projects;
	my $sql = qq[SELECT distinct s.project_id, p.project_name FROM db_projects p, sample_id_mapping s, tracking_database d where p.project_id = s.project_id and s.db_id = p.db_id and p.db_id = d.db_id and d.db_name = ?];
	my $sth = $vrtrack->{_dbh}->prepare($sql);
	if ($sth->execute($db)) {
		my ($id, $name);
		$sth->bind_columns(\($id, $name));
		while ($sth->fetch) {
			$projects{$name} = $id;
		}
	}
	return %projects;
}

sub fetchProjectName
{
	my $vrtrack = shift;
	my $pid = shift;
	my $dbid = shift;
    my $pname;
	my $sql = qq[SELECT project_name FROM db_projects where project_id = ? and db_id =?];
	my $sth = $vrtrack->{_dbh}->prepare($sql);
	if ($sth->execute($pid, $dbid)) {
		my ($name);
		$sth->bind_col(1, \$name);
		while ($sth->fetch) {
			$pname = $name;
		}
	}
	return $pname;
}

sub getDatabaseID 
{
	my $vrtrack = shift;
	my $db = shift;
	my $db_id;
    my $sql = qq[SELECT db_id FROM tracking_database where db_name = ?];
	my $sth = $vrtrack->{_dbh}->prepare($sql);
	if ($sth->execute($db)) {
		my ($id);
		$sth->bind_col(1, \$id);
		while ($sth->fetch) {
			$db_id = $id;
		}
	}
	return $db_id;
}

sub getSampleMappings 
{
	my $vrtrack = shift;
	my $db_id = shift;
	my $projectID = shift; 
	my %mappings;
	my $sql = qq[SELECT supplier_name, accession_number, sanger_sample_name FROM sample_id_mapping where db_id = ? and project_id = ?];
	my $sth = $vrtrack->{_dbh}->prepare($sql);
	if ($sth->execute($db_id, $projectID)) {
		my ($supp, $acc, $sang);
		$sth->bind_columns(\($supp, $acc, $sang));
		while ($sth->fetch) {
			push @{ $mappings{$sang} }, ($supp, $acc);
		}
	}
	return %mappings;	
}

sub downloadMappingsTSV
{
	my ($cgi, $vrtrack, $db, $projectID) = @_;
	my $db_id = getDatabaseID ($vrtrack, $db);
	my $pname = fetchProjectName($vrtrack, $projectID, $db_id);
	print $cgi->header(-type=>'text/tsv',  -attachment=>"$pname.tsv");
	print join ("\t","Sanger Sample Name","Supplier Name","Accession"), "\n";
	my $sql = qq[SELECT supplier_name, accession_number, sanger_sample_name FROM sample_id_mapping where db_id = ? and project_id = ?];
	my $sth = $vrtrack->{_dbh}->prepare($sql);
	if ($sth->execute($db_id, $projectID)) {
		my ($supp, $acc, $sang);
		$sth->bind_columns(\($supp, $acc, $sang));
		while ($sth->fetch) {
		    print (join ("\t",$sang, $supp, $acc), "\n"); 
		}
	}
}

sub downloadMappingsCSV
{
	my ($cgi, $vrtrack, $db, $projectID) = @_;
	my $db_id = getDatabaseID ($vrtrack, $db);
	my $pname = fetchProjectName($vrtrack, $projectID, $db_id);
	print $cgi->header(-type=>'text/csv',  -attachment=>"$pname.csv");
	print join (",","Sanger Sample Name","Supplier Name","Accession"), "\n";
	my $sql = qq[SELECT supplier_name, accession_number, sanger_sample_name FROM sample_id_mapping where db_id = ? and project_id = ?];
	my $sth = $vrtrack->{_dbh}->prepare($sql);
	if ($sth->execute($db_id, $projectID)) {
		my ($supp, $acc, $sang);
		$sth->bind_columns(\($supp, $acc, $sang));
		while ($sth->fetch) {
		    print (join (",",$sang, $supp, $acc), "\n"); 
		}
	}
}