package VertRes::QCGrind::Util;
# QC Grind common variables and modules
use strict;

sub new {
    my $self={};

	$self->{SCRIPTS} = {
		DATABASES_VIEW => 'databases_view.pl',
		PROJECTS_VIEW	=> 'projects_view.pl',
		STUDY_LANES	=> 'study_lanes.pl',
		SAMPLES_VIEW	=> 'samples_view.pl',
		LANE_VIEW	=> 'lane_view.pl',
	};

	#possible states
	$self->{STATES} = {
		PASSED => 'passed',
		PENDING => 'pending',
		NO_QC => 'no_qc',
		FAILED => 'failed',
		GT_PENDING => 'gt_pending',
		INVESTIGATE => 'investigate',
		CLOSE_LIBRARY => 'close',
		OPEN_LIBRARY => 'open',
	};

	# Use SSO to authenticate, and then authorise against a list of people allowed to update the database.
	$self->{AUTH_USERS} = {
		'jws' => 'Jim Stalker',
		'tk2' => 'Thomas Keane',
		'pd3' => 'Petr Danecek',
		'sb10'=> 'Sendu Bala',
		'rd'  => 'Richard Durbin',
		'kb1' => 'Karen McLaren',
		'ylx' => "Yali Xue from Chris's group",
		'ak6' => 'Anja Kolb-Kokocinski (kuusamo)',
		'cj5' => 'chris joyce',
	};

    $self->{CSS} = <<CSS ;

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

    bless $self;
    return $self;
}

sub displayError
{
	my ($self,$message,$sw) = @_;
    
    print qq[<h2>A problem occurred</h2>\n];
    print qq[<p class="error1">$message</p>\n];
    print $sw->footer();
    exit;
}

sub connectToDatabase
{
	my ($self,$db) = @_;

    my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $db, mode => 'rw');
    return $vrtrack;
}

sub redirectErrorScreen
{
	my ($self,$cgi,$SCRIPT_NAME) = @_;

    my $error = 'An unknown error occurred';
    if( @_ == 2 )
    {
        $error = $_[ 1 ];
    }
    
	my $mode = $self->{MODES}{ERROR_DISPLAY};
    my $location = "$SCRIPT_NAME?mode=$mode&amp;error_msg=$error";
    #print $cgi->header(); 
    #print "window.location=\"$location\";\n\n";
    print $cgi->redirect( -URL => $location, -method   => 'GET', -status   => 302 );
    #print "Location: $location";
    exit;
}

# returns CSS colour for a QC status
sub get_colour_for_status
{
	my ($self,$status) = @_;

	my $status_colour;

	if( $status eq $self->{STATES}{CLOSE_LIBRARY}) {
		$status_colour="#FFBF00";
	}
	elsif ($status eq $self->{STATES}{NO_QC}) {
		$status_colour="#FFFFFF";
	}
	elsif ($status eq $self->{STATES}{PASSED}) {
		$status_colour="#CCFFCC";
	}
	elsif ($status eq $self->{STATES}{FAILED}) {
		$status_colour="#FFC0C0";
	}
	else {
		$status_colour="#F5F5F5";
	}
	return $status_colour;
}

# add commas in to big numbers
sub commify { 
	my $self = shift;
    my $number = shift;

    $number =  reverse $number;
    $number =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g; 
    return scalar reverse $number; 
}

sub printFullImage
{
	my $self = shift;
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
	my $self = shift;
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
	my $self = shift;
	my $db = shift;

	my @dbs = VertRes::Utils::VRTrackFactory->databases(1);
	foreach( @dbs ){if( $db eq $_ ){return 1;}}
	return 0;
}

sub bp_to_nearest_unit
{
	my ($self,$bp,$dp) = @_;
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
1;
