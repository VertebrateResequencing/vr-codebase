package VertRes::QCGrind::ViewUtil;
use base qw(VertRes::QCGrind::Util);
#Extension of VertRes::QCGrind::Util for other tracking database view webpages (map_view, pending_view, sample_id_mapping)

use strict;

sub new {
    my $self={};

	$self->{SCRIPTS} = {
		MAP_DB_VIEW => 'map_database.pl',
		MAP_PROJ_VIEW	=> 'map_project.pl',
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

    bless $self;
    return $self;
}

sub displayDatabasesPage 
{
   	my ($self, $title, $script) = @_;
   	print qq[ <h2 align="center" style="font: normal 900 1.5em arial">$title</h2> ];
    
    my @main_dbs = qw (vrtrack_human_wgs vrtrack_human_wes vrtrack_mouse_wgs vrtrack_mouse_wes);
    print qq[
       <div class="centerFieldset">
       <fieldset style="width: 500px">
       <legend>Main Databases</legend>
    ];
    foreach( @main_dbs ) {
   		print $cgi->p("<a href='$script?db=$_'> $_ </a>");
    }
    print qq[ </fieldset> </div> ];

	my @uk10k_dbs = qw (vrtrack_uk10k_cohort vrtrack_uk10k_neuro vrtrack_uk10k_obesity vrtrack_uk10k_rare);
    print qq[
        <div class="centerFieldset">
        <fieldset style="width: 500px">
        <legend>UK10K Databases</legend>
    ];
    foreach( @uk10k_dbs ) {
		print $cgi->p("<a href='$script?db=$_'> $_ </a>");
    }
    print qq[ </fieldset> </div> ];
    
    my %done;
	push (@main_dbs, @uk10k_dbs);
	foreach (@main_dbs) {
		$done{$_}++;
	}
    
    my @dbs = fetchTrackingDatabases();
    print qq[
        <div class="centerFieldset">
        <fieldset style="width: 500px">
        <legend>Other Databases</legend>
    ];
    foreach( @dbs ) {
		next if $done{$_};
		print $cgi->p("<a href='$script?db=$_'> $_ </a>");
    }
    print qq[ </fieldset> </div> ];
}

sub fetchTrackingDatabases
{
	my ($self) = @_;
	my @dbs;
	my $web_db = 'vrtrack_web_index';
	my $vrtrack = $self->connectToDatabase($web_db);
	$self->displayError( "Failed to connect to web database: $web_db" ) unless defined( $vrtrack );
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


1;
