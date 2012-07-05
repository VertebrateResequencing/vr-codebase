package VertRes::QCGrind::ViewUtil;
use base qw(VertRes::QCGrind::Util);
# QC Grind common variables and modules
use strict;

sub new {
    my $self={};

	$self->{SCRIPTS} = {
	    MAP_VIEW              => 'map_view.pl', 
		MAP_PROJECTS_VIEW	  => 'map_projects_view.pl',
		MAP_LANES_VIEW        => 'map_lanes_view.pl',
		PENDING_VIEW          => 'pending_view.pl',
		PENDING_REQ_VIEW	  => 'pending_requests_view.pl',
		SAMPLE_MAPPING        => 'sample_mapping.pl',
		SAMPMAP_PROJECTS_VIEW => 'sample_mapping_projects_view.pl',
		SAMPMAP_LANES_VIEW    => 'sample_mapping_lanes_view.pl',
		QCGRIND_LANE          => '../qc_grind/lane_view.pl',
		QCGRIND_SAMPLES       => '../qc_grind/samples_view.pl',
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

.coolfieldset, .coolfieldset.expanded{
	border:1px solid #aaa;   
}

.coolfieldset.collapsed{
	border:0;
	border-top:1px solid #aaa;
}

.coolfieldset legend{
	padding-left:13px;
	font-weight:bold;
	cursor:pointer;
}

.coolfieldset legend, .coolfieldset.expanded legend{
	background: transparent url(http://www.sanger.ac.uk/modelorgs/mousegenomes/images/expanded.gif) no-repeat center left;
}

.coolfieldset.collapsed legend{
	background: transparent url(http://www.sanger.ac.uk/modelorgs/mousegenomes/images/collapsed.gif) no-repeat center left;
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

tr:nth-child(2n+1) {
	background-color: #ecf1ef;
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


sub displayDatabasesPage {
    my ($self,$title,$cgi,$script,$alldb) = @_;

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
	if ($alldb) {
		my %done;
		push (@main_dbs, @uk10k_dbs);
		foreach (@main_dbs) {
			$done{$_}++;
		}
    	my @dbs = $self->fetchTrackingDatabases();
    	print qq[
        	<div class="centerFieldset">
        	<fieldset id="fieldset1" class="coolfieldset" style="width: 500px">
        	<legend>Other Databases (click to show/hide)</legend> <div>
    	];
    	foreach( @dbs ) {
			next if $done{$_};
			print $cgi->p("<a href='$script?db=$_'> $_ </a>");
    	}
    	print qq[ </div> </fieldset> </div>
    		<script>
			\$('#fieldset1').coolfieldset({collapsed:true});
			</script>];
	}
}

sub displayDatabasePage
{
    my ($self,$title,$cgi,$vrtrack,$db,$init_script,$lanes_script) = @_;

    print qq[
         <h2 align="center" style="font: normal 900 1.5em arial"><a href="$init_script">$title</a></h2>
         <h3 align="center" style="font: normal 700 1.5em arial">Database : $db</h3>
        <div class="centerFieldset">
        <fieldset style="width: 500px">
        <legend>Select study to View</legend>
    ];
    my @projects = sort {$a->name cmp $b->name} @{$vrtrack->projects()};
    foreach my $project (@projects)
    {
        my $id = $project->id();
        print qq[<p><a href="$lanes_script?db=$db&amp;proj_id=$id">].$project->name().qq[</a></p>];
    }
    print qq[
        </fieldset>
        </div>
    ];
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

sub getDatabaseID 
{
	my ($self, $db) = @_;
	my $db_id;
    my $web_db = 'vrtrack_web_index';
	my $vrtrack = $self->connectToDatabase($web_db);
	$self->displayError( "Failed to connect to web database: $web_db" ) unless defined( $vrtrack );
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
	my ($self, $db_id, $pid) = @_;
	my %mappings;
    my $web_db = 'vrtrack_web_index';
	my $vrtrack = $self->connectToDatabase($web_db);
	$self->displayError( "Failed to connect to web database: $web_db" ) unless defined( $vrtrack );
	my $sql = qq[SELECT supplier_name, accession_number, sanger_sample_name FROM sample_id_mapping where db_id = ? and project_id = ?];
	my $sth = $vrtrack->{_dbh}->prepare($sql);	
	if ($sth->execute($db_id, $pid)) {
		my ($supp, $acc, $sang);
		$sth->bind_columns(\($supp, $acc, $sang));
		while ($sth->fetch) {
			push @{ $mappings{$sang} }, ($supp, $acc);
		}
	}
	return %mappings;	
}

sub fetchProjectName
{
	my ($self, $db_id, $pid) = @_;
    my $pname;
    my $web_db = 'vrtrack_web_index';
	my $vrtrack = $self->connectToDatabase($web_db);
	$self->displayError( "Failed to connect to web database: $web_db" ) unless defined( $vrtrack );
	my $sql = qq[SELECT project_name FROM db_projects where project_id = ? and db_id =?];
	my $sth = $vrtrack->{_dbh}->prepare($sql);
	if ($sth->execute($pid, $db_id)) {
		my ($name);
		$sth->bind_col(1, \$name);
		while ($sth->fetch) {
			$pname = $name;
		}
	}
	return $pname;
}

1;
