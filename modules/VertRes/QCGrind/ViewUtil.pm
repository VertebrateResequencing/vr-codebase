package VertRes::QCGrind::ViewUtil;
use base qw(VertRes::QCGrind::Util);
use VertRes::Utils::VRTrackFactory;
use DBI;
# Common variables and modules for view webpages (i.e. non-QCGrind)
use strict;

sub new {
	
    $ENV{VRTRACK_HOST} = 'mcs10';
    $ENV{VRTRACK_PORT} = 3306;
    $ENV{VRTRACK_RO_USER} = 'vreseq_ro';
    $ENV{VRTRACK_RW_USER} = 'vreseq_rw';
    $ENV{VRTRACK_PASSWORD} = 't3aml3ss';	
	
    my $self={};

	$self->{SCRIPTS} = {
	    MAP_VIEW                  => 'map_view.pl', 
		MAP_PROJECTS_VIEW	      => 'map_projects_view.pl',
		MAP_LANES_VIEW            => 'map_lanes_view.pl',
		PENDING_VIEW              => 'pending_view.pl',
		PENDING_REQ_VIEW	      => 'pending_requests_view.pl',
		SAMPLE_MAPPING            => 'sample_mapping.pl',
		SAMPMAP_PROJECTS_VIEW     => 'sample_mapping_projects_view.pl',
		SAMPMAP_LANES_VIEW        => 'sample_mapping_lanes_view.pl',
		DISK_USAGE                => 'vrpipe_disk_usage.pl',
		DISK_USAGE_DIRECTORY      => 'directory_disk_usage.pl',
		DISK_USAGE_PIPELINE       => 'pipeline_disk_usage.pl',
		DISK_USAGE_INDIV_DIR      => 'individual_directory_usage.pl',
		DISK_USAGE_INDIV_PIPE     => 'individual_pipeline_usage.pl',		
		DISK_USAGE_FILES          => 'file_disk_usage.pl',
		QCGRIND_LANE              => '../qc_grind/lane_view.pl',
		QCGRIND_SAMPLES           => '../qc_grind/samples_view.pl',
		INDEX_PAGE                => '../index.pl',
		HIPSCI_COHORTS_VIEW       => 'cohorts_view.pl',
		HIPSCI_SUMMARY_VIEW       => 'cohort_summary_view.pl',
		HIPSCI_DETAILED_VIEW      => 'cohort_detailed_view.pl',		
	};
	
	
	$self->{SQL} = {
	    PENDING_LIB              => "select distinct p.name, s.name, r.prep_status, r.changed, r.ssid from study y, latest_project p, latest_sample s, latest_library_request r where y.study_id = p.study_id and p.project_id = s.project_id and s.sample_id = r.sample_id and r.prep_status in ('pending', 'started') and y.name = ? order by p.name, r.ssid",
	    PENDING_SEQ              => "select distinct p.name, s.name, r.seq_status, r.changed, r.ssid from study y, latest_project p, latest_sample s, latest_library l, latest_seq_request r where y.study_id = p.study_id and p.project_id = s.project_id and s.sample_id = l.sample_id and l.library_id = r.library_id and r.seq_status in ('pending', 'started') and y.name = ? order by p.name, r.ssid",
	    PENDING_MULTIPLEX_SEQ    => "select distinct p.name, s.name, r.seq_status, r.changed, r.ssid from study y, latest_project p, latest_sample s, latest_library l, library_multiplex_pool m, latest_seq_request r where y.study_id = p.study_id and p.project_id = s.project_id and s.sample_id = l.sample_id and l.library_id = m.library_id and m.multiplex_pool_id = r.multiplex_pool_id and r.seq_status in ('pending', 'started') and y.name = ?",
	    DISK_USAGE_DIRECTORY     => "select path_root, file_type, s, s_gb, file_count from vrpipe_root_top_level where s_gb > 0 and top_level_display = 1 order by s_gb desc",
	    DISK_USAGE_PIPELINE      => "select ps_id, ps_name, ps_user, ps_type, total_s_gb, file_count from vrpipe_usage_total where total_s_gb > 0 order by total_s_gb desc",
	    DISK_USAGE_INDIV_DIR     => "select ps_id, ps_name, ps_user, ps_type, file_type, path_root, s_gb, file_count from vrpipe_usage_top_level where path_root = ? and s_gb > 0 order by s_gb desc", 
	    DISK_USAGE_INDIV_PIPE    => "select ps_id, ps_name, ps_user, ps_type, file_type, path_root, s_gb, file_count from vrpipe_usage_top_level where ps_id = ? and s_gb > 0 order by s_gb desc", 
	    DISK_USAGE_FILES         => "select ps_id, ps_name, ps_user, ps_type, s, type, path from vrpipe_file_info where ps_id = ? and path_root = ? and type = ? order by s desc",
	    HIPSCI_FETCH_COHORT      => "SELECT distinct i.name from individual i, latest_sample s, latest_library b, latest_lane l where i.individual_id = s.individual_id and s.sample_id = b.sample_id and b.library_id = l.library_id and l.withdrawn is null",
	    HIPSCI_FETCH_SAMPLES     => "SELECT s.name, s.note_id from latest_sample s, individual i where s.individual_id = i.individual_id and i.name = ?",
	    HIPSCI_FETCH_CNV_TOTALS  => "SELECT distinct s.name, m.raw_bases, m.clip_bases, m.gt_expected from latest_mapstats m, latest_lane l, latest_library b, latest_sample s, individual i where m.lane_id=l.lane_id and l.library_id=b.library_id and b.sample_id=s.sample_id and  s.individual_id = i.individual_id and i.name = ?",
	    HIPSCI_FETCH_LANE_ID     => "SELECT min(l.lane_id) from latest_lane l, latest_library b, latest_sample s, individual i where l.library_id = b.library_id and b.sample_id = s.sample_id and s.individual_id = i.individual_id and i.name = ?",
    };
    
    $self->{VRTRACK_DATABASES} = {
		GENOTYPING => 'vrtrack_hipsci_qc1_genotyping',
		EXPRESSION => 'vrtrack_hipsci_qc1_expression',
		WEB_TABLES => 'vrtrack_web_tables',
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
top: -600px;
left: -480px; /*position where enlarged image should offset horizontally */

}
CSS

    bless $self;
    return $self;
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
        <span><img src="$uri" alt="$caption" style="border:5px solid #83A4C3;"/></span></a>
    ];
}

sub displayDatabasesPage {
    my ($self,$title,$cgi,$script,$alldb,$pending) = @_;
	my $index = $self->{SCRIPTS}{INDEX_PAGE};
	print qq[ <h4 align="center" style="font: arial"><i><a href="$index">Team 145</a></i> : $title</h4> ];
    my $pending_db = 'vrtrack_pending_requests';
    my @main_dbs = qw (vrtrack_mouse_wes_GRCm38 vrtrack_mouse_wgs_GRCm38 vrtrack_human_wes_v5 vrtrack_human_wgs vrtrack_human_wes vrtrack_mouse_wgs vrtrack_mouse_wes vrtrack_scerevisiae_wgs);
    print qq[
        <div class="centerFieldset">
        <fieldset style="width: 500px">
        <legend>Main Databases</legend>
    ];
    my $db_href = $pending ? "<a href='$script?db=$pending_db&amp;dbpend=" : "<a href='$script?db=";
    foreach( @main_dbs ) {
		print $cgi->p($db_href."$_'> $_ </a>");
    }
    print qq[ </fieldset> </div> ];

	my @uk10k_dbs = qw (vrtrack_uk10k_obesity_wga vrtrack_uk10k_rare_wga vrtrack_uk10k_cohort vrtrack_uk10k_neuro vrtrack_uk10k_obesity vrtrack_uk10k_rare);
    print qq[
        <div class="centerFieldset">
        <fieldset style="width: 500px">
        <legend>UK10K Databases</legend>
    ];
    foreach( @uk10k_dbs ) {
		print $cgi->p($db_href."$_'> $_ </a>");
    }
    print qq[ </fieldset> </div> ];
	if ($alldb) {
		my %done;
		push (@main_dbs, @uk10k_dbs);
		foreach (@main_dbs) {
			$done{$_}++;
		}
        my $dbh = $self->nonVrtrackConnection('information_schema');
        my @dbs = @{$dbh->selectcol_arrayref('show databases')};
		@dbs = grep(!/information_schema/, @dbs);
    	@dbs = grep(!/test/, @dbs);
		@dbs = grep(!/jm23/, @dbs);
		@dbs = grep(!/tttt/, @dbs);
		@dbs = grep(!/dump/, @dbs);
		@dbs = grep(!/irods/, @dbs);
		@dbs = grep(!/kuusamo/, @dbs);
		@dbs = grep(!/requests/, @dbs);
		@dbs = grep(!/web_tables/, @dbs);
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
    my ($self,$title,$cgi,$vrtrack,$db,$init_script,$lanes_script,$all_flag) = @_;
    my $index = $self->{SCRIPTS}{INDEX_PAGE};
    print qq[
        <h4 align="center" style="font: arial"><i><a href="$index">Team 145</a> : <a href="$init_script">$title</a></i> :  $db</h4>
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
    print qq[<p><a href="$lanes_script?db=$db&amp;proj_id=99999">SHOW_ALL_PROJECTS</a></p>] if $all_flag;
    print qq[
        </fieldset>
        </div>
    ];
}

sub getDatabaseID 
{
	my ($self, $db) = @_;
	my $db_id;
	my $vrtrack = $self->connectToDatabase($self->{VRTRACK_DATABASES}{WEB_TABLES});
	$self->displayError( "Failed to connect to vrtrack web database" ) unless defined( $vrtrack );
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
	my $vrtrack = $self->connectToDatabase($self->{VRTRACK_DATABASES}{WEB_TABLES});
	$self->displayError( "Failed to connect to vrtrack web database" ) unless defined( $vrtrack );
	my $sql = qq[SELECT supplier_name, accession_number, sanger_sample_name FROM sample_id_mapping where db_id = ?];
	$sql = $sql." and project_id = ?" unless $pid == 99999;
	my $sth = $vrtrack->{_dbh}->prepare($sql);	
	if ( $pid == 99999 ) {
		$sth->execute($db_id);
	}
	else {
		$sth->execute($db_id, $pid);
	}
	while (my ($supp, $acc, $sang) = $sth->fetchrow_array()) {
		push @{ $mappings{$sang} }, ($supp, $acc);
	}
	return %mappings;	
}

sub fetchProjectName
{
	my ($self, $db_id, $pid) = @_;
    my $pname;
	my $vrtrack = $self->connectToDatabase($self->{VRTRACK_DATABASES}{WEB_TABLES});
	$self->displayError( "Failed to connect to vrtrack web database" ) unless defined( $vrtrack );
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

sub fetchUsageTotals
{
	my ($self, $sql) = @_;
	my $vrtrack = $self->connectToDatabase($self->{VRTRACK_DATABASES}{WEB_TABLES});
	$self->displayError( "Failed to connect to vrtrack web database" ) unless defined( $vrtrack );
	my $sth = $vrtrack->{_dbh}->prepare($sql);
	my @totals;
	if ($sth->execute()) {
		while (my @rset = $sth->fetchrow_array()) {
			push @totals, [@rset];
		}
    }
    return @totals;
} 

sub fetchIndividualUsage
{
	my ($self, $sql, $param1, $param2, $param3) = @_;
	my $vrtrack = $self->connectToDatabase($self->{VRTRACK_DATABASES}{WEB_TABLES});
	$self->displayError( "Failed to connect to vrtrack web database" ) unless defined( $vrtrack );	
	my $sth = $vrtrack->{_dbh}->prepare($sql);
	my @totals;
	if ( $param2 && $param3 ) {
		$sth->execute($param1, $param2, $param3);
	}
	else {
		$sth->execute($param1);
	}
	while (my (@rset) = $sth->fetchrow_array()) {
		push @totals, [@rset];
	}
    return @totals;
}

sub nonVrtrackConnection
{
	my ($self,$db) = @_;
    my $dbh = DBI->connect("dbi:mysql:$db;host=$ENV{VRTRACK_HOST};port=$ENV{VRTRACK_PORT}", $ENV{VRTRACK_RO_USER}, undef, { 'RaiseError' => 1 } );
    $self->displayError( "Failed to connect to non-VRTrack database: $db" ) unless defined( $dbh );
    return $dbh;
}

sub getHipsciCohorts
{
	my ($self, $db) = @_;
    my $vrtrack = $self->connectToDatabase($db);
	$self->displayError( "Failed to connect to database: $db" ) unless defined( $vrtrack );
	my $sql = $self->{SQL}{HIPSCI_FETCH_COHORT};
	my $sth = $vrtrack->{_dbh}->prepare($sql);
	my @cohorts;
	if ($sth->execute()) {
		while (my ($cohort) = $sth->fetchrow_array()) {
			push @cohorts, $cohort;
		}
    }
    return @cohorts;
}

sub getGenotypingSamples 
{
	my ($self, $db, $cohort) = @_;
	my %samples;
	my $vrtrack = $self->connectToDatabase($db);
	$self->displayError( "Failed to connect to database: $db" ) unless defined( $vrtrack );
	my $sql = $self->{SQL}{HIPSCI_FETCH_SAMPLES};
	my $sth = $vrtrack->{_dbh}->prepare($sql);	
	$sth->execute($cohort);
	while (my ($sample, $control) = $sth->fetchrow_array()) {
		$samples{$sample} = $control;
	}
	return %samples;	
}

sub getExpressionLaneID
{
	my ($self, $db, $cohort) = @_;
	my $vrtrack = $self->connectToDatabase($db);
	$self->displayError( "Failed to connect to database: $db" ) unless defined( $vrtrack );
	my $sql = $self->{SQL}{HIPSCI_FETCH_LANE_ID};
	my $sth = $vrtrack->{_dbh}->prepare($sql);	
	$sth->execute($cohort);
	my $lane_id;
	while (my ($id) = $sth->fetchrow_array()) {
		$lane_id = $id;
	}
    return $lane_id;	
}

sub getGenotypingCNVTotals 
{
	my ($self, $db, $cohort) = @_;
	my %cnv;
	my $vrtrack = $self->connectToDatabase($db);
	$self->displayError( "Failed to connect to database: $db" ) unless defined( $vrtrack );
	my $sql = $self->{SQL}{HIPSCI_FETCH_CNV_TOTALS};
	my $sth = $vrtrack->{_dbh}->prepare($sql);	
	$sth->execute($cohort);
	while (my ($sample, $intersection, $diff, $cnv_type) = $sth->fetchrow_array()) {
		push @{ $cnv{$sample}{$cnv_type} }, ($intersection, $diff);
	}
	return %cnv;	
}

sub getGenotypingCNVLeastValue
{
	my ($self, $db, $cohort, $cnv_flag) = @_;
	my %cnv;
	my $vrtrack = $self->connectToDatabase($db);
	$self->displayError( "Failed to connect to database: $db" ) unless defined( $vrtrack );
	my $sql = $self->{SQL}{HIPSCI_FETCH_CNV_TOTALS};
	my $sth = $vrtrack->{_dbh}->prepare($sql);	
	$sth->execute($cohort);
	my $least = -1;
	while (my ($sample, $intersection, $diff, $cnv_type) = $sth->fetchrow_array()) {
		$least = $diff if $cnv_type eq $cnv_flag && ( $least < 0 || $least > $diff);
	}
    return $least;	
}

1;
