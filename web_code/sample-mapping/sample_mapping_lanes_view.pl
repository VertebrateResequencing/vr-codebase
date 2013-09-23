#!/usr/local/bin/perl -T

use strict;
use warnings;
use lib '../../lib';
use URI;

use SangerPaths qw(core team145);
use SangerWeb;
use VRTrack::VRTrack;
use VRTrack::Project;
use VertRes::Utils::VRTrackFactory;
use VertRes::QCGrind::ViewUtil;

my $utl = VertRes::QCGrind::ViewUtil->new();

my $title = 'Sample ID Mapping';

my $sw  = SangerWeb->new({
    'title'   => $title,
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
    'jsfile'  => ['http://jsdev.sanger.ac.uk/prototype.js','http://jsdev.sanger.ac.uk/toggle.js','ttp://jsdev.sanger.ac.uk/scriptaculous/scriptaculous.js','http://jsdev.sanger.ac.uk/sidebar.js','http://jsdev.sanger.ac.uk/urchin.js','http://jsdev.sanger.ac.uk/zebra.js','http://js.sanger.ac.uk/sorttable_v2.js',],    
    'style'   => $utl->{CSS},
});

my $cgi = $sw->cgi();

my $db = $cgi->param('db');
my $pid = $cgi->param('proj_id');

unless ($db) {
    print $sw->header();
	$utl->displayError( "No database ID",$sw );
}
my $vrtrack = $utl->connectToDatabase($db);

if ($cgi->param('download')) {
	my $pid = $cgi->param('proj_id');
    if( defined( $pid ) ) {
        if ($cgi->param('download') eq 'TSV') {
    		&downloadMappings($cgi, "\t", 'tsv', $db, $pid);
    		exit;
    	}
    	elsif ($cgi->param('download') eq 'CSV (Excel)') {
    		&downloadMappings($cgi, "\t", 'tsv', $db, $pid);
    		exit;
    	}
    }	
}

$utl->displayError( "No Project ID",$sw ) unless $pid;

print $sw->header();
displayProjectMappingsPage($cgi,$db,$pid);
print $sw->footer();
exit;

#######

sub displayProjectMappingsPage 
{
    my ($cgi,$database, $projectID) = @_;
    
    my $init_script = $utl->{SCRIPTS}{SAMPLE_MAPPING};
    my $proj_script = $utl->{SCRIPTS}{SAMPMAP_PROJECTS_VIEW};
    
    my $db_id = $utl->getDatabaseID ($database);
     
    my $pname = $projectID == 99999 ? "All $database projects" : $utl->fetchProjectName($db_id, $projectID);
    
    my $index = $utl->{SCRIPTS}{INDEX_PAGE};
    
    print qq[
    <h4 align="center" style="font: arial"><i><a href="$index">Team 145</a> : <a href="$init_script">$title</a> :  <a href="$proj_script?db=$database">$database</a></i> : $pname</h4>
    ];
    
    my %sampleMappings = $utl->getSampleMappings ($db_id, $projectID);
    
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

sub downloadMappings
{
	my ($cgi, $sep, $type, $db, $projectID) = @_;
	my $db_id = $utl->getDatabaseID ($db);
	my $pname = $projectID == 99999 ? "all_$db" :$utl->fetchProjectName($db_id, $projectID);
  my $vrtrack = $utl->connectToDatabase($utl->{VRTRACK_DATABASES}{WEB_TABLES});
  print $cgi->header(-type=>'text/$type',  -attachment=>"$pname.$type");
	print join ($sep,"Sanger Sample Name","Supplier Name","Accession"), "\n";
	my $sql = qq[SELECT supplier_name, accession_number, sanger_sample_name FROM sample_id_mapping where db_id = ?];
	$sql = $sql." and project_id = ?" unless $projectID == 99999;
	my $sth = $vrtrack->{_dbh}->prepare($sql);	
	if ( $projectID == 99999 ) {
		$sth->execute($db_id);
	}
	else {
		$sth->execute($db_id, $projectID);
	}
	while (my ($supp, $acc, $sang) = $sth->fetchrow_array()) {
		print (join ($sep, $sang, $supp, $acc), "\n"); 
	}
}

#~ sub downloadMappingsCSV
#~ {
	#~ my ($cgi, $db, $projectID) = @_;
	#~ my $db_id = $utl->getDatabaseID ($db);
	#~ my $pname = $utl->fetchProjectName($db_id, $projectID);
	#~ my $vrtrack = $utl->connectToDatabase('vrtrack_web_index');
	#~ print $cgi->header(-type=>'text/csv',  -attachment=>"$pname.csv");
	#~ print join (",","Sanger Sample Name","Supplier Name","Accession"), "\n";
	#~ my $sql = qq[SELECT supplier_name, accession_number, sanger_sample_name FROM sample_id_mapping where db_id = ? and project_id = ?];
	#~ my $sth = $vrtrack->{_dbh}->prepare($sql);
	#~ if ($sth->execute($db_id, $projectID)) {
		#~ my ($supp, $acc, $sang);
		#~ $sth->bind_columns(\($supp, $acc, $sang));
		#~ while ($sth->fetch) {
		    #~ print (join (",",$sang, $supp, $acc), "\n"); 
		#~ }
	#~ }
#~ }
