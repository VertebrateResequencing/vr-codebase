#!/usr/local/bin/perl -T

use strict;
use warnings;
use URI;

use lib '/var/www/lib';
use SangerPaths qw(core team145);
use SangerWeb;
use VertRes::QCGrind::ViewUtil;

my $title = "Individual Pipeline Usage";

my $utl = VertRes::QCGrind::ViewUtil->new();

my $sw  = SangerWeb->new({
    'title'   => $title,
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
    'jsfile'  => ['/js.sanger.ac.uk/prototype.js','/js.sanger.ac.uk/toggle.js','/js.sanger.ac.uk/scriptaculous/scriptaculous.js','/js.sanger.ac.uk/sidebar.js','/js.sanger.ac.uk/urchin.js','/js.sanger.ac.uk/zebra.js','/js.sanger.ac.uk/sorttable_v2.js',],    
    'style'   => $utl->{CSS},
});

my $cgi = $sw->cgi();

my $ps_id = $cgi->param('setup');

unless ($ps_id) {
    print $sw->header();
	$utl->displayError( "No pipeline setup ID",$sw );
}

my $db = $utl->{VRTRACK_DATABASES}{WEB_TABLES};

if ($cgi->param('download')) {
	my $sql = $utl->{SQL}{DISK_USAGE_INDIV_PIPE};
    if ($cgi->param('download') eq 'TSV') {
    	&downloadMappings($cgi, $sql, "\t", 'tsv', $ps_id);
    	exit;
    }
    elsif ($cgi->param('download') eq 'CSV (Excel)') {
    	&downloadMappings($cgi, $sql, ',', 'csv', $ps_id);
    	exit;
    }	
}

print $sw->header();
displayPipelineOverview($cgi, $title, $ps_id);
print $sw->footer();
exit;

#######

sub displayPipelineOverview
{
    my ($cgi,$title,$ps_id) = @_;
    
    my $index = $utl->{SCRIPTS}{INDEX_PAGE};
    my $init_script = $utl->{SCRIPTS}{DISK_USAGE};
    my $file_script = $utl->{SCRIPTS}{DISK_USAGE_FILES};
    my $pipeline_script = $utl->{SCRIPTS}{DISK_USAGE_PIPELINE};
    
	print qq[
        <h4 align="center" style="font: arial"><i><a href="$index">Team 145</a> : <a href="$init_script">Vrpipe Disk Usage</a> : <a href="$pipeline_script">Pipeline Disk Usage</a></i> : Setup $ps_id</h4>
		<div class="centerFieldset">
    ];
	my @pipelineusage = $utl->fetchIndividualUsage($utl->{SQL}{DISK_USAGE_INDIV_PIPE}, $ps_id, 0, 0);
    
    print qq[<div class="centerFieldset">
    <fieldset style = "border:none">
    <table width="80%">];
    
    print $cgi->start_form;
    print $cgi->hidden("setup","$ps_id");
    print qq[
    <tr></tr>
    <tr>
    <td colspan = 6 align='right'>Download:];
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
    <legend>Pipeline disk usage overview</legend>
    <table class='sortable' width="80%">
    <tr>
    <th>Pipeline ID</th>
    <th>Pipeline name</th>
    <th>Pipeline type</th>
    <th>File type</th>
    <th>Path root</th>
    <th>Disk usage (GB)</th>
    <th>No. of files</th>
    <th>Owner</th>
    </tr>
    ];
    
    foreach my $pipeline(@pipelineusage) {
		my @pipeline_data = @{ $pipeline };
        my $filelink = "<a href='$file_script?setup=$pipeline_data[0]\&root=$pipeline_data[5]\&type=$pipeline_data[4]\&refer=pipe'> $pipeline_data[6] </a>";
		print qq[<tr><td align="center">$pipeline_data[0]</td><td align="center">$pipeline_data[1]</td><td align="center">$pipeline_data[3]</td><td align="center">$pipeline_data[4]</td><td align="center">$pipeline_data[5]</td><td align="center">$filelink</td><td align="center">$pipeline_data[7]</td><td align="center">$pipeline_data[2]</td></tr>];
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
	my ($cgi, $sql, $sep, $type, $ps_id) = @_;
	print $cgi->header(-type=>"text/$type",  -attachment=>"disk_usage_overview_setup_$ps_id.$type");
	print join ($sep,"Pipeline ID", "Pipeline name", "Pipeline owner", "Pipeline type", "File type", "Path root", "Disk usage (GB)", "No. of files"), "\n";
	my @pipeline_usage = $utl->fetchIndividualUsage($sql, $ps_id, 0, 0);
	foreach my $pipeline(@pipeline_usage) {
		my @pipeline_data = @{ $pipeline };
		$pipeline_data[1] =~ s/,//g ;
		print (join $sep, @pipeline_data, "\n");
    }	
}
