#!/usr/local/bin/perl -T

BEGIN {
    $ENV{VRTRACK_HOST} = 'mcs10';
    $ENV{VRTRACK_PORT} = 3306;
    $ENV{VRTRACK_RO_USER} = 'vreseq_ro';
    $ENV{VRTRACK_RW_USER} = 'vreseq_rw';
    $ENV{VRTRACK_PASSWORD} = 't3aml3ss';
};

use strict;
use warnings;
use URI;

use SangerPaths qw(core team145);
use SangerWeb;
use VertRes::QCGrind::ViewUtil;

my $utl = VertRes::QCGrind::ViewUtil->new();

my $title = 'Pipeline Disk Usage';

my $sw  = SangerWeb->new({
    'title'   => $title,
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
    'jsfile'  => ['http://jsdev.sanger.ac.uk/prototype.js','http://jsdev.sanger.ac.uk/toggle.js','ttp://jsdev.sanger.ac.uk/scriptaculous/scriptaculous.js','http://jsdev.sanger.ac.uk/sidebar.js','http://jsdev.sanger.ac.uk/urchin.js','http://jsdev.sanger.ac.uk/zebra.js','http://js.sanger.ac.uk/sorttable_v2.js',],    
    'style'   => $utl->{CSS},
});

my $db = 'vrtrack_web_index';
my $cgi = $sw->cgi();

if ($cgi->param('download')) {
	my $sql = $utl->{SQL}{DISK_USAGE_PIPELINE};
    if ($cgi->param('download') eq 'TSV') {
    	&downloadMappings($cgi, $sql, "\t", 'tsv');
    	exit;
    }
    elsif ($cgi->param('download') eq 'CSV (Excel)') {
    	&downloadMappings($cgi, $sql, ',', 'csv');
    	exit;
    }	
}

print $sw->header();
displayPipelineTotals($cgi, $title);
print $sw->footer();
exit;

#######

sub displayPipelineTotals
{
    my ($cgi,$title) = @_;
            
	my $index = $utl->{SCRIPTS}{INDEX_PAGE};
    my $init_script = $utl->{SCRIPTS}{DISK_USAGE};
    my $pipeline_script = $utl->{SCRIPTS}{DISK_USAGE_INDIV_PIPE};
    	
	print qq[ <h4 align="center" style="font: arial"><i><a href="$index">Team 145</a> : <a href="$init_script">Vrpipe Disk Usage</a></i> : $title</h4> ];
    my @diskusagetotals = $utl->fetchUsageTotals($utl->{SQL}{DISK_USAGE_PIPELINE});
    
    print qq[<div class="centerFieldset">
    <fieldset style = "border:none">
    <table width="80%">];
    
    print $cgi->start_form;
    print qq[
    <tr></tr>
    <tr>
    <td colspan = 4 align='right'>Download:];
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
    <legend>Disk usage totals</legend>
    <table class='sortable' width="80%">
    <tr>
    <th>Pipeline ID</th>
    <th>Pipeline name</th>
    <th>Pipeline type</th>
    <th>Disk usage (GB)</th>
    <th>No. of files</th>
    <th>Owner</th>
    </tr>
    ];
    
    my $href = "<a href='pipeline_script?setup=";
    foreach my $pipeline(@diskusagetotals) {
		my @pipeline_data = @{ $pipeline };
		my $pslink = "<a href='$pipeline_script?setup=$pipeline_data[0]'> $pipeline_data[1] </a>";
		print qq[<tr><td align="center">$pipeline_data[0]</td><td align="center">$pslink</td><td align="center">$pipeline_data[3]</td><td align="center">$pipeline_data[4]</td><td align="center">$pipeline_data[5]</td><td align="center">$pipeline_data[2]</td></tr>];
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
	my ($cgi, $sql, $sep, $type) = @_;
	print $cgi->header(-type=>"text/$type",  -attachment=>"disk_usage_vrpipe_pipeline_totals.$type");
	print join ($sep,"Pipeline ID", "Pipeline name", "Pipeline owner", "Pipeline type", "Disk usage (GB)", "No. of files"), "\n";
	my @usage = $utl->fetchUsageTotals($sql);
	foreach my $pipeline(@usage) {
		my @pipeline_data = @{ $pipeline };
		$pipeline_data[1] =~ s/,/;/g ;
		print join ($sep, @pipeline_data), "\n";
	}
}
