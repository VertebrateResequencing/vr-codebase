#!/usr/local/bin/perl -T

use strict;
use warnings;
use URI;

use lib '../../lib';
use SangerPaths qw(core team145);
use SangerWeb;
use VertRes::QCGrind::ViewUtil;

my $utl = VertRes::QCGrind::ViewUtil->new();

my $title = 'Directory Disk Usage';

my $sw  = SangerWeb->new({
    'title'   => $title,
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
    'jsfile'  => ['http://jsdev.sanger.ac.uk/prototype.js','http://jsdev.sanger.ac.uk/toggle.js','http://jsdev.sanger.ac.uk/scriptaculous/scriptaculous.js','http://jsdev.sanger.ac.uk/sidebar.js','http://jsdev.sanger.ac.uk/urchin.js','http://jsdev.sanger.ac.uk/zebra.js','http://js.sanger.ac.uk/sorttable_v2.js',],    
    'style'   => $utl->{CSS},
});

my $cgi = $sw->cgi();

if ($cgi->param('download')) {
	my $sql = $utl->{SQL}{DISK_USAGE_DIRECTORY};
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
displayRootTotals($cgi, $title);
print $sw->footer();
exit;

#######

sub displayRootTotals
{
    my ($cgi,$title) = @_;
    
 	my $index = $utl->{SCRIPTS}{INDEX_PAGE};
 
    my $init_script = $utl->{SCRIPTS}{DISK_USAGE};
    my $directory_script = $utl->{SCRIPTS}{DISK_USAGE_INDIV_DIR};
        
    my @rootusagetotals = $utl->fetchUsageTotals($utl->{SQL}{DISK_USAGE_DIRECTORY});
		
	print qq[ <h4 align="center" style="font: arial"><i><a href="$index">Team 145</a> : <a href="$init_script">Vrpipe Disk Usage</a></i> : $title</h4> ];

    print qq[<div class="centerFieldset">
    <fieldset style = "border:none">
    <table width="80%">];
    
    print $cgi->start_form;

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
    <legend>Disk usage totals</legend>
    <table class='sortable' width="80%">
    <tr>
    <th align="center">Directory</th>
    <th align="center">File types</th>
    <th align="center">Disk usage (GB)</th>
    <th align="center">No. of files</th>
    </tr>
    ];
    
    foreach my $pipeline(@rootusagetotals) {
		my @pipeline_data = @{ $pipeline };
		my $pslink = "<a href='$directory_script?root=$pipeline_data[0]'> $pipeline_data[0] </a>";
		print qq[<tr><td align="center">$pslink</td><td align="center">$pipeline_data[1]</td><td align="center">$pipeline_data[3]</td><td align="center">$pipeline_data[4]</td></tr>];
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
	print $cgi->header(-type=>"text/$type",  -attachment=>"directory_disk_usage.$type");
	print join ($sep,"Directory","File types","Disk usage (GB)","No. of files"), "\n";
	my @usage = $utl->fetchUsageTotals($sql);
	foreach my $pipeline(@usage) {
		my @pipeline_data = @{ $pipeline };
		$pipeline_data[1] =~ s/,/;/g ;
		print join ($sep, $pipeline_data[0], $pipeline_data[1], $pipeline_data[3], $pipeline_data[4]), "\n";
	}
}
