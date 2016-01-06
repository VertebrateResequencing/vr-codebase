#!/usr/local/bin/perl -T

use strict;
use warnings;
use URI;

use lib '/var/www/lib';
use SangerPaths qw(core team145);
use SangerWeb;
use VertRes::QCGrind::ViewUtil;

my $title = "Individual Directory Usage";

my $utl = VertRes::QCGrind::ViewUtil->new();

my $sw  = SangerWeb->new({
    'title'   => $title,
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
    'jsfile'  => ['/js.sanger.ac.uk/prototype.js','/js.sanger.ac.uk/toggle.js','/js.sanger.ac.uk/scriptaculous/scriptaculous.js','/js.sanger.ac.uk/sidebar.js','/js.sanger.ac.uk/urchin.js','/js.sanger.ac.uk/zebra.js','/js.sanger.ac.uk/sorttable_v2.js',],    
    'style'   => $utl->{CSS},
});

my $cgi = $sw->cgi();

my $root = $cgi->param('root');

unless ($root) {
    print $sw->header();
	$utl->displayError( "Incomplete info",$sw );
}

if ($cgi->param('download')) {
	my $sql = $utl->{SQL}{DISK_USAGE_INDIV_DIR};
    if ($cgi->param('download') eq 'TSV') {
    	&downloadMappings($cgi, $sql, "\t", 'tsv', $root);
    	exit;
    }
    elsif ($cgi->param('download') eq 'CSV (Excel)') {
    	&downloadMappings($cgi, $sql, ',', 'csv', $root);
    	exit;
    }	
}

print $sw->header();
displayIndividualDirectory($cgi, $title, $root);
print $sw->footer();
exit;

#######

sub displayIndividualDirectory
{
    my ($cgi,$title,$root) = @_;
    
    my $index = $utl->{SCRIPTS}{INDEX_PAGE};
    my $init_script = $utl->{SCRIPTS}{DISK_USAGE};
    my $directory_script = $utl->{SCRIPTS}{DISK_USAGE_DIRECTORY};
    my $file_script = $utl->{SCRIPTS}{DISK_USAGE_FILES};
    
	print qq[
        <h4 align="center" style="font: arial"><i><a href="$index">Team 145</a> : <a href="$init_script">Vrpipe Disk Usage</a> : <a href="$directory_script">Directory Disk Usage</a></i> : $root</h4>
		<div class="centerFieldset">
    ];
	my @directory_usage = $utl->fetchIndividualUsage($utl->{SQL}{DISK_USAGE_INDIV_DIR}, $root, 0, 0);
    print qq[<div class="centerFieldset">
    <fieldset style = "border:none">
    <table width="80%">];
    
    print $cgi->start_form;
    print $cgi->hidden("root","$root");
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
    <legend>Directory disk usage overview</legend>
    <table class='sortable' width="80%">
    <tr>
    <th align="center">Pipeline ID</th>
    <th align="center">Pipeline name</th>
    <th align="center">Pipeline type</th>
    <th align="center">File type</th>
    <th align="center">Directory</th>
    <th align="center">Disk usage (GB)</th>
    <th align="center">No. of files</th>
    <th align="center">Owner</th>
    </tr>
    ];
    
    foreach my $pipeline(@directory_usage) {
		my @pipeline_data = @{ $pipeline };
        my $filelink = "<a href='$file_script?setup=$pipeline_data[0]\&root=$pipeline_data[5]\&type=$pipeline_data[4]\&refer=dir'> $pipeline_data[1] </a>";
		print qq[<tr><td align="center">$pipeline_data[0]</td><td align="center">$filelink</td><td align="center">$pipeline_data[3]</td><td align="center">$pipeline_data[4]</td><td align="center">$pipeline_data[5]</td><td align="center">$pipeline_data[6]</td><td align="center">$pipeline_data[7]</td><td align="center">$pipeline_data[2]</td></tr>];
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
	my ($cgi, $sql, $sep, $type, $root) = @_;
	my $file_name = $root;
	$file_name =~ s/\//_/g ;
	my @directory_usage = $utl->fetchIndividualUsage($sql, $root, 0, 0);
	print $cgi->header(-type=>"text/$type",  -attachment=>"disk_usage_overview$file_name.$type");
	print join ($sep, "Pipeline ID", "Pipeline name", "Pipeline owner",  "Pipeline type", "File type", "Directory", "Disk usage (GB))", "\n");
	foreach my $pipeline(@directory_usage) {
		my @pipeline_data = @{ $pipeline };
		$pipeline_data[1] =~ s/,//g ;
		print (join $sep, @pipeline_data, "\n");
    }
}
