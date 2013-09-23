#!/usr/local/bin/perl -T

use strict;
use warnings;
use URI;

use lib '/var/www/lib';
use SangerPaths qw(core team145);
use SangerWeb;
use VertRes::QCGrind::ViewUtil;

my $title = "Files Disk Usage";

my $utl = VertRes::QCGrind::ViewUtil->new();

my $sw  = SangerWeb->new({
    'title'   => $title,
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
    'jsfile'  => ['http://jsdev.sanger.ac.uk/prototype.js','http://jsdev.sanger.ac.uk/toggle.js','http://jsdev.sanger.ac.uk/scriptaculous/scriptaculous.js','http://jsdev.sanger.ac.uk/sidebar.js','http://jsdev.sanger.ac.uk/urchin.js','http://jsdev.sanger.ac.uk/zebra.js','http://js.sanger.ac.uk/sorttable_v2.js',],    
    'style'   => $utl->{CSS},
});

my $cgi = $sw->cgi();

my $ps_id = $cgi->param('setup');
my $root = $cgi->param('root');
my $type = $cgi->param('type');
my $refer = $cgi->param('refer');

unless ($ps_id && $root && $type && $refer) {
    print $sw->header();
	$utl->displayError( "Incomplete information supplied",$sw );
}


if ($cgi->param('download')) {
	my $sql = $utl->{SQL}{DISK_USAGE_FILES};
    if ($cgi->param('download') eq 'TSV') {
    	&downloadMappings($cgi, $sql, "\t", 'tsv', $ps_id, $root, $type);
    	exit;
    }
    elsif ($cgi->param('download') eq 'CSV (Excel)') {
    	&downloadMappings($cgi, $sql, ',', 'csv', $ps_id, $root, $type);
    	exit;
    }	
}

print $sw->header();
displayFileView($cgi, $title, $ps_id, $root, $type, $refer);
print $sw->footer();
exit;

#######

sub displayFileView
{
    my ($cgi,$title,$ps_id,$root,$type,$refer) = @_;
    
    my $index = $utl->{SCRIPTS}{INDEX_PAGE};
    my $init_script = $utl->{SCRIPTS}{DISK_USAGE};
       
    my $directory_script = $refer eq 'pipe' ? $utl->{SCRIPTS}{DISK_USAGE_PIPELINE} : $utl->{SCRIPTS}{DISK_USAGE_DIRECTORY};
    my $individual_script = $refer eq 'pipe' ? $utl->{SCRIPTS}{DISK_USAGE_INDIV_PIPE} : $utl->{SCRIPTS}{DISK_USAGE_INDIV_DIR};
    my $second_level = $refer eq 'pipe' ? 'Pipeline Disk Usage' : 'Directory Disk Usage';
    my $init = $refer eq 'pipe' ? 'Setup $ps_id' : '';
    my $third_level = $refer eq 'pipe' ? "setup=$ps_id\">Setup $ps_id" : "root=$root\">$root";
    my $fourth_level = $refer eq 'pipe' ? "$root ($type files)" : "setup $ps_id";
    
    print qq[
        <h4 align="center" style="font: arial"><i><a href="$index">Team 145</a> : <a href="$init_script">Vrpipe Disk Usage</a> : <a href="$directory_script">$second_level</a> : <a href="$individual_script?$third_level</a></i> : $fourth_level</h4>
		<div class="centerFieldset">
    ];
	my @fileusage = $utl->fetchIndividualUsage($utl->{SQL}{DISK_USAGE_FILES}, $ps_id, $root, $type);
	my $file_count = scalar @fileusage;
    print qq[<div class="centerFieldset">
    <fieldset style = "border:none">
    <table width="80%">];
    
    print $cgi->start_form;
    print $cgi->hidden("setup","$ps_id");
    print $cgi->hidden("root","$root");
    print $cgi->hidden("type","$type");
    print $cgi->hidden("refer","$refer");
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
    <fieldset > ];
    
    if ( $file_count > 50 ) {
		print qq[<legend>Pipeline disk usage : Top 50 (out of $file_count) largest $type files</legend>];
	}
	else {
		print qq[<legend>Pipeline disk usage</legend>];
	}
    my $limit = $file_count > 50 ? 49 : $file_count-1;
    #need to check end file size - if reads 0GB, output in bytes
    my $check_gb_size = int ( @{ $fileusage[$limit] }[3] / 1073741824);
    my $size_units = $check_gb_size > 0 ? 'GB' : 'Bytes';
	print qq[<table class='sortable' width="80%">
    <tr>
    <th>Pipeline ID</th>
    <th>Owner</th>
    <th>File path</th>
    <th>File type</th>
    <th>Disk usage ($size_units)</th>
    </tr>
    ];

    foreach my $pipeline(@fileusage[0 .. $limit]) {
		my @pipeline_data = @{ $pipeline };
		my $usage = $check_gb_size > 0 ? int($pipeline_data[4] / 1073741824) : $pipeline_data[4];
		print qq[<tr><td align="center">$pipeline_data[0]</td><td align="center">$pipeline_data[2]</td><td align="center">$pipeline_data[6]</td><td align="center">$pipeline_data[5]</td><td align="center">$usage</td></tr>];
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
	my ($cgi, $sql, $sep, $type, $ps_id, $root, $root_file_type) = @_;
	print $cgi->header(-type=>"text/$type",  -attachment=>"disk_usage_pipeline_$ps_id\_$root_file_type\_files.$type");
	print join ($sep, "Pipeline ID", "Pipeline name", "Pipeline owner", "File path", "File type", "Disk usage (Bytes))", "\n");
	my @usage = $utl->fetchIndividualUsage($sql, $ps_id, $root, $root_file_type);
	foreach my $pipeline(@usage) {
		my @pipeline_data = @{ $pipeline };
		$pipeline_data[1] =~ s/,//g ;
		#my $usage = int($pipeline_data[3] / 1073741824);
		print (join $sep, $ps_id, $pipeline_data[1], $pipeline_data[2], $pipeline_data[5], $pipeline_data[4], $pipeline_data[3], "\n");
	}
}
