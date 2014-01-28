#!/usr/local/bin/perl -T

use strict;
use warnings;
use URI;

use lib '/var/www/lib';
use SangerPaths qw(core team145);
use SangerWeb;
use VertRes::QCGrind::ViewUtil;

my $utl = VertRes::QCGrind::ViewUtil->new();

my $title = 'Vrpipe Disk Usage';

my $sw  = SangerWeb->new({
    'title'   => $title,
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
    'jsfile'  => ['http://jsdev.sanger.ac.uk/prototype.js','http://jsdev.sanger.ac.uk/toggle.js','http://jsdev.sanger.ac.uk/scriptaculous/scriptaculous.js','http://jsdev.sanger.ac.uk/sidebar.js','http://jsdev.sanger.ac.uk/urchin.js','http://jsdev.sanger.ac.uk/zebra.js','http://js.sanger.ac.uk/sorttable_v2.js',],    
    'style'   => $utl->{CSS},
});

my $db = $utl->{VRTRACK_DATABASES}{WEB_TABLES};
my $cgi = $sw->cgi();

print $sw->header();
displayTopTenTotals($cgi, $title);
print $sw->footer();
exit;

#######

sub displayTopTenTotals
{
    my ($cgi,$title) = @_;
    
    my $index = $utl->{SCRIPTS}{INDEX_PAGE};
    
    my $directory_script = $utl->{SCRIPTS}{DISK_USAGE_DIRECTORY};
    my $indiv_dir_script = $utl->{SCRIPTS}{DISK_USAGE_INDIV_DIR};
        
    my $pipeline_script = $utl->{SCRIPTS}{DISK_USAGE_PIPELINE};
    my $indiv_pipe_script = $utl->{SCRIPTS}{DISK_USAGE_INDIV_PIPE};

		
	print qq[ <h4 align="center" style="font: arial"><i><a href="$index">Team 145</a></i> : $title</h4> ];
	
    my @directories = $utl->fetchUsageTotals($utl->{SQL}{DISK_USAGE_DIRECTORY});
    my @projects = $utl->fetchUsageTotals($utl->{SQL}{DISK_USAGE_PIPELINE});
     
    print qq[<div class="centerFieldset">
    <fieldset style = "border:none">
    <table width="80%">];
    
    print qq[
    <div class="centerFieldset">
    <fieldset > 
    <legend>Disk usage - top 10 directories <h1><a href="$directory_script">Show All Directories</a></h1></legend>
    <table class='sortable' width="80%">
    <tr>
    <th align="center">Directory root</th>
    <th align="center">File types</th>
    <th align="center">Disk usage (GB)</th>
    <th align="center">No. of files</th>
    </tr>
    ];
    
    foreach my $pipeline(@directories[0 .. 9]) {
		my @pipeline_data = @{ $pipeline };
		my $pslink = "<a href='$indiv_dir_script?root=$pipeline_data[0]'> $pipeline_data[0] </a>";
		print qq[<tr><td align="center">$pslink</td><td align="center">$pipeline_data[1]</td><td align="center">$pipeline_data[3]</td><td align="center">$pipeline_data[4]</td></tr>];
	}
    print qq[
        </table>
        </fieldset>
        </div>
    ];
    
    print qq[<div class="centerFieldset">
    <fieldset style = "border:none">
    <table width="80%">];
    
    print qq[
    <div class="centerFieldset">
    <fieldset > 
    <legend>Disk usage - top 10 projects <h1><a href="$pipeline_script">Show All Pipelines</a></h1></legend>
    <table class='sortable' width="80%">
    <tr>
    <th>Pipeline ID</th>
    <th>Pipeline name</th>
    <th>Pipeline type</th>
    <th>Disk usage (GB)</th>
    <th align="center">No. of files</th>
    <th align="center">Owner</th>
    </tr>
    ];
    
    foreach my $pipeline(@projects[0 .. 9]) {
		my @pipeline_data = @{ $pipeline };
		my $pslink = "<a href='$indiv_pipe_script?setup=$pipeline_data[0]'> $pipeline_data[1] </a>";
		print qq[<tr><td align="center">$pipeline_data[0]</td><td align="center">$pslink</td><td align="center">$pipeline_data[3]</td><td align="center">$pipeline_data[4]</td><td align="center">$pipeline_data[5]</td><td align="center">$pipeline_data[2]</td></tr>];
	}
    print qq[
        </table>
        </fieldset>
        </div>
    ];
    
}
