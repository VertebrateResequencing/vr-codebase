#!/usr/local/bin/perl -T
# Display QC Grind projects view

use strict;
use warnings;
use URI;

use lib '/var/www/lib';
use SangerPaths qw(core team145);
use SangerWeb;
use VRTrack::Project;
use VertRes::Utils::VRTrackFactory;
use VertRes::QCGrind::Util;

my $title = 'Projects View';
my $sw  = SangerWeb->new({
    'title'   => $title,
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
});
my $utl = VertRes::QCGrind::Util->new();

my $cgi = $sw->cgi();

my $db = $cgi->param('db');
unless ($db) {
    print $sw->header();
	$utl->displayError( "No database ID",$sw );
}
my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $db, mode => 'rw');

print $sw->header();
displayProjectsPage($cgi, $vrtrack, $db);
print $sw->footer();

exit;

#####
sub displayProjectsPage
{
    my ($cgi, $vrtrack, $database ) = @_;
    
    my @projects = sort {$a->name cmp $b->name} @{$vrtrack->projects()};
    
	my $main_script = $utl->{SCRIPTS}{DATABASES_VIEW};
	my $st_lanes_script = $utl->{SCRIPTS}{STUDY_LANES};
	my $st_samples_script = $utl->{SCRIPTS}{SAMPLES_VIEW};
    print qq[
    <h2 align="center" style="font: normal 900 1.5em arial"><a href="$main_script">QC Grind</a> $title</h2>
    <h3 align="center" style="font: normal 700 1.5em arial">Database : $database</h3>
    ];
    
    print qq[
        <div class="centerFieldset">
        <fieldset style="width: 700px">
        <legend>Projects</legend>
        <table id="projects" class="zebra" width="100%" cellpadding=4>
        <tr>
        <th>Name</th>
        <th>Accession</th>
        <th></th>
        <th></th>
        <th></th>
        </tr>
    ];
    
    foreach my $project ( @projects ) {

        my $study = $project->study();
        my $acc = '-';
        if( defined( $study ) ) {
            $acc = $study->acc();
        }
        my $pid = $project->id();
        my $projname = $project->name;

        print qq[
            <tr>
                <td>$projname</td>
                <td>$acc</td>
                <td><a href="study.pl?db=$database&amp;proj_id=$pid">Lanes</a></td>
                <td><a href="$st_lanes_script?db=$database&amp;proj_id=$pid">Lanelets</a></td>
                <td><a href="$st_samples_script?db=$database&amp;proj_id=$pid">Samples</a></td>
            </tr>
        ];
    }
    print qq[
        </table>
        </fieldset>
        </div>
    ];
}
