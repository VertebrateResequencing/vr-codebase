#!/usr/local/bin/perl -T
# Displays projects view

BEGIN {
    $ENV{VRTRACK_HOST} = 'mcs4a';
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
use VRTrack::Project;
use VertRes::Utils::VRTrackFactory;
use VertRes::QCGrind::Util;

my $title = 'QC Grind Projects';
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
	my $st_samples_script = $utl->{SCRIPTS}{STUDY_SAMPLES};
    print qq[
    <h2 align="center" style="font: normal 900 1.5em arial"><a href="$main_script">QC Grind</a></h2>
    <h3 style="font: normal 700 1.5em arial">].ucfirst($database).qq[</h3>
    ];
    
    my $t = ucfirst( $database );
    print qq[
        <div class="centerFieldset">
        <fieldset style="width: 700px">
        <legend>Projects</legend>
        <table width="90%">
        <tr>
        <th>Name</th>
        <th>Accession</th>
        <th></th>
        </tr>
    ];
    
    foreach( @projects )
    {
        my $project = $_;
        my $study = $project->study();
        my $acc = '-';
        if( defined( $study ) ) {
            $acc = $study->acc();
        }
        my $pid = $project->id();
        my $name = $project->name;

        print qq[
            <tr>
                <td><a href="$st_lanes_script?db=$database&amp;proj_id=$pid">$name</a></td>
                <td>$acc</td>
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
