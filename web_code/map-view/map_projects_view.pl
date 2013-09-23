#!/usr/local/bin/perl -T

use strict;
use warnings;
use URI;

use lib '../../lib';
use SangerPaths qw(core team145);
use SangerWeb;
use VRTrack::VRTrack;
use VRTrack::Project;
use VertRes::Utils::VRTrackFactory;
use VertRes::QCGrind::ViewUtil;

my $utl = VertRes::QCGrind::ViewUtil->new();

my $title = 'Map View';

my $sw  = SangerWeb->new({
    'title'   => $title,
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
});

my $cgi = $sw->cgi();

my $db = $cgi->param('db');
unless ($db) {
    print $sw->header();
	$utl->displayError( "No database ID",$sw );
}
my $vrtrack = $utl->connectToDatabase($db);
unless ( defined $vrtrack ) {
	print $sw->header();
	$utl->displayError( "No database connection to $db",$sw );
}

my $init_script = $utl->{SCRIPTS}{MAP_VIEW};
my $lanes_script = $utl->{SCRIPTS}{MAP_LANES_VIEW};

print $sw->header();
$utl->displayDatabasePage($title,$cgi,$vrtrack,$db,$init_script,$lanes_script,0);
print $sw->footer();
exit;
