#!/usr/local/bin/perl -T

# Displays QC information for Sanger short-read sequencing to allow lanes
# to be passed/failed

use strict;
use warnings;
use URI;

use lib '/var/www/lib';
use SangerPaths qw(core team145);
use SangerWeb;
use VertRes::Utils::VRTrackFactory;
use VertRes::QCGrind::ViewUtil;

my $utl = VertRes::QCGrind::ViewUtil->new();

my $title = 'Map View';

my $sw  = SangerWeb->new({
    'title'   => $title,
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
    'jsfile'  => ['http://jsdev.sanger.ac.uk/prototype.js','http://jsdev.sanger.ac.uk/jquery-1.4.2.min.js','http://www.sanger.ac.uk/modelorgs/mousegenomes/jquery.coolfieldset.js'],
    'style'   => $utl->{CSS},
});

my $cgi = $sw->cgi();

my $script = $utl->{SCRIPTS}{MAP_PROJECTS_VIEW};

print $sw->header();
$utl->displayDatabasesPage($title,$cgi,$script,1,0);
print $sw->footer();
exit;

