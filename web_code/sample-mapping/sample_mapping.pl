#!/usr/local/bin/perl -T

use strict;
use warnings;
use lib '/var/www/lib';
use URI;

use SangerPaths qw(core team145);
use SangerWeb;
use VertRes::Utils::VRTrackFactory;
use VertRes::QCGrind::ViewUtil;

my $utl = VertRes::QCGrind::ViewUtil->new();

my $title = 'Sample ID Mapping';

my $sw  = SangerWeb->new({
    'title'   => $title,
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
    'jsfile'  => ['/js.sanger.ac.uk/prototype.js','/js.sanger.ac.uk/jquery-1.4.2.min.js','http://www.sanger.ac.uk/modelorgs/mousegenomes/jquery.coolfieldset.js'],
    'style'   => $utl->{CSS},
});

my $cgi = $sw->cgi();

my $script = $utl->{SCRIPTS}{SAMPMAP_PROJECTS_VIEW};

print $sw->header();
$utl->displayDatabasesPage($title,$cgi,$script,1,0);
print $sw->footer();
exit;
