#!/usr/local/bin/perl -T

use strict;
use warnings;
use URI;

use lib '/var/www/lib';
use SangerPaths qw(core team145);
use SangerWeb;
use VertRes::Utils::VRTrackFactory;
use VertRes::QCGrind::ViewUtil;

my $utl = VertRes::QCGrind::ViewUtil->new();

my $title = 'Pending View';

my $sw  = SangerWeb->new({
    'title'   => $title,
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
    'jsfile'  => ['/js.sanger.ac.uk/prototype.js','/js.sanger.ac.uk/jquery-1.4.2.min.js','http://www.sanger.ac.uk/modelorgs/mousegenomes/jquery.coolfieldset.js'],
    'style'   => $utl->{CSS},
});

my $cgi = $sw->cgi();

my $script = $utl->{SCRIPTS}{PENDING_REQ_VIEW};

print $sw->header();
$utl->displayDatabasesPage($title,$cgi,$script,1,1);
print $sw->footer();
exit;
