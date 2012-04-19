#!/usr/local/bin/perl -T
# Displays vtrack databases

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

my $title = 'QC Grind Databases';
my $sw  = SangerWeb->new({
    'title'   => $title,
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
});
my $utl = VertRes::QCGrind::Util->new();

my $cgi = $sw->cgi();

print $sw->header();
displayDatabasesPage();
print $sw->footer();

exit;

#####

sub displayDatabasesPage {

    print qq[
        <h2 align="center" style="font: normal 900 1.5em arial">$title</h2>
        <div class="centerFieldset">
        <fieldset style="width: 500px">
        <legend>Databases</legend>
    ];
    
    my @dbs = VertRes::Utils::VRTrackFactory->databases(1);

	my $script = $utl->{SCRIPTS}{PROJECTS_VIEW};
	
    foreach( @dbs ) {
		print $cgi->p("<a href='$script?db=$_'> $_ </a>");
    }
    print qq[ </fieldset> </div> ];
}
