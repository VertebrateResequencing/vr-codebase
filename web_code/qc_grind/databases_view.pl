#!/usr/local/bin/perl -T
# Displays vtrack databases

use strict;
use warnings;
use URI;
use CGI::Carp qw(fatalsToBrowser);

use SangerPaths qw(core team145);
use SangerWeb;
use VRTrack::Project;
use VertRes::Utils::VRTrackFactory;
use VertRes::QCGrind::Util;

my $title = 'QC Grind Databases';
my $sw  = SangerWeb->new({
    'title'   => $title,
    'jsfile'  => ['http://code.jquery.com/jquery-latest.js','/Teams/Team145/js/qc.js','/Teams/Team145/js/jquery.tablesorter.min.js'],
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

	my $script = $utl->{SCRIPTS}{PROJECTS_VIEW};
	print qq[ <h2 align="center" style="font: normal 900 1.5em arial">$title</h2> ];

    # project search
    print $cgi->p({-align=>"right"},'Find Project <input type="text" id="projFind" />');
    print $cgi->p('<span id="projHint"></span>');

	my @main_dbs = qw (vrtrack_human_wgs vrtrack_human_wes_v5 vrtrack_mouse_wgs_GRCm38 vrtrack_mouse_wes_GRCm38 vrtrack_scerevisiae_wgs);
    print qq[
        <fieldset style="width: 500px">
        <legend>Main Databases</legend>
    ];
    foreach( @main_dbs ) {
		print $cgi->p("<a href='$script?db=$_'> $_ </a>");
    }
    print qq[ </fieldset> ];

	my @uk10k_dbs = qw (vrtrack_uk10k_cohort vrtrack_uk10k_neuro vrtrack_uk10k_obesity vrtrack_uk10k_rare);
    print qq[
        <fieldset style="width: 500px">
        <legend>UK10K Databases</legend>
    ];
    foreach( @uk10k_dbs ) {
		print $cgi->p("<a href='$script?db=$_'> $_ </a>");
    }
    print qq[ </fieldset> ];

	my %done;
	push (@main_dbs, @uk10k_dbs);
	foreach (@main_dbs) {
		$done{$_}++;
	}
    my @dbs = VertRes::Utils::VRTrackFactory->databases(1);
    print qq[
        <fieldset style="width: 500px">
        <legend>Other Databases</legend>
    ];
    foreach( @dbs ) {
		next if $done{$_};
		print $cgi->p("<a href='$script?db=$_'> $_ </a>");
    }
    print qq[ </fieldset> ];

}
