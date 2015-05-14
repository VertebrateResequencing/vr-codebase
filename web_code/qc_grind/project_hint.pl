#!/usr/local/bin/perl -T
# get project name hint, called from javascript qc.js
use strict;
use warnings;
use CGI;
use CGI::Carp qw(fatalsToBrowser);

use lib '/var/www/lib';
use SangerPaths qw(core team145);
use VRTrack::Project;
use VertRes::Utils::VRTrackFactory;
use VertRes::QCGrind::Util;

my $utl = VertRes::QCGrind::Util->new();
my $cgi = new CGI;
my $q = $cgi->param('q');
my @dbs = VertRes::Utils::VRTrackFactory->databases(1);

print $cgi->header;
print '<fieldset>';
print "<table border=0>";

my $script = $utl->{SCRIPTS}{PROJECTS_VIEW};

foreach my $db ( sort @dbs ) {

    my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $db, mode => 'r');
    next unless $vrtrack;

    my @projects = sort {$a->name cmp $b->name} @{$vrtrack->projects()};
    foreach my $project ( @projects ) {
        my $projname = $project->name;
        if ( $projname =~ /$q/i) {
            print $cgi->Tr("<td>$projname</td><td><a href='$script?db=$db'> $db </a></td>");
        }
    }
}
print "</table>";
print "</fieldset>";

exit;
