#!/usr/local/bin/perl -T

use strict;
use warnings;
use URI;

use lib '/var/www/lib';
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
    'jsfile'  => ['http://jsdev.sanger.ac.uk/prototype.js','http://jsdev.sanger.ac.uk/toggle.js','ttp://jsdev.sanger.ac.uk/scriptaculous/scriptaculous.js','http://jsdev.sanger.ac.uk/sidebar.js','http://jsdev.sanger.ac.uk/urchin.js','http://jsdev.sanger.ac.uk/zebra.js','http://js.sanger.ac.uk/sorttable_v2.js',],
    'style'   => $utl->{CSS},
});

my $cgi = $sw->cgi();

my $db = $cgi->param('db');
my $pid = $cgi->param('proj_id');

unless ($db) {
    print $sw->header();
	$utl->displayError( "No database ID",$sw );
}
my $vrtrack = $utl->connectToDatabase($db);

$utl->displayError( "No Project ID",$sw ) unless $pid;

print $sw->header();
displayProjectLanesPage($cgi,$vrtrack,$db,$pid);
print $sw->footer();
exit;

#######

sub displayProjectLanesPage 
{
    my ($cgi, $vrtrack, $database, $projectID) = @_;
    
    my $init_script = $utl->{SCRIPTS}{MAP_VIEW};
    my $proj_script = $utl->{SCRIPTS}{MAP_PROJECTS_VIEW};
    my $index = $utl->{SCRIPTS}{INDEX_PAGE};
    my $qcgrind_script = $utl->{SCRIPTS}{QCGRIND_LANE};
    
    my $project = VRTrack::Project->new( $vrtrack, $projectID );
    displayError( "Cant get project: $projectID" ) unless $project;
    
    my $samples = $project->samples();
    displayError( "Cant get samples for project: $projectID" ) unless $samples;
    
    my $pname = $project->name;

    print qq[
        <h4 align="center" style="font: arial"><i><a href="$index">Team 145</a> : <a href="$init_script">$title</a> :  <a href="$proj_script?db=$database">$database</a></i> : $pname</h4>
    ];
    
    print qq[
    
    <div class="centerFieldset">
    <fieldset > 
    <legend>Lane data</legend>
    <table width="60%">
    <tr>
    <th>Library</th>
    <th>Improved</th>
    <th>Called</th>
    <th>Name</th>
    ];
    
    my @lanes;
    my %mappers;
    foreach( sort { $a->ssid() <=> $b->ssid() } @$samples)
    {
        my @libraries = sort {$a->name cmp $b->name} @{$_->libraries()};
        foreach( @libraries ) 
        {
            my $library = $_;
            my  @lanes_ = @{ $library->lanes() };
            foreach my $lane ( @lanes_ ) 
            {
                push( @lanes, $lane );
                my @mapstats = @{ $lane->mappings() };
                foreach my $mapstat (@mapstats)
                {
                    if( $mapstat->mapper() )
                    {
                        my $mapper = $mapstat->mapper();
                        $mappers{ $mapper->name().qq[ v].$mapper->version() } = 1;
                    }
                }
            }
        }
    }
    
    foreach my $mname ( sort( keys( %mappers ) ) )
    {
        print qq[<th>$mname</th>];
    }
    print qq[</tr>];
    
    foreach my $lane ( @lanes )
    {
        print qq[<tr>];
        my $id = $lane->id();
        my $library = VRTrack::Library->new($vrtrack, $lane->library_id());
        my $libName = $library->name();
        my $improved = $lane->is_processed('improved') ? 'yes' : 'no';
        my $called = $lane->is_processed('snp_called') ? 'yes' : 'no';
        print qq[<td>$libName</td><td>$improved</td><td>$called</td><td><a href="$qcgrind_script?lane_id=$id&db=$database">].$lane->name().qq[</a></td>];
        my @mappings = @{ $lane->mappings() };
        my %lane_mappers;
        foreach my $mapstat ( @mappings )
        {
            if( $mapstat->mapper() && $mapstat->raw_bases() && $mapstat->raw_bases() > 0 )
            {
                my $mapper = $mapstat->mapper();
                $lane_mappers{ $mapper->name().qq[ v].$mapper->version() } = 1;
            }
        }
        foreach my $mapper ( sort( keys( %mappers ) ) )
        {
            if( $lane_mappers{ $mapper } )
            {
                print qq[<td>yes</td>];
            }
            else{print qq[<td>no</td>];}
        }
        print qq[</tr>];
    }
    print qq[
        </table>
        </fieldset>
        </div>
    ];
}
