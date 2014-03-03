#!/usr/local/bin/perl -T

use strict;
use warnings;
use URI;

use lib '/var/www/lib';
use SangerPaths qw(core team145);
use SangerWeb;
use VRTrack::Project;
use VertRes::Utils::VRTrackFactory;
use VertRes::QCGrind::ViewUtil;
use CGI::Carp qw(fatalsToBrowser);

my $title = 'HipSci Cohorts Viewer';
my $utl = VertRes::QCGrind::ViewUtil->new();
my $sw  = SangerWeb->new({
    'title'   => $title,
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
    'style'   => $utl->{CSS},
});

my $main_script = $utl->{SCRIPTS}{DATABASES_VIEW};

my $proj_view_script = $utl->{SCRIPTS}{PROJECTS_VIEW};

my $study_lanes_script = $utl->{SCRIPTS}{STUDY_LANES};

my $genotyping_db = $utl->{VRTRACK_DATABASES}{GENOTYPING};
my $expression_db = $utl->{VRTRACK_DATABASES}{EXPRESSION};

my $cgi = $sw->cgi();
my $SCRIPT_NAME = $cgi->url(-relative=>1);

print $sw->header();
displaySamplesPage($cgi, $genotyping_db, $expression_db);
print $sw->footer();
exit;

sub displaySamplesPage
{
    my ($cgi, $genotyping_db, $expression_db) = @_;
    
    # will need to add biosample ids here when avaliable -> poss consider
    # converting these arrays into hashes to contain (cohort_id, biosample_id) pairs?
    my @genotyping_cohorts = $utl->getHipsciCohorts ($genotyping_db);
    my @expression_cohorts = $utl->getHipsciCohorts ($expression_db);
    
    my $summary_view_script = $utl->{SCRIPTS}{HIPSCI_SUMMARY_VIEW};
	my $detailed_view_script = $utl->{SCRIPTS}{HIPSCI_DETAILED_VIEW};
	my $index = $utl->{SCRIPTS}{INDEX_PAGE};
	
	my %hipsci_cohorts;
	for my $cohort (@expression_cohorts, @genotyping_cohorts){
		$hipsci_cohorts{$cohort} = 1;
	}
	
	print qq[ <h4 align="center" style="font: arial"><i><a href="$index">Team 145</a></i> : $title</h4><br/> ];

    print qq[
        <div class="centerFieldset">
        <fieldset style="width: 700px">
        <legend>Donor information</legend>
        <table RULES=GROUPS width="100%">
        <tr>
        <th>Cohort identifier</th>
        <th>Control sample</th>
        <th></th>
        </tr>
    ];
    
    # my $biosample_id = "N/A"; *** what is this supposed to be? There used to
    # be a pointless BioSample ID column in the output that was hardcoded to be
    # N/A??
    
    for my $cohort( sort( keys( %hipsci_cohorts ) ) ) {
        my $control_sample = $utl->getControlSample($genotyping_db, $cohort);
        
        print qq[
          	<tr>
            	<td>$cohort</td>
                <td>$control_sample</td>
                <td><a href="$detailed_view_script?cohort=$cohort">Detailed view</a></td>                
            </tr>
        ];
    }
    print qq[
        </table>
        </fieldset>
        </div>
    ];
}
