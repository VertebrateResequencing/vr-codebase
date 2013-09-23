#!/usr/local/bin/perl -T
# Displays Qc Grind Samples View

use strict;
use warnings;
use URI;

use lib '../../lib';
use SangerPaths qw(core team145);
use SangerWeb;
use VRTrack::Project;
use VertRes::Utils::VRTrackFactory;
use VertRes::QCGrind::Util;
use CGI::Carp qw(fatalsToBrowser);

my $pending_view = "../pending-view/pending_view.pl";

my $title = 'Samples View';
my $sw  = SangerWeb->new({
    'title'   => $title,
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
});
my $utl = VertRes::QCGrind::Util->new();
my $main_script = $utl->{SCRIPTS}{DATABASES_VIEW};
my $proj_view_script = $utl->{SCRIPTS}{PROJECTS_VIEW};
my $study_lanes_script = $utl->{SCRIPTS}{STUDY_LANES};

my $USER = $sw->username();
my $cgi = $sw->cgi();
my $SCRIPT_NAME = $cgi->url(-relative=>1);

my $db = $cgi->param('db');
unless ($db) {
    print $sw->header();
	$utl->displayError( "No database ID",$sw );
}
my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $db, mode => 'rw');

my $projectID = $cgi->param('proj_id');
unless ($projectID) {
    print $sw->header();
	$utl->displayError( "No Project ID",$sw );
}
my $project = VRTrack::Project->new( $vrtrack, $projectID );

print $sw->header();
displaySamplesPage($cgi, $vrtrack, $db, $project);
print $sw->footer();
exit;

sub displaySamplesPage
{
    my ($cgi, $vrtrack, $database, $project) = @_;
    
    my $samples = $project->samples();
    $utl->displayError( "Cant get samples for project: $projectID",$sw ) unless $samples;
    
    my %individuals2Samples;
    foreach( @{ $samples } ) {
        my $sample = $_;
        my $ind = $sample->individual();
        my $indname = $ind->name();
        if( $individuals2Samples{ $indname } ) {
            push( @{ $individuals2Samples{ $indname } }, $sample );
        }
        else{$individuals2Samples{ $indname } = [ $sample ];}
    }
    
    my $pname = $project->name;

	print $cgi->h2({-align=>"center", -style=>"font: normal 900 1.5em arial"},"<a href='$main_script'>QC Grind</a> $title");
	print $cgi->h3({-align=>"center", -style=>"font: normal 700 1.5em arial"},"Database : <a href='$proj_view_script?db=$database'>$database</a>");
	print $cgi->h3({-align=>"center", -style=>"font: normal 700 1.5em arial"},"Project : $pname");
	print $cgi->h5({-style=>"font: arial"},"<a href='$pending_view?mode=2&amp;db=$database'>Pending Requests</a>");
	print qq[<div align=right><h5><a href="$study_lanes_script?;proj_id=$projectID&amp;db=$database">Lane View</a></h5></div>];

    print $cgi->br;

    print qq[
        <div class="centerFieldset">
        <fieldset style="width: 700px">
        <legend>Samples</legend>
        <table RULES=GROUPS width="100%">
        <tr>
        <th>Individual</th>
        <th>Sanger</th>
        <th>Library</th>
        <th>Lanes</th>
        <th>Passed</th>
        <th title="Passed bases mapped with dups removed">Pass seq</th>
        <th>Depth</th>
        <th>Pending</th>
        </tr>
    ];
    
    my $grandTotalLanes = 0;
    my $grandTotalPassed = 0;
    my $grandTotalDepth = 0;
    my $grandTotalNoQC = 0;
    my $grandTotalSamples = 0;

    foreach( sort( keys( %individuals2Samples ) ) ) {
        my $iname = $_;
        my @samples = @{$individuals2Samples{ $_ }};
        
        print qq[ <tr> ];
        
        print qq[
            <td>$iname</td>
            <td></td>
        ];

        my $sampleLanes = 0;
        my $sampleDepth = 0;
        my $samplePassed = 0;
        my $samplePassseq = 0;
        my $sample_no_qcLanes = 0;
        my $firstL = 1;

        foreach my $sample ( @samples ) {

			my @libraries = sort {$a->name cmp $b->name} @{$sample->libraries()};

			foreach my $library ( @libraries ) {
				my $lname = $library->name;
				my $lid = $library->id();
				my $lanes = $library->lanes();
				my $passedLanes = 0;
				my $passedBases = 0;
				my $sampledBases = 0;
				my $total_passed_bases = 0;
				my $lib_no_qcLanes = 0;
                my $final_net = 0;

				foreach my $lane( @$lanes ) {
					$sampleLanes ++;
					if( $lane->qc_status() eq $utl->{STATES}{PASSED} ) {
						$passedLanes++;
						my $mapping = $lane->latest_mapping();
						if( $mapping ) {
							$passedBases += $mapping->rmdup_bases_mapped;
							$sampledBases += $mapping->raw_bases();
							$total_passed_bases += $lane->raw_bases();

							$samplePassed ++;

                            # subtract Overlap duplicate bases from bases_mapped to get net total bases
                            my $overlap_dup = 0;
                            my @autoqc_statuses = @{ $mapping->autoqcs() };
                            foreach my $autoqc (@autoqc_statuses) {

                                next unless $autoqc->test =~ /Overlap duplicate base percent/;

                                # eg "The percent of bases duplicated due to reads of a pair overlapping (2.8) is smaller than or equal to 4."
                                $autoqc->reason =~ /The percent of bases duplicated due to reads of a pair overlapping \((.*)\) /;
                                $overlap_dup = sprintf("%.2f",$1);
                            }
                            $final_net += int($mapping->rmdup_bases_mapped - $mapping->rmdup_bases_mapped * $overlap_dup / 100);
						}
					}
					elsif( $lane->qc_status() eq $utl->{STATES}{PENDING} ) {
						$lib_no_qcLanes ++;
					}
				}

				my $depth = 0;
				my $pass_seq = 0;
				if( $total_passed_bases > 0 ) {
					$depth = ( ( $passedBases / $sampledBases ) * $total_passed_bases ) / 3000000000;
					$depth = sprintf("%.2f", $depth);
					$sampleDepth += $depth;
				}
                $pass_seq = $utl->bp_to_nearest_unit($final_net, 1);
                $samplePassseq += $final_net;

				my $colour = $utl->get_colour_for_status( $library->open() ? $library->qc_status() : $utl->{STATES}{CLOSE_LIBRARY} );
				my $numLanes = @$lanes;

				if( ! $firstL ){print qq[<tr><td></td><td></td>];$firstL=0;}

				print qq[ <td style="background-color:$colour;"><a href="$study_lanes_script?;proj_id=$projectID&amp;db=$database&amp;lib=$lname">$lname</a></td> ];

				print $numLanes > 0 ? qq[<td>$numLanes</td><td>$passedLanes</td>] : qq[<td></td><td></td>];
				print $total_passed_bases > 0 ? qq[<td>$pass_seq</td>] : qq[<td></td>];
				print $depth > 0 ? qq[<td>$depth].qq[x</td>] : qq[<td></td>];

				if( $lib_no_qcLanes > 0 ) {
					print qq[
						<td>$lib_no_qcLanes</td>
						</tr>
						];
				}
				else{print qq[<td></td></tr>];}
				$firstL = 0;
				$sample_no_qcLanes += $lib_no_qcLanes;
			}
		}
        
        $samplePassseq = $utl->bp_to_nearest_unit($samplePassseq, 1);
        print qq[<tr><th></th><th></th><th></th><th>$sampleLanes</th><th>$samplePassed</th><th>$samplePassseq</th><th>$sampleDepth].qq[x</th>];
        
        if( $sample_no_qcLanes > 0 ) {
            print qq[<th>$sample_no_qcLanes</a></th>];
        }
        else {
            print qq[<th></th>];
        }
        
        print qq[</tr>];
        $grandTotalLanes += $sampleLanes;
        $grandTotalPassed += $samplePassed;
        $grandTotalDepth += $sampleDepth;
        $grandTotalNoQC += $sample_no_qcLanes;
        $grandTotalSamples ++;
    }
    
    print qq[<tr><tfoot><th>$grandTotalSamples</th><th></th><th></th><th>$grandTotalLanes</th><th>$grandTotalPassed</th><th>$grandTotalDepth].qq[x];
    
    if( $grandTotalNoQC > 0 ) {
        print qq[<th>$grandTotalNoQC</th>];
    }
    else {
        print qq[<th></th>];
    }
    
    print qq[
        </tfoot></tr>
        </table>
        </fieldset>
        </div>
    ];
}
