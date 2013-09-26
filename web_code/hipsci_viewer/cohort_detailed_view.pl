#!/usr/local/bin/perl -T

use strict;
use warnings;
use lib '/var/www/lib';
use URI;

use SangerPaths qw(core team145);
use SangerWeb;

use VRTrack::Project;
use VertRes::Utils::VRTrackFactory;
use VertRes::QCGrind::ViewUtil;
use CGI::Carp qw(fatalsToBrowser);

my $utl = VertRes::QCGrind::ViewUtil->new();
my $title = ' : Detailed View';
my $sw  = SangerWeb->new({
    'title'   => $title,
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
    'style'   => $utl->{CSS},
});

my $main_script = $utl->{SCRIPTS}{HIPSCI_COHORTS_VIEW};
my $genotyping_db = $utl->{VRTRACK_DATABASES}{GENOTYPING};
my $expression_db = $utl->{VRTRACK_DATABASES}{EXPRESSION};

my $cgi = $sw->cgi();
print $sw->header();

my $cohort = $cgi->param('cohort');
$utl->displayError( "No Cohort identifier",$sw ) unless $cohort;

displayDetailedInformation( $cgi, $genotyping_db, $expression_db, $cohort );
print $sw->footer();
exit;

#####

sub displayDetailedInformation
{
    my ( $cgi, $genotyping_db, $expression_db, $cohort) = @_;
    
    my $penncnv_flag = 'penncnv';
    my $quantisnp_flag = 'quantisnp';
    
    my $lane_id = $utl->getExpressionLaneID($expression_db, $cohort);
    my ( $vrtrack, $lane, $mapstats );
    if ( $lane_id ) {
		$vrtrack = $utl->connectToDatabase($expression_db);
		$lane = VRTrack::Lane->new( $vrtrack, $lane_id );
		displayError("Unable to retrieve lane $lane_id\n") unless $lane;
		$mapstats = $lane->latest_mapping();
	}
    
    my %geno_sample_controls = $utl->getGenotypingSamples($genotyping_db, $cohort);
    my %geno_cnv_totals = $utl->getGenotypingCNVTotals($genotyping_db, $cohort);
    my $penncnv_least = $utl->getGenotypingCNVLeastValue($genotyping_db, $cohort, $penncnv_flag);
    my $quancnv_least = $utl->getGenotypingCNVLeastValue($genotyping_db, $cohort, $quantisnp_flag);
    my @geno_samples = keys %geno_sample_controls;
    
    # Is there any chance we can have > 1 control?? Use array JUST IN CASE.....
    my @geno_control;
    foreach ( @geno_samples ) {
		push @geno_control, $_ unless $geno_sample_controls{$_} == 2;
	}
	my $control_string = join(', ', @geno_control);
    
    print $cgi->h2({-align=>"center", -style=>"font: normal 900 1.5em arial"},"<a href='$main_script'>HipSci Cohort Viewer</a> $title");
	print qq[ <h3 align="center" style="font: normal 700 1.5em arial"> Cohort identifier : $cohort </h3>];
	print qq[ <h3 align="center" style="font: normal 700 1.5em arial"> Control sample:  $control_string</h3>];
	
    if ($mapstats) {
		print qq[
            <div class="centerFieldset">
            <fieldset style="width: 800px">
            <legend>Cell line Pluritest plots</legend>
            <table width="80%">
            <tr>
        ];
		
		my $imglist = $mapstats->images;
        my %images = map {$_->name => $_} @$imglist;
        my @big_images = qw (pluritest_image03.png pluritest_image03c.png);
        foreach (@big_images) {
			my $img = $images{$_};
			if ($img) {
                print qq[<td>];$utl->printFullImage($img);print qq[</td>];
                delete $images{$_};
            }
        }
        
        print qq[
            </tr>
            </table>
            <table>
            <tr>
        ];
        
        my @i = sort keys %images;
        foreach( @i ) {
            print qq[<td>];$utl->printPreviewImage($images{ $_ });print qq[</td>];
        }
        
        print qq[
            </tr>
            </table>
            </fieldset>
            </div>
        ];
    }
    else {
		print qq[
            <div class="centerFieldset">
            <fieldset style="width: 800px">
            <legend>Cell line Pluritest plots</legend>
            <table width="80%">
            <tr><td><i>No images found at present.</i></td></tr>
            </table>
            </fieldset>
            </div>
        ];
	}
	
	my @cnv_samples = keys %geno_cnv_totals;
	if ( @cnv_samples > 0 ) {
		print qq[
            <div class="centerFieldset">
            <fieldset style="width: 800px">
        	<legend>Cell line CNV data</legend>
        	<table RULES=GROUPS width="100%">
        	<tr>
        	<th>Cell line identifier</th>
        	<th>PennCNV total</th>
        	<th>PennCNV diff</th>
	        <th>QuantiSNP total</th>
        	<th>QuantiSNP diff</th>
        	</tr>
    	];

    	foreach ( @cnv_samples ) {
			my @penn_values = @{$geno_cnv_totals{$_}{$penncnv_flag}};
			my @quanti_values = @{$geno_cnv_totals{$_}{$quantisnp_flag}};
			print qq[
          		<tr>
            		<td>$_</td>
            		<td>$penn_values[0]</td>];
        	if ( $penn_values[1] == $penncnv_least ) {
				print qq[	
			    	<td bgcolor="LawnGreen">$penn_values[1]</td>
			    	<td>$quanti_values[0]</td>];
			}
			else {
				print qq[
            		<td>$penn_values[1]</td>
            		<td>$quanti_values[0]</td>];
			}
        	if ( $quanti_values[1] == $quancnv_least ) {
				print qq[	
			    	<td bgcolor="LawnGreen">$quanti_values[1]</td>];
			}
			else {
				print qq[
            		<td>$quanti_values[1]</td>                
            		</tr>];
			}
		}
    
    	print qq[
        	</table>
        	</fieldset>
        	</div>
    	];
    }
    else {
		print qq[
            <div class="centerFieldset">
            <fieldset style="width: 800px">
            <legend>Cell line CNV data</legend>
            <table width="80%">
            <tr><td><i>No CNV data found at present.</i></td></tr>
            </table>
            </fieldset>
            </div>
        ];
	}
        
}

