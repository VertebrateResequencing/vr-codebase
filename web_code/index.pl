#!/usr/local/bin/perl -T
use strict;
use warnings;
use URI;

use SangerPaths qw(core team145);
use SangerWeb;
use VertRes::QCGrind::ViewUtil;

my $utl = VertRes::QCGrind::ViewUtil->new();

my $title = 'Team 145 Web Tools';

my $sw  = SangerWeb->new({
    'title'   => $title,
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
    'style'   => $utl->{CSS},
});

my $cgi = $sw->cgi();

print $sw->header();

print qq[ <h2 align="center" style="font: normal 900 1.5em arial">Team 145 Web Tools</a></h2> ];

print qq[<div class="centerFieldset">
        <fieldset style="width: 550px">
            <legend>Select web tool</legend>
    
            <p><a href="qc_grind/qc_grind.pl">QC Grind</a><br/>Display and/or download QC data for project lanes</p>
        
            <p><a href="sample-mapping/sample_mapping.pl">Sample ID Mapping</a><br/>Download and/or view mappings of sample/supplier names and accession numbers</p>
        
            <p><a href="map-view/map_view.pl">Map View</a><br/>Display QC information for short-read sequencing to allow lanes to be passed/failed</p>
        
            <p><a href="pending-view/pending_view.pl">Pending View</a><br/>Display pending/started requests for libraries and sequences</p>
           
        </fieldset>
        </div>];

print $sw->footer();
exit;

