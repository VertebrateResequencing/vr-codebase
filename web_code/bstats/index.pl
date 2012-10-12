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
            <legend>Hello World</legend>
    
            <p>Nothing fancy in here yet.</p>
           
        </fieldset>
        </div>];

print $sw->footer();
exit;

