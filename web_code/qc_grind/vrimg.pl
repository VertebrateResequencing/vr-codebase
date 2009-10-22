#!/usr/local/bin/perl -T
# vrimg.pl
# 
# Returns image from tracking databases
#
# Author:        jws
# Maintainer:    jws
# Created:       2009-08-10

use strict;
use warnings;
no warnings 'uninitialized';

use SangerPaths qw(core team145);
use VRTrack::VRTrack;
use VRTrack::Image;
use SangerWeb;

my $sw  = SangerWeb->new({
    'title'   => q(vrimg),
});

my $cgi = $sw->cgi();

$|++;

my $spp		= $cgi->param('spp');
my $img_id	= $cgi->param('img');
my %db_for_spp = (  'g1k'   => 'g1k_track',
                    'mouse' => 'mouse_reseq_track',
                    #'mouse' => 'mouse_cancer_track',
                );

my $db = $db_for_spp{$spp} or die "$spp not recognised";

my $vrtrack = VRTrack::VRTrack->new({ host => 'mcs4a',
                                    port => 3306,
                                    user => 'vreseq_ro',
                                    database => $db,
                                   });

unless ($vrtrack){
    die "Couldn't connect to vrtrack database";
}

my $image = VRTrack::Image->new($vrtrack,$img_id);
unless ($image){
    die "Can't retrieve image $img_id";
}
$image->name =~ /\.(\w+)$/;
my $type = $1;
print $cgi->header("image/$type");
print $image->image;
exit 1;
