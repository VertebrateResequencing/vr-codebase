#!/usr/bin/env perl
#
# track_to_si.pl my_track_db > my.sequence.index
#
# Pulls out all the QC-passed VRTrack tracking database lanes and forms a
# sequence.index from the information, which would be suitable for use with
# update_vrmeta.pl for creating a corresponding VRTrack meta database.
#
# Author: Sendu Bala <bix@sendu.me.uk>
#

use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use VertRes::Utils::VRTrackFactory;
use VertRes::Utils::Hierarchy;

# get user input
my ($help);
my $spinner = 0;
GetOptions('h|help'      => \$help);

my $database = shift;

my $missing_opts = 0;
unless ($database) {
    $missing_opts = 1;
}

($help || $missing_opts) and die <<USAGE;
Creates a sequence.index from a VRTrack tracking database.

Usage: $0 tracking_database_name > sequence.index
        --help        <this message>

USAGE

# setup database connection
my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $database,
                                                          mode => 'r');
unless ($vrtrack) {
    die "DB connection failed: ".$DBI::errstr."\n";;
}

my $hu = VertRes::Utils::Hierarchy->new();

# update_vrtrack.pl relies on a samples file to supply information that isn't
# stored in sequence.index files, and the sample file is keyed on first column
# to get the other info. However, the mouse samples file has first column
# entries that are not in the mouse tracking database... we just have a hack
# hash lookup for now...
my %sample_lookup = (
    'CAST/EiJ' => 'CAST',
    'PWK/PhJ' => 'PWK',
    'SPRET/EiJ' => 'SPRET',
    'C57BL/6NJ' => 'C57BL',
    'C3H/HeJ' => 'C3H',
    'CBA/J' => 'CBA',
    'A/J' => 'A_J',
    'AKR/J' => 'AKR',
    '129P2/OlaHsd' => '129P2',
    'NOD/ShiLtJ' => 'NOD',
    '129S1/SvImJ' => '129S1',
    'LP/J' => 'LP_J',
    'WSB/EiJ' => 'WSB',
    '129S5/SvEvBrd' => '129S5',
    'NZO/HlLtJ' => 'NZO',
    'BALB/cJ' => 'BALB',
    'DBA/2J' => 'DBA');

# loop through all QC-passed lanes *** this method auto-excludes withdrawn
# lanes, which we actually wanted theoretically, but we'll leave it ok for now
my $passed_lanes = $vrtrack->qc_filtered_lane_hnames('passed');
foreach my $lane_name (@{$passed_lanes}) {
    my %lane_info = $hu->lane_info($lane_name, vrtrack => $vrtrack);
    next unless $lane_info{imported};
    my $lane = $lane_info{vrlane};
    
    my @si;
    #[0]  FASTQ_FILE
    #[1]  MD5
    #[2]  RUN_ID
    #[3]  STUDY_ID
    #[4]  STUDY_NAME
    #[5]  CENTER_NAME
    #[6]  SUBMISSION_ID
    #[7]  SUBMISSION_DATE
    #[8]  SAMPLE_ID
    #[9]  SAMPLE_NAME
    #[10] POPULATION
    #[11] EXPERIMENT_ID
    #[12] INSTRUMENT_PLATFORM
    #[13] INSTRUMENT_MODEL
    #[14] LIBRARY_NAME
    #[15] RUN_NAME
    #[16] RUN_BLOCK_NAME
    #[17] INSERT_SIZE
    #[18] LIBRARY_LAYOUT
    #[19] PAIRED_FASTQ
    #[20] WITHDRAWN
    #[21] WITHDRAWN_DATE
    #[22] COMMENT
    #[23] READ_COUNT
    #[24] BASE_COUNT
    $si[0] = undef;
    $si[1] = undef;
    $si[2] = $lane_name;
    $si[3] = $lane_info{project};
    $si[4] = '?';
    $si[5] = $lane_info{centre};
    $si[6] = '?';
    $si[7] = '?';
    $si[8] = $lane_info{individual_acc};
    $si[9] = $sample_lookup{$lane_info{individual}} || $lane_info{individual};
    $si[10] = $lane_info{population};
    $si[11] = '?';
    $si[12] = $lane_info{technology};
    $si[13] = '?';
    $si[14] = $lane_info{library_raw};
    $si[15] = '?';
    $si[16] = '?';
    $si[17] = $lane_info{insert_size};
    $si[18] = '?';
    $si[19] = undef;
    $si[20] = $lane_info{withdrawn};
    $si[21] = '';
    $si[22] = '';
    $si[23] = undef;
    $si[24] = undef;
    
    my @files = @{$lane->files};
    my $paired = @files > 1 ? 1 : 0;
    foreach my $file (@files) {
        if ($paired) {
            # sequencescape can give single files for paired sequenceing, where
            # the filename doesn't match the lane name. Petr's internal import
            # pipeline splits these files into forward/reverse fastqs and stores
            # their names in the database; we skip the original and only take
            # the splits
            my ($type) = $file->name =~ /${lane_name}_(\d)/;
            next unless $type;
        }
        
        $si[0] = $lane_info{hierarchy_path}.'/'.$file->name.'.gz';
        $si[1] = $file->md5;
        $si[19] = $paired;
        $si[23] = $file->raw_reads;
        $si[24] = $file->raw_bases;
        
        print join("\t", @si), "\n";
    }
}

exit;
