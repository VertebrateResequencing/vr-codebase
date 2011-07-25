#!/usr/bin/env perl
#
# lane_bams.pl --root /abs/META --db g1k_meta --all_samples \
#              --slx_mapper bwa --454_mapper ssaha --assembly_name NCBI37 \
#              --ignore_platform SOLID --improved
#              --project_regex => "low_coverage" > lane_bams.fofn
#
#
# Author: Sendu Bala <bix@sendu.me.uk>
#

use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use File::Spec;
use VertRes::Utils::Hierarchy;
use VertRes::Utils::VRTrackFactory;
use VRTrack::Lane;
use VRTrack::History;

# get user input
my ($help, $database, $all_samples, @ignore_platforms, @samples, $root, $qc,
    $improved, $date, $mapper_slx, $mapper_454, @mapper_alias_slx, @mapper_alias_454, 
    $assembly, $project_regex, $coverage_threshold, $fai);
my $spinner = 0;
GetOptions('db=s'              => \$database,
           'all_samples'       => \$all_samples,
           'samples:s{1,}'     => \@samples,
           'root=s'            => \$root,
           'slx_mapper=s'      => \$mapper_slx,
           '454_mapper=s'      => \$mapper_454,
           'slx_mapper_alias:s{1,}'      => \@mapper_alias_slx,
           '454_mapper_alias:s{1,}'      => \@mapper_alias_454,
           'assembly_name=s'   => \$assembly,
           'ignore_platform=s' => \@ignore_platforms,
           'project_regex=s'   => \$project_regex,
           'improved'          => \$improved,
           'qc'                => \$qc,
           'coverage=f'        => \$coverage_threshold,
           'fai=s'             => \$fai,
           'date=s'            => \$date,
           'h|help'            => \$help);

my $missing_opts = 0;
unless ($database && ($all_samples || @samples) && $root && $mapper_slx && $mapper_454 && $assembly) {
    $missing_opts = 1;
}

($help || $missing_opts) and die <<USAGE;
Prints out the paths to all the lane bams that would be needed to build
sample or platform-level merged bams. This file can be used as input to the
MergeUp pipeline. It will only output lanes for a sample if all that sample's
lanes have been mapped (ignoring lanes matching an --ignore_platform option).

Usage: $0 [required args] [optional args]
Required:
        --db               <database name>
        --all_samples      get all lane paths in the database
        (or --samples      <sample name> [<sample name> ...])
        --slx_mapper       <mapper name> eg. bwa
        --454_mapper       <mapper name> eg. ssaha
        --assembly_name    <assembly name> eg. NCBI37

Optional:
        --slx_mapper_alias alternate mapper names eg. bwa_aln (can supply more 
                                           than one - separate by space)
        --454_mapper_alias as with slx_mapper_alias
        --improved         choose bams that have been run through the
                           BamImprovement pipeline
        --qc               consider only lanes that are set as qc passed
        --coverage         <float> minimum coverage per sample to consider 
                             (requires fai option to also be set)
        --fai              <fai_file> path to reference fai index file 
                             (required of --coverage option set)
        --ignore_platform  <SLX|454|SOLID> ignore lanes sequenced with this
                                           platform (this whole option can be
                                           specified more than once)
        --project_regex    "regex" ignore lanes that belong to projects that do
                           not match the given regular expression
        --date             '2010-01-04 10:49:10' view the database as it was
                                                 immediately prior to this date
        --help             this message

USAGE

if ($coverage_threshold && !$fai) {
    die "Must supply fai option if coverage option is set\n";
}

if ($fai && !(-s $fai)) {
    die "Could not find reference index fai, $fai\n";
}

if ($all_samples && @samples) {
    warn "Both --all_samples and --samples were set; ignoring the --all_samples request\n";
    $all_samples = 0;
}

# travel back in time?
my $hist = VRTrack::History->new();
$hist->time_travel($date) if $date;

# get all the lanes for our desired samples, excluding ignored platforms
my @platforms;
my %ignore_platforms = map { $_ => 1 } @ignore_platforms;
foreach my $platform ('SLX', '454', 'SOLID') {
    next if exists $ignore_platforms{$platform};
    push(@platforms, $platform);
}

my %cd = VertRes::Utils::VRTrackFactory->connection_details('r');
$cd{database} = $database;

my $hu = VertRes::Utils::Hierarchy->new();
my @lanes = $hu->get_lanes(db => { %cd },
                           $all_samples ? () : (sample => [@samples]),
                           $project_regex ? (project_regex => $project_regex) : (),
                           platform => [@platforms]);

if ( !@lanes ) { die "No lanes found??\n"; }

# get the paths of all lanes, and determine which samples are fully mapped
# (now, not at $date)
$hist->time_travel('latest');
my $vrtrack = $lanes[0]->vrtrack;

my %bad_samples;
my %lanes_by_sample;
my %seen_paths;
foreach my $lane (@lanes) {
    if ($date) {
        $lane = VRTrack::Lane->new($vrtrack, $lane->id);
    }
    my %objs = $hu->lane_hierarchy_objects($lane);
    my $sample = $objs{sample}->name;
    
    next if exists $bad_samples{$sample};
    next if( $qc && ( !$lane->qc_status() || $lane->qc_status() ne 'passed' ) );
    
    # are we mapped/improved?
    unless ($lane->is_processed('mapped') && $improved ? $lane->is_processed('improved') : 1) {
        my $lane_name = $lane->name;
        my $problem = $improved ? 'mapped/improved' : 'mapped';
        warn "$lane_name was not $problem!\n";
        $bad_samples{$sample} = 1;
        delete $lanes_by_sample{$sample};
        warn "Not all the lanes under sample '$sample' were $problem; they will all be excluded\n";
        next;
    }
    
    my $path = File::Spec->catdir($root, $vrtrack->hierarchy_path_of_lane($lane));
    if (exists $seen_paths{$path}) {
        warn "already saw $path...\n";
        next;
    }
    $seen_paths{$path} = 1;
    push(@{$lanes_by_sample{$sample}}, $path);
}

# samples exceed coverage threshold?
if ($coverage_threshold) {
    my $genome_size = genome_size($fai);
    foreach my $sample (keys %lanes_by_sample) {
        my $coverage = $hu->hierarchy_coverage(sample => [$sample], genome_size => $genome_size, qc_passed => $qc, vrtrack => $vrtrack);
        next unless ($coverage < $coverage_threshold);
        $bad_samples{$sample} = 1;
        delete $lanes_by_sample{$sample};
        warn "'$sample' did not exceed coverage threshold ($coverage < $coverage_threshold) and will be excluded\n";
    }
}

# output good lane bams
while (my ($sample, $lanes) = each %lanes_by_sample) {
    foreach my $path (@{$lanes}) {
        # Allow to skip bad lanes
        my @bams;
        eval { @bams = $hu->lane_bams($path,
                                  vrtrack => $vrtrack,
                                  slx_mapper => $mapper_slx,
                                  '454_mapper' => $mapper_454,
                                  slx_mapper_alias => [@mapper_alias_slx],
                                  '454_mapper_alias' => [@mapper_alias_454],
                                  assembly_name => $assembly);
        };
        if ( $@ ) {
            warn "$path:\n\t$@\n";
            next;
        }
        
        foreach my $bam (@bams) {
            print $bam, "\n";
        }
    }
}

exit;


sub genome_size {
    my ($fai) = @_;
    
    my $genome_size = 0;
    open my $fh,'<',$fai || die("$fai: $!"); 
    while (<$fh>)
    {
        my (undef, $length, undef) = split /\t/;
        next unless $length;
        $genome_size += $length;
    }
    close $fh;
    return $genome_size;
}

