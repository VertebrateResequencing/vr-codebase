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
use VRTrack::Multiplex_pool;

# get user input
my ($help, $database, $all_samples, @ignore_platforms, @samples, $root, $qc, $auto_qc, $gt, 
    $improved, $date, $mapper_slx, $mapper_454, @mapper_alias_slx, @mapper_alias_454, 
    $assembly, $project_regex, $seq_threshold, $strict, $no_requests, $verbose);
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
           'auto_qc'           => \$auto_qc,
           'gt'                => \$gt,
           'strict'            => \$strict,
           'no_requests'       => \$no_requests,
           'sequence=f'        => \$seq_threshold,
           'date=s'            => \$date,
           'verbose'           => \$verbose,
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
        --auto_qc          consider only lanes that have been auto_qc passed
                            (however, will not exclude if qc_status is passed)
        --gt               consider only lanes that have passed genotype check
        --sequence         <float> minimum amount of raw sequence per sample
        --strict           if any lane for a sample fails one of the improved/qc/auto_qc/gt 
                             checks, then the whole sample is excluded
        --ignore_platform  <SLX|454|SOLID> ignore lanes sequenced with this
                                           platform (this whole option can be
                                           specified more than once)
        --project_regex    "regex" ignore lanes that belong to projects that do
                           not match the given regular expression
        --date             '2010-01-04 10:49:10' view the database as it was
                                                 immediately prior to this date
        --verbose          be extra verbose
        --help             this message

USAGE

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
    my $project = $objs{project}->name;
    my $sample = $objs{sample}->name;
    
    next if exists $bad_samples{"$project/$sample"};
    
    if( $qc && ( !$lane->qc_status() || $lane->qc_status() ne 'passed' ) ) {
        my $lane_name = $lane->name;
        warn "$lane_name has not passed qc!\n" if $verbose;
        if ($strict) {
            $bad_samples{"$project/$sample"} = 1;
            delete $lanes_by_sample{"$project/$sample"};
            warn "Not all the lanes under sample '$project/$sample' have passed qc; they will all be excluded\n";
        }
        next;
    }
    
    if( $auto_qc && ( !$lane->auto_qc_status() || $lane->auto_qc_status() ne 'passed' ) ) {
        my $lane_name = $lane->name;
        if ( $lane->qc_status() eq 'passed' ) {
            warn "$lane_name has not passed auto_qc, but is marked as qc passed... keeping\n" if $verbose;
        } else {
            warn "$lane_name has not passed auto_qc!\n" if $verbose;
            if ($strict) {
                $bad_samples{"$project/$sample"} = 1;
                delete $lanes_by_sample{"$project/$sample"};
                warn "Not all the lanes under sample '$project/$sample' have passed auto_qc; they will all be excluded\n";
            }
            next;
        }
    }
    
    if( $gt && ( !$lane->genotype_status() || $lane->genotype_status() ne 'confirmed' ) ) {
        my $lane_name = $lane->name;
        if ( $lane->qc_status() eq 'passed' ) {
            warn "$lane_name genotype not confirmed, but is marked as qc passed... keeping\n" if $verbose;
        } else {
            warn "$lane_name genotype not confirmed!\n" if $verbose;
            if ($strict) {
                $bad_samples{"$project/$sample"} = 1;
                delete $lanes_by_sample{"$project/$sample"};
                warn "Not all the lanes under sample '$project/$sample' have confirmed genotype; they will all be excluded\n";
            }
            next;
        }
    }
    
    # are we mapped/improved?
    unless ($lane->is_processed('mapped') && ($improved ? $lane->is_processed('improved') : 1)) {
        my $lane_name = $lane->name;
        my $problem = $improved ? 'mapped/improved' : 'mapped';
        warn "$lane_name was not $problem!\n" if $verbose;
        if ($strict) {
            $bad_samples{"$project/$sample"} = 1;
            delete $lanes_by_sample{"$project/$sample"};
            warn "Not all the lanes under sample '$project/$sample' were $problem; they will all be excluded\n";
        }
        next;
    }
    
    my $path = File::Spec->catdir($root, $vrtrack->hierarchy_path_of_lane($lane));
    if (exists $seen_paths{$path}) {
        warn "already saw $path...\n";
        next;
    }
    $seen_paths{$path} = 1;
    push(@{$lanes_by_sample{"$project/$sample"}}, $path);
}

if ($seq_threshold) {
    while (my ($project_sample, $lanes) = each %lanes_by_sample) {
        my ($project, $sample) = split(/\//, $project_sample);
        my $sequence = 0;
        foreach my $lane_path (@{$lanes}) {
            my $lane_name = basename($lane_path);
            my $lane = VRTrack::Lane->new_by_hierarchy_name($vrtrack, $lane_name);
            $sequence += $lane->raw_bases()/1e9;
        }
        next unless ($sequence < $seq_threshold);
        delete $lanes_by_sample{"$project/$sample"};
        warn "'$project/$sample' did not exceed sequencing threshold ($sequence Gb < ${seq_threshold} Gb) and will be excluded\n";
    }
}

# Skip samples with pending library requests or incomplete sequencing requests
if ($no_requests) {
    # status: ('unknown', 'pending', 'started', 'passed', 'failed', 'cancelled', 'hold')
    my @incomplete = ('pending', 'started', 'hold');
    SAMPLE: foreach my $project_sample (keys %lanes_by_sample) {
        my ($project, $sample) = split /\//, $project_sample;
        my $vr_project = $vrtrack->get_project_by_name($project);
        my $vr_sample = $vr_project->get_sample_by_name($sample);
        
        my $library_requests = $vr_sample->library_requests();
        foreach my $librequest ( @{$library_requests} ) {
            my $status = $librequest->prep_status();
            if ($status eq 'pending') { 
                delete $lanes_by_sample{"$project/$sample"};
                warn "There are pending library requests for sample '$project/$sample', it will be excluded\n";
                next SAMPLE;
            }
        }
        
        my $libraries = $vr_sample->libraries();
        foreach my $library ( @{$libraries} ) {
            my $seq_requests = $library->seq_requests();
            foreach my $seqrequest ( @{$seq_requests} ) {
                my $status = $seqrequest->seq_status();
                if (grep /$status/, @incomplete) { 
                    delete $lanes_by_sample{"$project/$sample"};
                    warn "Incomplete library sequencing requests ($status) for sample '$project/$sample', it will be excluded\n";
                    next SAMPLE;
                }
            }
            
            foreach ( @{$library->library_multiplex_pools}){
                my $mplex = VRTrack::Multiplex_pool->new($vrtrack, $_->multiplex_pool_id);
                foreach my $seqrequest ( @{ $mplex->seq_requests } ) {
                    my $status = $seqrequest->seq_status();
                    if (grep /$status/, @incomplete) { 
                        delete $lanes_by_sample{"$project/$sample"};
                        warn "Incomplete multiplex sequencing requests ($status) for sample '$project/$sample', it will be excluded\n";
                        next SAMPLE;
                    }
                }
            }
        }
    }
}

# output good lane bams
while (my ($sample, $lanes) = each %lanes_by_sample) {
    my $error;
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
            ++$error;
            next;
        }
        
        foreach my $bam (@bams) {
            print $bam, "\n";
        }
    }
    die "$error errors in locating lane bams" if ($error);
}

exit;
