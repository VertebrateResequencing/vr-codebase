=head1 NAME

VertRes::Utils::Hierarchy - hierarchy utility functions

=head1 SYNOPSIS

use VertRes::Utils::Hierarchy;

my $hierarchy_util = VertRes::Utils::Hierarchy->new();

$hierarchy_util->;

=head1 DESCRIPTION

General utility functions for working on or with team145's data/mapping/release
hierarchy directory structure.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Utils::Hierarchy;

use strict;
use warnings;
use VertRes::IO;
use VertRes::Utils::FileSystem;
use File::Basename;
use File::Spec;
use File::Path;
use File::Copy;
use Cwd qw(abs_path cwd);
use VertRes::Parser::sequence_index;
use VertRes::Wrapper::samtools;
use VertRes::Parser::sam;
use VRTrack::VRTrack;
use VRTrack::Lane;
use VRTrack::Library;
use VRTrack::Sample;
use VRTrack::Project;
use VRTrack::Study;
use VRTrack::History;

use base qw(VertRes::Base);

our %platform_aliases = (ILLUMINA => 'SLX',
                         Illumina => 'SLX',
                         LS454 => '454');

our $DEFAULT_DB_SETTINGS = {host => $ENV{VRTRACK_HOST},
                            port => $ENV{VRTRACK_PORT},
                            user => $ENV{VRTRACK_RO_USER},
                            database => 'g1k_meta'};

our $nfs_disc_basename = '/nfs/vertreseq';


=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Utils::Hierarchy->new();
 Function: Create a new VertRes::Utils::Hierarchy object.
 Returns : VertRes::Utils::Hierarchy object
 Args    : n/a

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    return $self;
}

=head2 parse_lane

 Title   : parse_lane
 Usage   : my %info = $obj->parse_lane('/path/to/lane');
 Function: Extract information about a lane based on its location in the
           hierarchy directory structure.
 Returns : hash with keys study, sample, platform, library and lane.
           (where 'study' is actually the population and "analysis group", eg.
           "CEU_low_coverage")
 Args    : a directory path

=cut

sub parse_lane {
    my ($self, $lane_path) = @_;
    
    my @dirs = File::Spec->splitdir($lane_path);
    @dirs >= 5 || $self->throw("lane path '$lane_path' wasn't valid");
    
    my ($study, $sample, $platform, $library, $lane) = @dirs[-5..-1];
    
    return (study => $study, sample => $sample, platform => $platform,
            library => $library, lane => $lane);
}

=head2 lane_info

 Title   : lane_info
 Usage   : my $path = $obj->lane_info('lane_name');
 Function: Get information about a lane from the VRTrack meta database.
 Returns : hash of information, with keys:
           hierarchy_path => string,
           study          => string, (the true project code)
           project        => string, (may not be the true project code)
           sample         => string,
           individual     => string,
           individual_acc => string,
           individual_coverage => float, (the coverage of this lane's individual)
           population     => string,
           technology     => string, (aka platform, the way DCC puts it, eg.
                                      'ILLUMINA' instead of 'SLX')
           seq_tech       => string, (aka platform, the way Sanger puts it, eg.
                                      'SLX' instead of 'ILLUMINA')
           library        => string, (the hierarchy name, which is most likely
                                      similar to the true original library name)
           library_raw    => string, (the name stored in the database, which may
                                      be a uniquified version of the original
                                      library name)
           library_true   => string, (an attempt at getting the true original
                                      library name, as it was before it was
                                      munged in various ways to create library
                                      and library_raw)
           lane           => string, (aka read group)
           centre         => string, (the sequencing centre name)
           insert_size    => int, (can be undef if this lane is single-ended)
           withdrawn      => boolean,
           imported       => boolean,
           mapped         => boolean,
           vrlane         => VRTrack::Lane object
           (returns undef if lane name isn't in the database)
 Args    : lane name (read group) OR a VRTrack::Lane object.
           Optionally, a hash with key db OR vrtrack to provide the database
           connection info (defaults depend on the VRTRACK_* environment
           variables, as per VertRes::Utils::VRTrackFactory):
           db => {
            host => 'xxx',
            port => xxx,
            user => 'xxx',
            password => undef,
            database => 'g1k_meta'
           }
           -or-
           vrtrack => VRTrack::VRTrack object

           optionally, pre_swap => 1 to get info applicable to the lane in its
           state immediately prior to the last time is_processed('swapped', 1)
           was called on it.

           optionally, the optional args understood by individual_coverage() to
           configure how individual_coverage will be calculated. Or supply the
           special 'no_coverage => 1' option to not calculate
           individual_coverage, which can be quite slow.

=cut

sub lane_info {
    my ($self, $lane, %args) = @_;
    
    my $hist = VRTrack::History->new();
    my $orig_time_travel = $hist->time_travel;
    
    my ($rg, $vrlane, $vrtrack);
    if (ref($lane) && $lane->isa('VRTrack::Lane')) {
        $vrlane = $lane;
        $vrtrack = $vrlane->vrtrack;
        $rg = $vrlane->hierarchy_name;
        $lane = $rg;
    }
    else {
        if ($args{vrtrack}) {
            $vrtrack = $args{vrtrack};
        }
        else {
            my $db = $args{db} || $DEFAULT_DB_SETTINGS;
            $vrtrack = VRTrack::VRTrack->new($db);
        }
        
        $vrlane = VRTrack::Lane->new_by_hierarchy_name($vrtrack, $lane);
        $rg = $lane;
    }
    
    return unless ($rg && $vrlane && $vrtrack);
    
    my $datetime = 'latest';
    if ($args{pre_swap}) {
        $datetime = $hist->was_processed($vrlane, 'swapped');
    }
    # make sure we've got a lane of the correct time period
    $hist->time_travel($datetime);
    $vrlane = VRTrack::Lane->new_by_hierarchy_name($vrtrack, $lane) || $self->throw("Could not get a vrlane with name $lane prior to $datetime");
    
    my %info = (lane => $rg, vrlane => $vrlane);
    
    $info{hierarchy_path} = $vrtrack->hierarchy_path_of_lane($vrlane);
    $info{withdrawn} = $vrlane->is_withdrawn;
    $info{imported} = $vrlane->is_processed('import');
    $info{mapped} = $vrlane->is_processed('mapped');
    
    my %objs = $self->lane_hierarchy_objects($vrlane);
    
    $info{insert_size} = $objs{library}->insert_size;
    $info{library} = $objs{library}->hierarchy_name || $self->throw("library hierarchy_name wasn't known for $rg");
    my $lib_name = $objs{library}->name || $self->throw("library name wasn't known for $rg");
    $info{library_raw} = $lib_name;
    ($lib_name) = split(/\|/, $lib_name);
    $info{library_true} = $lib_name;
    $info{centre} = $objs{centre}->name || $self->throw("sequencing centre wasn't known for $rg");
    my $seq_tech = $objs{platform}->name || $self->throw("sequencing platform wasn't known for $rg");
    $info{seq_tech} = $seq_tech;
    if ($seq_tech =~ /illumina|slx/i) {
        $info{technology} = 'ILLUMINA';
    }
    elsif ($seq_tech =~ /solid/i) {
        $info{technology} = 'ABI_SOLID';
    }
    elsif ($seq_tech =~ /454/) {
        $info{technology} = 'LS454';
    }
    $info{sample} = $objs{sample}->name || $self->throw("sample name wasn't known for $rg");
    $info{individual} = $objs{individual}->name || $self->throw("individual name wasn't known for $rg");
    $info{individual_acc} = $objs{individual}->acc; # || $self->throw("sample accession wasn't known for $rg");
    unless ($args{no_coverage}) {
        $info{individual_coverage} = $self->hierarchy_coverage(individual => [$info{individual}],
                                                               vrtrack => $vrtrack,
                                                               $args{genome_size} ? (genome_size => $args{genome_size}) : (),
                                                               $args{gt_confirmed} ? (gt_confirmed => $args{gt_confirmed}) : (),
                                                               $args{qc_passed} ? (qc_passed => $args{qc_passed}) : (),
                                                               $args{mapped} ? (mapped => $args{mapped}) : ());
    }
    $info{population} = $objs{population}->name;
    $info{project} = $objs{project}->name;
    $info{study} = $objs{study} ? $objs{study}->acc : $info{project};
    
    $hist->time_travel($orig_time_travel);
    
    return %info;
}

=head2 lane_hierarchy_objects

 Title   : lane_hierarchy_objects
 Usage   : my %objects = $obj->lane_hierarchy_objects($lane);
 Function: Get all the parent objects of a lane, from the library up to the
           project.
 Returns : hash with these key and value pairs:
           study => VRTrack::Study object
           project => VRTrack::Project object
           sample => VRTrack::Sample object
           individual => VRTrack::Individual object
           population => VRTrack::Population object
           platform => VRTrack::Seq_tech object
           centre => VRTrack::Seq_centre object
           library => VRTrack::Library object
 Args    : VRTrack::Lane object

=cut

sub lane_hierarchy_objects {
    my ($self, $vrlane) = @_;
    
    my $vrtrack = $vrlane->vrtrack;
    my $lib = VRTrack::Library->new($vrtrack, $vrlane->library_id);
    my $sc = $lib->seq_centre;
    my $st = $lib->seq_tech;
    my $sample = VRTrack::Sample->new($vrtrack, $lib->sample_id);
    my $individual = $sample->individual;
    my $pop = $individual->population;
    my $project_obj = VRTrack::Project->new($vrtrack, $sample->project_id);
    my $study_obj = VRTrack::Study->new($vrtrack, $project_obj->study_id) if $project_obj->study_id;
    
    return (study => $study_obj,
            project => $project_obj,
            sample => $sample,
            individual => $individual,
            population => $pop,
            platform => $st,
            centre => $sc,
            library => $lib);
}

=head2 hierarchy_coverage

 Title   : hierarchy_coverage
 Usage   : my $coverage = $obj->hierarchy_coverage(sample => ['NA19239'],
                                                   genome_size => 3e9);
 Function: Discover the sequencing coverage calculated over certain lanes.
 Returns : float
 Args    : At least one hierarchy level as a key, and an array ref of names
           as values, eg. sample => ['NA19239'], platform => ['SLX', '454'].
           Valid key levels are project, sample, individual, population,
           platform, centre and library. (With no options at all, coverage will
           be calculated over all lanes in the database)
           -OR-
           A special mode can be activated by supplying a single lane name with
           the key lane, and a desired level with the level key, eg.:
           lane => 'lane_name', level => 'individual'. This would calculate the
           coverage of all the lanes that belong to the individual that the
           supplied lane belonged to.

           plus optional hash:
           genome_size => int (total genome size in bp; default 3e9)
           gt_confirmed => boolean (only consider genotype confirmed lanes;
                                    default false)
           qc_passed => boolean (only consider qc passed lanes; default false)
           mapped => boolean (coverage of mapped bases; default false: coverage
                              of total bases)

           Optionally, a hash with key db OR vrtrack to provide the database
           connection info (defaults depend on the VRTRACK_* environment
           variables, as per VertRes::Utils::VRTrackFactory):
           db => {
            host => 'xxx',
            port => xxx,
            user => 'xxx',
            password => undef,
            database => 'g1k_meta'
           }
           -or-
           vrtrack => VRTrack::VRTrack object

=cut

sub hierarchy_coverage {
    my ($self, %args) = @_;
    my $genome_size = delete $args{genome_size} || 3e9;
    my $gt = delete $args{gt_confirmed} ? 1 : 0;
    my $qc = delete $args{qc_passed} ? 1 : 0;
    my $mapped = delete $args{mapped} ? 1 : 0;
    
    if (exists $args{lane} || exists $args{level}) {
        my $lane = delete $args{lane};
        my $level = delete $args{level};
        $self->throw("Both lane and level options must be supplied if either of them are") unless $lane && $level;
        
        my @levels = qw(project sample individual population platform centre library);
        foreach my $valid_level (@levels) {
            $self->throw("'$valid_level' option is mutually exclusive of lane&level") if exists $args{$valid_level};
        }
        my %levels = map { $_ => 1 } @levels;
        $self->throw("Supplied level '$level' wasn't valid") unless exists $levels{$level};
        
        my $db = $args{db} || $DEFAULT_DB_SETTINGS;
        my $vrtrack = VRTrack::VRTrack->new($db);
        my $vrlane = VRTrack::Lane->new_by_name($vrtrack, $lane) || $self->throw("Could not get a lane from the db with name '$lane'");
        
        my %objs = $self->lane_hierarchy_objects($vrlane);
        $self->throw("Could not get the $level of lane $lane") unless defined $objs{$level};
        
        $args{$level} = [$objs{$level}->name];
    }
    
    my @store = ($gt, $qc, $mapped);
    while (my ($key, $val) = each %args) {
        unless (ref($val)) {
            push(@store, $val);
        }
        else {
            if (ref($val) eq 'ARRAY') {
                push(@store, @{$val});
            }
            elsif (ref($val) eq 'HASH') {
                while (my ($sub_key, $sub_val) = each %{$val}) {
                    push(@store, $sub_val);
                }
            }
            elsif ($val->isa('VRTrack::VRTrack')) {
                my $db_params = $val->database_params;
                while (my ($sub_key, $sub_val) = each %{$db_params}) {
                    push(@store, $sub_val);
                }
            }
        }
    }
    my $store = join(",", sort @store);
    
    unless (defined $self->{_cover_bases}->{$store}) {
        my @lanes = $self->get_lanes(%args);
        @lanes || return 0;
        my $bps = 0;
        
        # sum raw bases for all the qc passed, gt confirmed and not withdrawn
        # lanes
        foreach my $lane (@lanes) {
            next if $lane->is_withdrawn;
            if ($gt) {
                next unless ($lane->genotype_status && $lane->genotype_status eq 'confirmed');
            }
            if ($qc) {
                next unless ($lane->qc_status && $lane->qc_status eq 'passed');
            }
            
            my $bp = $lane->raw_bases || 0;
            
            if ($mapped) {
                my $mapstats = $lane->latest_mapping;
                
                if ($mapstats && $mapstats->raw_bases){
                    if ($mapstats->genotype_ratio) {
                        # this is a QC mapped lane, so we make a projection
                        $bps += $bp * ($mapstats->rmdup_bases_mapped / $mapstats->raw_bases);
                    }
                    else {
                        # this is a fully mapped lane, so we know the real answer
                        $bps += $mapstats->bases_mapped;
                    }
                }
                else {
                    $bps += $bp * 0.9; # not sure what else to do here?
                }
            }
            else {
                $bps += $bp;
            }
        }
        
        $self->{_cover_bases}->{$store} = $bps;
    }
    
    return sprintf('%.2f', $self->{_cover_bases}->{$store} / $genome_size);
}

=head2 get_lanes

 Title   : get_lanes
 Usage   : my @lanes = $obj->get_lanes(sample => ['NA19239']);
 Function: Get all the lanes under certain parts of the hierarchy, excluding
           withdrawn lanes.
 Returns : list of VRTrack::Lane objects
 Args    : At least one hierarchy level as a key, and an array ref of names
           as values, eg. sample => ['NA19239'], platform => ['SLX', '454'].
           Valid key levels are project, sample, individual, population,
           platform, centre and library. (With no options at all, all active
           lanes in the database will be returned)

           Optionally, a hash with key db OR vrtrack to provide the database
           connection info (defaults depend on the VRTRACK_* environment
           variables, as per VertRes::Utils::VRTrackFactory):
           db => {
            host => 'xxx',
            port => xxx,
            user => 'xxx',
            password => undef,
            database => 'g1k_meta'
           }
           -or-
           vrtrack => VRTrack::VRTrack object

=cut

sub get_lanes {
    my ($self, %args) = @_;
    
    my $vrtrack;
    if ($args{vrtrack}) {
        $vrtrack = $args{vrtrack};
    }
    else {
        my $db = $args{db} || $DEFAULT_DB_SETTINGS;
        $vrtrack = VRTrack::VRTrack->new($db);
    }
    
    my @good_lanes;
    foreach my $project (@{$vrtrack->projects}) {
        my $ok = 1;
        if (defined ($args{project})) {
            $ok = 0;
            foreach my $name (@{$args{project}}) {
                if ($name eq $project->name || $name eq $project->hierarchy_name || $name eq $project->study->acc) {
                    $ok = 1;
                    last;
                }
            }
        }
        $ok || next;
        
        foreach my $sample (@{$project->samples}) {
            my $ok = 1;
            if (defined ($args{sample})) {
                $ok = 0;
                foreach my $name (@{$args{sample}}) {
                    if ($name eq $sample->name) {
                        $ok = 1;
                        last;
                    }
                }
            }
            $ok || next;
            
            my %objs;
            $objs{individual} = $sample->individual;
            $objs{population} = $objs{individual}->population;
            
            my ($oks, $limits) = (0, 0);
            foreach my $limit (qw(individual population)) {
                if (defined $args{$limit}) {
                    $limits++;
                    my $ok = 0;
                    foreach my $name (@{$args{$limit}}) {
                        if ($name eq $objs{$limit}->name || ($objs{$limit}->can('hierarchy_name') && $name eq $objs{$limit}->hierarchy_name)) {
                            $ok = 1;
                            last;
                        }
                    }
                    $oks += $ok;
                }
            }
            next unless $oks == $limits;
            
            foreach my $library (@{$sample->libraries}) {
                my $ok = 1;
                if (defined ($args{library})) {
                    $ok = 0;
                    foreach my $name (@{$args{library}}) {
                        if ($name eq $library->name || $name eq $library->hierarchy_name) {
                            $ok = 1;
                            last;
                        }
                    }
                }
                $ok || next;
                
                my %objs;
                $objs{centre} = $library->seq_centre;
                $objs{platform} = $library->seq_tech;
                
                my ($oks, $limits) = (0, 0);
                foreach my $limit (qw(centre platform)) {
                    if (defined $args{$limit}) {
                        $limits++;
                        my $ok = 0;
                        foreach my $name (@{$args{$limit}}) {
                            if ($name eq $objs{$limit}->name) {
                                $ok = 1;
                                last;
                            }
                        }
                        $oks += $ok;
                    }
                }
                next unless $oks == $limits;
                
                push(@good_lanes, @{$library->lanes});
            }
        }
    }
    
    # filter out withdrawn lanes
    my @active;
    foreach my $lane (@good_lanes) {
        next if $lane->is_withdrawn;
        push(@active, $lane);
    }

    return @active;
}

=head2 check_lanes_vs_database

 Title   : check_lanes_vs_database
 Usage   : my $ok = $obj->check_lanes_vs_database(['/lane/paths', ...],
                                                  $vrtrack);
 Function: Check that the given lanes reside in the correct part of the
           hierarchy by checking the information in the database.
 Returns : boolean (true if all lanes agree with the database)
 Args    : reference to a list of lane paths to check, VRTrack::VRTrack object,
           boolean to turn on optional checking to see if you're missing any
           lanes from any of the samples that your input lanes belong to

=cut

sub check_lanes_vs_database {
    my ($self, $lanes, $vrtrack, $check_for_missing) = @_;
    
    my $all_ok = 1;
    my @lane_ids;
    my @sample_names;
    foreach my $lane_path (@{$lanes}) {
        my %lane_info = $self->parse_lane($lane_path);
        my $lane_id = basename($lane_path);
        my $vrlane = VRTrack::Lane->new_by_hierarchy_name($vrtrack, $lane_id);
        
        # check this lane is even in the database
        unless ($vrlane) {
            $all_ok = 0;
            $self->warn("not in db: $lane_id ($lane_path)");
            next;
        }
        
        # check this lane hasn't been withdrawn
        if ($vrlane->is_withdrawn) {
            $self->warn("withdrawn: $lane_id ($lane_path)");
            $all_ok = 0;
            next;
        }
        
        push(@lane_ids, $lane_id);
        my %objs = $self->lane_hierarchy_objects($vrlane);
        push(@sample_names, $objs{sample}->name);
        
        
        my $sample_name = $objs{individual}->hierarchy_name;
        my $platform = $objs{platform}->name;
        if (exists $platform_aliases{$platform}) {
            $platform = $platform_aliases{$platform};
        }
        my $library = $objs{library}->hierarchy_name;
        
        # quick test against the hierarchy_path
        my $expected_path = $vrtrack->hierarchy_path_of_lane($vrlane);
        next if $lane_path =~ /$expected_path$/;
        
        # study swaps
        my $given_study = $lane_info{study};
        my $study = $objs{project}->hierarchy_name;
        unless ($study eq $given_study) {
            $self->warn("study swap: $study vs $given_study for $lane_id ($lane_path -> $expected_path)");
            $all_ok = 0;
        }
        
        # sample swaps
        unless ($sample_name eq $lane_info{sample}) {
            $self->warn("sample swap: $sample_name vs $lane_info{sample} for $lane_id ($lane_path -> $expected_path)");
            $all_ok = 0;
        }
        
        # platform swaps
        unless ($platform eq $lane_info{platform}) {
            $self->warn("platform swap: $platform vs $lane_info{platform} for $lane_id ($lane_path -> $expected_path)");
            $all_ok = 0;
        }
        
        # library swaps
        unless ($library eq $lane_info{library}) {
            $self->warn("library swap: $library vs $lane_info{library} for $lane_id ($lane_path -> $expected_path)");
            $all_ok = 0;
        }
    }
    
    if ($check_for_missing) {
        my @expected_lanes = $self->get_lanes(sample => \@sample_names);
        
        # uniquify
        my %expected_lanes = map { $_->hierarchy_name => 1 } @expected_lanes;
        @expected_lanes = sort keys %expected_lanes;
        
        my %actual_lanes = map { $_ => 1 } @lane_ids;
        
        foreach my $lane (@expected_lanes) {
            unless (exists $actual_lanes{$lane}) {
                $self->warn("missing: $lane was in the database but not in the supplied list of lanes");
                $all_ok = 0;
            }
        }
    }
    
    return $all_ok;
}

=head2 check_lanes_vs_sequence_index

 Title   : check_lanes_vs_sequence_index
 Usage   : my $ok = $obj->check_lanes_vs_sequence_index(['/lane/paths', ...],
                                                        'sequence.index');
 Function: Check that the given lanes reside in the correct part of the
           hierarchy by checking the information in the sequence.index file.
 Returns : boolean (true if all lanes agree with the sequence index)
 Args    : reference to a list of lane paths to check, sequence.index filename,
           boolean to turn on optional checking to see if you're missing any
           lanes that you should have according to the sequence.index

=cut

sub check_lanes_vs_sequence_index {
    my ($self, $lanes, $sequence_index, $check_for_missing) = @_;
    
    my $sip = VertRes::Parser::sequence_index->new(file => $sequence_index,
                                                   verbose => $self->verbose);
    
    my $all_ok = 1;
    my @lane_ids;
    foreach my $lane_path (@{$lanes}) {
        my %lane_info = $self->parse_lane($lane_path);
        my $lane_id = $lane_info{lane};
        push(@lane_ids, $lane_id);
        
        my $sample_name = $sip->lane_info($lane_id, 'sample_name');
        my $platform = $sip->lane_info($lane_id, 'INSTRUMENT_PLATFORM');
        if (exists $platform_aliases{$platform}) {
            $platform = $platform_aliases{$platform};
        }
        my $library = $sip->lane_info($lane_id, 'LIBRARY_NAME');
        $library =~ s/\s/_/;
        my $expected_path = join('/', $sample_name, $platform, $library);
        
        # check this lane is even in the sequence.index; $sip will warn if not
        unless ($sample_name) {
            $all_ok = 0;
            next;
        }
        
        # check this lane hasn't been withdrawn
        my $withdrawn = $sip->lane_info($lane_id, 'WITHDRAWN');
        if ($withdrawn) {
            $self->warn("withdrawn: $lane_id ($lane_path)");
            $all_ok = 0;
            next;
        }
        
        # study swaps
        my $given_study = $lane_info{study};
        my $study = $sip->lane_info($lane_id, 'population').'_'.$sip->lane_info($lane_id, 'ANALYSIS_GROUP');
        $study =~ s/\s/_/g;
        unless ($study eq $given_study) {
            $self->warn("study swap: $study vs $given_study for $lane_id ($lane_path -> $expected_path)");
            $all_ok = 0;
        }
        
        # sample swaps
        unless ($sample_name eq $lane_info{sample}) {
            $self->warn("sample swap: $sample_name vs $lane_info{sample} for $lane_id ($lane_path -> $expected_path)");
            $all_ok = 0;
        }
        
        # platform swaps
        unless ($platform eq $lane_info{platform}) {
            $self->warn("platform swap: $platform vs $lane_info{platform} for $lane_id ($lane_path -> $expected_path)");
            $all_ok = 0;
        }
        
        # library swaps
        unless ($library eq $lane_info{library}) {
            $self->warn("library swap: $library vs $lane_info{library} for $lane_id ($lane_path -> $expected_path)");
            $all_ok = 0;
        }
    }
    
    if ($check_for_missing) {
        my @lanes = $sip->get_lanes(ignore_withdrawn => 1);
        
        # if we only ignore lines with withdrawn, it doesn't stop us picking up
        # lanes that were both withdrawn and not withdrawn. Filter those out as
        # well:
        my @expected_lanes;
        foreach my $lane (@lanes) {
            my $withdrawn = $sip->lane_info($lane, 'withdrawn');
            if ($withdrawn) {
                $self->warn("lane $lane was both withdrawn and not withdrawn - treating it as withdrawn");
            }
            else {
                push(@expected_lanes, $lane);
            }
        }
        
        # uniquify
        my %expected_lanes = map { $_ => 1 } @expected_lanes;
        @expected_lanes = sort keys %expected_lanes;
        
        my %actual_lanes = map { $_ => 1 } @lane_ids;
        
        foreach my $lane (@expected_lanes) {
            unless (exists $actual_lanes{$lane}) {
                $self->warn("missing: $lane was in the sequence.index but not in the supplied list of lanes");
                $all_ok = 0;
            }
        }
    }
    
    return $all_ok;
}

=head2 fix_simple_swaps

 Title   : fix_simple_swaps
 Usage   : my $swapped = $obj->fix_simple_swaps(['/lane/paths', ...],
                                                'sequence.index');
 Function: For lanes that check_lanes_vs_sequence_index() would complain
           suffered from a swap, moves the lane to the correct part of the
           hierarchy.
 Returns : int (the number of swaps fixed)
 Args    : reference to a list of lane paths to check, sequence.index filename

=cut

sub fix_simple_swaps {
    my ($self, $lanes, $sequence_index) = @_;
    
    my $sip = VertRes::Parser::sequence_index->new(file => $sequence_index,
                                                   verbose => 0);
    
    my $fixed = 0;
    my @lane_ids;
    foreach my $lane_path (@{$lanes}) {
        my %lane_info = $self->parse_lane($lane_path);
        my $lane_id = $lane_info{lane};
        push(@lane_ids, $lane_id);
        
        my $sample_name = $sip->lane_info($lane_id, 'sample_name');
        next unless $sample_name;
        next if $sip->lane_info($lane_id, 'WITHDRAWN');
        
        my $platform = $sip->lane_info($lane_id, 'INSTRUMENT_PLATFORM');
        if (exists $platform_aliases{$platform}) {
            $platform = $platform_aliases{$platform};
        }
        my $library = $sip->lane_info($lane_id, 'LIBRARY_NAME');
        $library =~ s/\s/_/;
        
        # can't handle study swaps, because we don't know which study directory
        # we'd need to move to
        
        # sample swaps
        unless ($sample_name eq $lane_info{sample}) {
            $fixed += $self->_fix_simple_swap($lane_path, $sample_name, $platform, $library);
            next;
        }
        
        # platform swaps
        unless ($platform eq $lane_info{platform}) {
            $fixed += $self->_fix_simple_swap($lane_path, $sample_name, $platform, $library);
            next;
        }
        
        # library swaps
        unless ($library eq $lane_info{library}) {
            $fixed += $self->_fix_simple_swap($lane_path, $sample_name, $platform, $library);
            next;
        }
    }
    
    return $fixed;
}

sub _fix_simple_swap {
    my ($self, $old, $new_sample, $new_platform, $new_library) = @_;
    
    $old = abs_path($old);
    my @dirs = File::Spec->splitdir($old);
    @dirs >= 5 || $self->throw("lane path '$old' wasn't valid");
    my ($sample, $platform, $library, $lane) = splice(@dirs, -4);
    
    # do the swap by moving the lane dir
    my $new = File::Spec->catdir(@dirs, $new_sample, $new_platform, $new_library, $lane);
    if (-d $new) {
        $self->warn("Wanted to swap $old -> $new, but the destination already exists");
        return 0;
    }
    else {
        my $parent = File::Spec->catdir(@dirs, $new_sample, $new_platform, $new_library);
        mkpath($parent);
        $self->throw("Could not create path $parent") unless -d $parent;
    }
    move($old, $new);
    
    # did we just empty any of the parent directories? if so, remove them
    my $parent = File::Spec->catdir(@dirs, $sample, $platform, $library);
    $self->_remove_empty_parent($parent);
    $parent = File::Spec->catdir(@dirs, $sample, $platform);
    $self->_remove_empty_parent($parent);
    $parent = File::Spec->catdir(@dirs, $sample);
    $self->_remove_empty_parent($parent);
    
    return 1;
}

sub _remove_empty_parent {
    my ($self, $parent) = @_;
    
    opendir(my $pfh, $parent) || $self->throw("Could not open dir $parent");
    my $things = 0;
    foreach (readdir($pfh)) {
        next if /^\.+$/;
        $things++;
    }
    
    if ($things == 0) {
        system("rm -fr $parent");
    }
}

=head2 create_release_hierarchy

 Title   : create_release_hierarchy
 Usage   : my $ok = $obj->create_release_hierarchy(\@abs_lane_paths,
                                                   '/path/to/release',
                                                   vrtrack => $vrtrack,
                                                   slx_mapper => 'bwa',
                                                   '454_mapper' => 'ssaha',
                                                   assembly_name => 'NCBI37');
 Function: Given a list of absolute paths to mapped lanes in a mapping
           hierarchy, creates a release hierarchy with mapped bams symlinked
           across.
           NB: You should probably call check_lanes_vs_database() on the
           the input lanes beforehand.
 Returns : list of absolute paths to the bam symlinks in the created release
           hierarchy
 Args    : reference to a list of mapped lane paths, base path you want to
           build the release hierarchy in, and a hash of the following required
           options:
           vrtrack => VRTrack::VRTrack object
           slx_mapper => string
           454_mapper => string
           assembly_name => string

=cut

sub create_release_hierarchy {
    my ($self, $lane_paths, $release_dir, %args) = @_;
    
    unless (-d $release_dir) {
        mkdir($release_dir) || $self->throw("Unable to create base release directory: $!");
    }
    my $io = VertRes::IO->new();
    my $fsu = VertRes::Utils::FileSystem->new();
    
    my @all_linked_bams;
    
    foreach my $lane_path (@{$lane_paths}) {
        # first get the mapstats object so we'll know the bam prefix
        my $lane_name = basename($lane_path);
        my $vrlane = VRTrack::Lane->new_by_hierarchy_name($args{vrtrack}, $lane_name) || $self->throw("Unable to get lane $lane_name from db");
        my %objs = $self->lane_hierarchy_objects($vrlane);
        my $platform = lc($objs{platform}->name);
        my $mappings = $vrlane->mappings();
        my $mapstats;
        if ($mappings && @{$mappings}) {
            # find the most recent mapstats that corresponds to our mapping
            my $highest_id = 0;
            foreach my $possible (@{$mappings}) {
                # we're expecting it to have the correct assembly and mapper
                my $assembly = $possible->assembly() || next;
                $assembly->name eq $args{assembly_name} || next;
                my $mapper = $possible->mapper() || next;
                $mapper->name eq $args{$platform.'_mapper'} || next;
                
                if ($possible->id > $highest_id) {
                    $mapstats = $possible;
                    $highest_id = $possible->id;
                }
            }
        }
        $mapstats || $self->throw("Could not get a mapstats for lane $lane_name");
        my $mapstats_prefix = $mapstats->id;
        
        # setup release lane
        my @dirs = File::Spec->splitdir($lane_path);
        @dirs >= 5 || $self->throw("lane path '$lane_path' wasn't valid");
        my $release_lane_path = $fsu->catfile($release_dir, @dirs[-5..-1]);
        mkpath($release_lane_path);
        -d $release_lane_path || $self->throw("Unable to make release lane '$release_lane_path'");
        
        # symlink bams; what ended are we supposed to have?
        my $files = $vrlane->files();
        my %ended;
        foreach my $file (@{$files}) {
            my $type = $file->type;
            if (! $type) {
                $ended{se} = 1;
            }
            else {
                $ended{pe} = 1;
            }
        }
        
        my @linked_bams;
        foreach my $ended (keys %ended) {
            # prefer the recalibrated bam if it was made to the unrecalibrated
            # raw bam
            my $type = 'recal';
            my $source = $fsu->catfile($lane_path, "$mapstats_prefix.$ended.$type.sorted.bam");
            unless (-s $source) {
                $type = 'raw';
                $source = $fsu->catfile($lane_path, "$mapstats_prefix.$ended.$type.sorted.bam");
            }
            
            if (-s $source) {
                my $destination = $fsu->catfile($release_lane_path, "$ended.$type.sorted.bam");
                symlink($source, $destination) || $self->throw("Couldn't symlink $source -> $destination");
                push(@linked_bams, $destination);
            }
            else {
                $self->warn("expected there to be a bam '$source' but it didn't exist!");
            }
        }
        
        # old mapping pipeline didn't create pe/se bams
        unless (@linked_bams) {
            my $legacy_bam_file = 'raw.sorted.bam';
            my $source = $fsu->catfile($lane_path, $legacy_bam_file);
            if (-s $source) {
                # figure out if it was generated from paired or single ended reads
                my $meta_file = $fsu->catfile($lane_path, 'meta.info');
                if (-s $meta_file) {
                    $io->file($meta_file);
                    my $fh = $io->fh;
                    my %reads;
                    while (<$fh>) {
                        if (/^read(\d)/) {
                            $reads{$1} = 1;
                        }
                    }
                    $io->close;
                    
                    my $ended;
                    if ($reads{1} && $reads{2}) {
                        $ended = 'pe';
                    }
                    elsif ($reads{0}) {
                        $ended = 'se';
                    }
                    else {
                        $self->throw("Couldn't determine if reads were single or paired end in lane '$lane_path'");
                    }
                    
                    my $bam_file = "${ended}_raw.sorted.bam";
                    my $destination = $fsu->catfile($release_lane_path, $bam_file);
                    symlink($source, $destination) || $self->throw("Couldn't symlink $source -> $destination");
                    push(@linked_bams, $destination);
                }
            }
        }
        
        unless (@linked_bams) {
            $self->throw("Mapping lane '$lane_path' contained no linkable bam files!");
        }
        
        push(@all_linked_bams, @linked_bams);
    }
    
    return @all_linked_bams;
}

=head2 dcc_filename

 Title   : dcc_filename
 Usage   : my $filename = $obj->dcc_filename('release.bam',
                                             '20100208',
                                             'sequence.index');
 Function: Get the DCC filename of a bam file. For this to work, the bam file
           must have RG lines in the header where ID, PL, LB, PI, SM and CN are
           all set, and DS is set to the SRP project code.
 Returns : string (filename without .bam suffix)
 Args    : path to release bam file, the release name date string (YYYYMMDD
           corresponding to the sequence.index the release was made from), and
           the path to the sequence.index file. Optionally a chrom string if the
           supplied bam is unsplit, but you want to work out what the dcc
           filename of a certain chromosmal split of that bam would be prior to
           actually making the split bam.

=cut

sub dcc_filename {
    my ($self, $file, $date_string, $sequence_index, $given_chrom) = @_;
    $date_string || $self->throw("release date string must be supplied");
    -s $sequence_index || $self->throw("sequence.index file must be supplied");
    
    # NAXXXXX.[chromN].technology.[center].algorithm.population.analysis_group.YYYYMMDD.bam
    # eg. NA12878.chrom1.LS454.ssaha.CEU.high_coverage.20091216.bam
    # http://1000genomes.org/wiki/doku.php?id=1000_genomes:dcc:metadata
    
    # view the bam header
    my $stw = VertRes::Wrapper::samtools->new(quiet => 1);
    $stw->run_method('open');
    my $view_fh = $stw->view($file, undef, H => 1);
    $view_fh || $self->throw("Failed to samtools view '$file'");
    
    my $ps = VertRes::Parser::sam->new(fh => $view_fh);
    
    my ($sample, $platform) = ('unknown_sample', 'unknown_platform');
    my %techs;
    my %readgroup_info = $ps->readgroup_info();
    my $example_rg;
    while (my ($rg, $info) = each %readgroup_info) {
        # there should only be one sample, so we just keep resetting it
        $sample = $info->{SM} || 'unknown_sample';
        $example_rg = $rg;
        
        # might be more than one of these if we're a sample-level bam. We
        # standardise on the DCC nomenclature for the 3 platforms; they should
        # be in this form anyway, so this is just-in-case
        $platform = $info->{PL};
        if ($platform =~ /illumina|slx/i) {
            $platform = 'ILLUMINA';
        }
        elsif ($platform =~ /solid/i) {
            $platform = 'ABI_SOLID';
        }
        elsif ($platform =~ /454/) {
            $platform = 'LS454';
        }
        $techs{$platform}++;
    }
    
    # instead of the srp (project code), we now have population and analysis
    # group in it's place. These things are not stored in the bam header, so
    # we must check the sequence.index file to figure this out. We assume that
    # all the readgroups have the same pop and ag, since they DCC only do
    # same-sample bams
    my $sip = VertRes::Parser::sequence_index->new(file => $sequence_index,
                                                   verbose => $self->verbose);
    my $pop = $sip->lane_info($example_rg, 'POPULATION') || 'unknown_population';
    my $ag = $sip->lane_info($example_rg, 'ANALYSIS_GROUP') || 'unknown_analysisgroup';
    $ag =~ s/\s/_/g;
    
    # if there's more than 1 tech, it doesn't appear in the filename
    if (keys %techs > 1) {
        $platform = '';
    }
    else {
        $platform .= '.';
    }
    
    my $bamname = basename($file);
    my $chrom = '';
    if ($given_chrom) {
        if ($given_chrom =~ /^(?:\d+|[XY]|MT)$/) {
            $given_chrom = 'chrom'.$given_chrom;
        }
        $chrom = "$given_chrom.";
    }
    elsif ($bamname =~ /\.?(chrom(?:\d+|[XY]|MT)|nonchrom|unmapped)\./) {
        $chrom = "$1.";
    }
    
    my $algorithm = $ps->program || 'unknown_algorithm';
    # picard merge can fuck with program names, converting them to unique numbers
    if ($algorithm eq 'unknown_algorithm' || $algorithm =~ /^\d+$/ || $algorithm =~ /GATK/) { 
        if ($platform =~ /ILLUMINA/) {
            $algorithm = 'bwa';
        }
        elsif ($platform =~ /454/) {
            $algorithm = 'ssaha2';
        }
        elsif ($platform =~ /SOLID/) {
            $algorithm = 'bfast';
        }
    }
    
    my $dcc_filename = "$sample.$chrom$platform$algorithm.$pop.$ag.$date_string";
    
    return $dcc_filename;
}

=head2 nfs_disks

 Title   : nfs_disks
 Usage   : my @disks = $obj->nfs_disks();
 Function: Find and prepare all of team145's nfs disks intended for storing
           lane directories.
 Returns : list of disk root directories
 Args    : n/a

=cut

sub nfs_disks {
    my $self = shift;
    
    # /nfs/vertreseq01 - /nfs/vertreseq16
    my @disks;
    my $cwd = cwd();
    foreach my $i (1..99) {
        my $disk_num = sprintf("%02d", $i);
        my $disk = $nfs_disc_basename.$disk_num;
        
        # spin the disk up
        chdir($disk);
        chdir($cwd);
        
        # does it exist?
        -d $disk || last;
        push(@disks, $disk);
    }
    
    return @disks;
}

=head2 nfs_disk

 Title   : nfs_disk
 Usage   : my $dir = $obj->nfs_disk();
 Function: Get the path to the root of one of the nfs_disks() - the one with
           the most available disk space.
 Returns : string path
 Args    : n/a

=cut

sub nfs_disk {
    my $self = shift;
    my $fsu = VertRes::Utils::FileSystem->new();
    
    my @disks = $self->nfs_disks;
    my $most_available = 0;
    my $best_disk;
    foreach my $disk (@disks) {
        my $available = $fsu->disk_available($disk);
        if ($available > $most_available) {
            $most_available = $available;
            $best_disk = $disk;
        }
    }
    
    return $best_disk;
}

=head2 lane_storage_path

 Title   : lane_storage_path
 Usage   : my $path = $obj->lane_storage_path($lane);
 Function: Get the absolute path to where a lane either is or should be stored
           on an nfs disk.
 Returns : path string
 Args    : VRTrack::Lane object

=cut

sub lane_storage_path {
    my ($self, $lane) = @_;
    my $storage_path = $lane->storage_path;
    
    unless ($storage_path) {
        my $vrtrack = $lane->vrtrack;
        my $hpath = $vrtrack->hierarchy_path_of_lane($lane);
        my $db_params = $vrtrack->database_params;
        my $fsu = VertRes::Utils::FileSystem->new();
        $storage_path = $fsu->catfile($self->nfs_disk, 'hashed_lanes', $db_params->{database}, $fsu->hashed_path($hpath));
    }
    
    return $storage_path;
}

=head2 store_lane

 Title   : store_lane
 Usage   : $obj->store_lane('/abs/path/to/hierarchy/root', $lane);
 Function: Move a lane directory from its current location to where it should
           be according to lane_storage_path(). Once successfully moved a
           symlink to the storage location will be left in the original
           location, and the database will be updated with the storage_path.
 Returns : boolean (true if the move was successful or not necessary because the
           lane has already been stored on nfs; false otherwise)
 Args    : absolute path to hierarchy root containing the lane directory (or the
           full path to the lane), VRTrack::Lane object

=cut

sub store_lane {
    my ($self, $hroot, $lane) = @_;
    
    my $fsu = VertRes::Utils::FileSystem->new();
    
    my $storage_path = $self->lane_storage_path($lane);
    my $storage_path_temp = $storage_path."_store_lane_temp";
    my $do_move = 1;
    if (-d $storage_path && $lane->is_processed('stored')) {
        return 1;
    }
    elsif (-l $hroot && -d $storage_path) {
        $do_move = 0;
    }
    elsif (-d $hroot && (-d $storage_path || -d $storage_path_temp)) {
        $self->warn("storage path '$storage_path' already exists; will delete it first");
        $fsu->rmtree($storage_path);
        $fsu->rmtree($storage_path_temp);
    }
    
    my $hpath = $lane->vrtrack->hierarchy_path_of_lane($lane);
    
    my $cwd = cwd();
    my $source_dir;
    if ($do_move) {
        $hroot =~ s/$hpath//;
        $source_dir = $fsu->catfile($hroot, $hpath);
        unless (-d $source_dir) {
            $self->warn("Lane '$source_dir' wasn't a directory, can't move it to store it on NFS");
            return 0;
        }
        
        # if we're in the source_dir, move() will chdir to parent, so store our
        # current directory and chdir back to it afterwards
        my $moved = $fsu->move($source_dir, $storage_path_temp);
        
        unless ($moved) {
            $self->warn("Failed to move $source_dir to storage path '$storage_path_temp'");
            return 0;
        }
    }
    
    my $vrtrack = $lane->vrtrack;
    $vrtrack->transaction_start();
    
    if ($do_move) {
        symlink($storage_path, $source_dir) || $self->throw("Failed to create symlink from $storage_path to $source_dir");
        chdir($cwd);
        
        File::Copy::move($storage_path_temp, $storage_path) || $self->throw("Could not rename $storage_path_temp to $storage_path");
    }
    
    $lane->storage_path($storage_path);
    $lane->is_processed('stored', 1);
    $lane->update || $self->throw("Could not update db to note that lane $hpath has been stored");
    
    $vrtrack->transaction_commit();
    
    return 1;
}

1;
