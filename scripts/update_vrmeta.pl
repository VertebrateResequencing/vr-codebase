#!/usr/bin/env perl
#
# update_vrmeta.pl --g1k
#
# Adds/updates our VRTrack g1k meta database with information from the DCC
# sequence.index file. Can also update any other meta database from fake
# sequence.index files.
#
# Author: Sendu Bala <bix@sendu.me.uk>
#

use strict;
use warnings;
use File::Spec;
use File::Basename;
use Cwd 'abs_path';
use File::Path;
use Getopt::Long;
use VertRes::Parser::sequence_index;
use VertRes::Utils::VRTrackFactory;
use VRTrack::Project;
use VRTrack::Sample;
use VRTrack::Library;
use VRTrack::Lane;
use VRTrack::Individual;
use VRTrack::VRTrack;

my $schema_version = VRTrack::VRTrack::SCHEMA_VERSION;

# get user input
my ($help, $do_major_updates, $g1k_mode, $si_file, $sample_to_population_map_file, $database, $ignore_diff_md5s, $only_update_diff_md5s, $merge_projects, $update_withdrawn, $verbose, $debug);
my $spinner = 0;
GetOptions('g1k'         => \$g1k_mode,
           'samples=s'   => \$sample_to_population_map_file,
           'index=s'     => \$si_file,
           'database=s'  => \$database,
           'merge_projects' => \$merge_projects,
           'ignore_changed_fastqs' => \$ignore_diff_md5s,
           'only_update_changed_fastqs' => \$only_update_diff_md5s,
           'all_updates' => \$do_major_updates,
           'update_withdrawn' => \$update_withdrawn,
           'status'      => \$spinner,
           'verbose'     => \$verbose,
           'debug=s'     => \$debug,
           'h|help'      => \$help);

if ($g1k_mode) {
    $si_file ||= 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/sequence.index';
    $sample_to_population_map_file ||= '/lustre/scratch105/conf/g1k_samples';
    $database ||= 'g1k_meta';
    $merge_projects = 1;
}

my $missing_opts = 0;
unless ($si_file && $database) {
    $missing_opts = 1;
}

if ($debug) {
    $do_major_updates = 0;
}

($help || $missing_opts) and die <<USAGE;
Updates a VRTrack meta database from a sequence.index file.

Usage: $0 --samples sample.info --index sequence.index --database my_vrtrack_meta
        --g1k         sets samples, index, database to values suitable for
                      the main g1k project; merge_projects is turned on
        --samples     <file of sample info>
        --index       <sequence.index file>
        --database    <database name>
        --merge_projects (when specified, different projects that share the same
                          project style (eg. low coverage) and population will
                          have their disc locations specified the same)
        --ignore_changed_fastqs (when a fastq seems to have changed because the
                                 md5 is different comparing sequence.index and
                                 that already stored in the db for a file of a
                                 given name, don't worry about it. Useful for
                                 tracking -> meta database conversions, where
                                 the tracking md5 is for the uncompressed fastq)
        --only_update_changed_fastqs (when a fastq seems to have changed because
                                      the md5 is different, update the md5 in
                                      the db, but assume that the fastq data
                                      didn't really change - so don't trigger a
                                      reimport and remapping of this fastq)
        --all_updates (update lanes even when there are major changes,
                       like sample swaps and outright lane deletions)
        --update_withdrawn (in combination with --all_updates, will update major
                            changes to withdrawn lanes, which might be dangerous
                            if the details for withdrawn lanes are wrong and
                            clash with active lanes)
        --status      print a spinner to indicate activity
        --verbose     be extra verbose about what it is doing
        --help        <this message>

Sample info file is tab-delimited with columns:
lookup name (used to match the sequencescape sample name)
acc
individual name
alias (used for genotyping)
population name
species name
taxon id
sex

Sample info file is optional if all the necessary information has already been
added to the database with eg. load_vrtrack_allocations.pl

USAGE

my %platform_to_tech = (ILLUMINA => 'SLX',
                        LS454 => '454',
                        ABI_SOLID => 'SOLID',
                        SLX => 'SLX',
                        '454' => '454',
                        SOLID => 'SOLID');

my %sample_data;
if ($sample_to_population_map_file) {
    open(my $stpmfh, $sample_to_population_map_file) || die "Could not open $sample_to_population_map_file\n";
    while (<$stpmfh>) {
        chomp;
        # NA19062 SRS000752       NA19062 NA19062 JPT     Homo sapiens    9606    M
        my ($lookup, $acc, $ind, $alias, $pop, $spp, $taxon, $sex) = split("\t", $_);
        $lookup = uc($lookup);
        $sample_data{$lookup}->{acc} = $acc;
        $sample_data{$lookup}->{ind} = $ind;
        $sample_data{$lookup}->{alias} = $alias;
        $sample_data{$lookup}->{pop} = $pop;
        $sample_data{$lookup}->{spp} = $spp;
        $sample_data{$lookup}->{taxon} = $taxon;
        $sample_data{$lookup}->{sex} = uc($sex);
    }
    close($stpmfh);
}

my %object_cache;

# setup database connection
my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $database,
                                                          mode => 'rw');
unless ($vrtrack) {
    die "DB connection failed: ".$DBI::errstr."\n";;
}

# setup read of the sequence.index
my $sip = VertRes::Parser::sequence_index->new(file => $si_file, verbose => $verbose);
my $result_holder = $sip->result_holder;

my $new_lanes_count = 0;
my %si_lanes;
my %withdrawn_statuses;
my %paired_statuses;
my %issued_warnings;
my $spin = 0;
local $| = 1;
print STDERR " " if $spinner;
my $si_lines = 0;
my $ok_si_lines = 0;

# parse through the sequence.index, file by file, creating/updating lanes
while ($sip->next_result) {
    main_loop($result_holder);
}

# did the lane withdrawn and paired status really change?
while (my ($lane_name, $file_status) = each %withdrawn_statuses) {
    my $lane = VRTrack::Lane->new_by_name($vrtrack, $lane_name);
    
    my $changed = 0;
    
    # withdrawn
    my $current_status = $lane->is_withdrawn();
    
    # we can make a lane withdrawn, but we don't unwithdraw lanes;
    # a lane can have multiple withdrawn states in the same
    # sequence.index file, and we withdraw if any of its fastqs are
    # withdrawn
    my $new_status = 0;
    while (my ($file, $status) = each %{$file_status}) {
        $new_status = 1 if $status;
    }
    
    if (defined $current_status) {
        if ($current_status == 0 && $new_status == 1) {
            issue_warning("Lane $lane_name just got withdrawn");
            # we allow lanes to be withdrawn without special option
            $lane->is_withdrawn($new_status);
            $changed = 1;
        }
        elsif ($current_status == 1 && $new_status == 0) {
            issue_warning("Lane $lane_name was withdrawn in the database, but is no longer withdrawn in the sequence index!");
            # we only allow lanes to be unwithdrawn if special option provided
            if ($do_major_updates) {
                $lane->is_withdrawn($new_status);
                issue_warning("Lane $lane_name was unwithdrawn");
                $changed = 1;
            }
        }
    }
    else {
        # new lane, set the status
        $lane->is_withdrawn($new_status);
        $changed = 1;
    }
    
    
    # paired
    undef $current_status;
    $current_status = $lane->is_paired();
    
    # we can have both a single-ended fastq, and a pair of paired
    # fastqs for a lane; we'll call the lane paired if it contains a
    # pair
    $new_status = 0;
    while (my ($file, $status) = each %{$paired_statuses{$lane_name}}) {
        $new_status++ if (! $status || $status ne '');
    }
    if ($new_status == 1) {
        # if sequence.index claims this is paired, but there's only 1 fastq file
        # we disregard
        $new_status = 0;
    }
    elsif ($new_status > 1) {
        $new_status = 1;
    }
    
    if (defined $current_status) {
        if ($current_status == 0 && $new_status == 1) {
            issue_warning("Lane $lane_name just became paired");
            # we allow lanes to be paired without special option
            $lane->is_paired($new_status);
            $changed = 1;
        }
        elsif ($current_status == 1 && $new_status == 0) {
            issue_warning("Lane $lane_name was paired in the database, but is no longer paired in the sequence index!");
            # we only allow lanes to be unpaired if special option provided
            if ($do_major_updates) {
                $lane->is_paired($new_status);
                issue_warning("Lane $lane_name was unpaired");
                $changed = 1;
            }
        }
    }
    else {
        # new lane, set the status
        $lane->is_paired($new_status);
        $changed = 1;
    }
    
    if ($changed) {
        $lane->update || die "update failed";
    }
}

# are there lanes in our database no longer present in the sequence.index?
my $all_lanes = $vrtrack->qc_filtered_lane_names();
foreach my $lane_name (@{$all_lanes}) {
    unless (exists $si_lanes{$lane_name}) {
        issue_warning("Lane $lane_name is in the database, but no longer present in sequence.index!");
        if ($do_major_updates) {
            my $lane = VRTrack::Lane->new_by_name($vrtrack, $lane_name);
            $lane->is_processed('deleted', 1);
            $lane->update || die "update failed";
            issue_warning("Lane $lane_name was scheduled for deletion");
        }
    }
}

# say how many lanes we added before exiting
if ($verbose) {
    print STDERR "Parsed $si_lines lines in the sequence.index, of which $ok_si_lines were considered\n";
}
print STDERR "Added $new_lanes_count new lanes\n";

exit;

sub add_file {
    my ($lane, $path, $md5, $pair_path, $raw_reads, $raw_bases) = @_;
    my $fastq_filename = basename($path);
    
    # when withdrawn, the filename changes to remove the filt and the md5
    # becomes invalid; just skip
    return if $md5 && $md5 eq '................................';
    
    # since dcc can change the path of a file without changing the file, we
    # need to check all files of the lane and match on md5
    my $file;
    foreach my $obj (@{$lane->files || []}) {
        my $this_md5 = $obj->md5;
        if ($this_md5 eq $md5 || $obj->hierarchy_name eq $fastq_filename) {
            $file = $obj;
            last;
        }
    }
    
    if ($file) {
        # check for inconsistencies
        my $lane_name = $lane->hierarchy_name();
        my $this_md5 = $file->md5;
        my $hname = $file->hierarchy_name();
        
        my $name = $file->name;
        if ($path ne $name) {
            issue_warning("Lane $lane_name has a file with md5 $this_md5, but has path $name in the database but $path in sequence.index!");
            if ($do_major_updates) {
                $file->name($path);
                $file->update || die "update failed";
                issue_warning("Lane ${lane_name}'s file with md5 $this_md5 had its path updated");
            }
        }
        
        if ($hname ne $fastq_filename) {
            issue_warning("Lane $lane_name has a file with md5 $this_md5, but it is called $hname in the database but $fastq_filename in sequence.index!");
            if ($do_major_updates) {
                $file->hierarchy_name($fastq_filename);
                $file->update || die "update failed";
                issue_warning("Lane ${lane_name}'s file with md5 $this_md5 had its hierarchy name updated");
            }
        }
        
        my $counts_same = 1;
        my $these_raw_reads = $file->raw_reads() || 'not available';
        $raw_reads ||= 'not available';
        if ("$these_raw_reads" ne "$raw_reads" && "$these_raw_reads" ne 'not available' && $raw_reads ne 'not available') {
            $counts_same = 0;
            issue_warning("Lane $lane_name has a file with md5 $this_md5, but it has $these_raw_reads reads in the database but $these_raw_reads in sequence.index!");
            if ($do_major_updates) {
                $file->raw_reads($raw_reads);
                $file->update || die "update failed";
                my $existing = $lane->raw_reads;
                $existing ||= 0;
                $lane->raw_reads($existing - $these_raw_reads + $raw_reads);
                $lane->update || die "update failed";
                issue_warning("Lane ${lane_name}'s file with md5 $this_md5 had its reads updated");
            }
        }
        
        my $these_raw_bases = $file->raw_bases() || 'not available';
        $raw_bases ||= 'not available';
        if ("$these_raw_bases" ne "$raw_bases" && "$these_raw_bases" ne 'not available' && $raw_bases ne 'not available') {
            $counts_same = 0;
            issue_warning("Lane $lane_name has a file with md5 $this_md5, but it has $these_raw_bases bases in the database but $raw_bases in sequence.index!");
            if ($do_major_updates) {
                $file->raw_bases($raw_bases);
                $file->update || die "update failed";
                my $existing = $lane->raw_bases;
                $existing ||= 0;
                $lane->raw_bases($existing - $these_raw_bases + $raw_bases);
                $lane->update || die "update failed";
                issue_warning("Lane ${lane_name}'s file with md5 $this_md5 had its bases updated");
            }
        }
        
        # if the filename is the same and md5 has changed, we'll need to
        # reimport (and remap) it
        if (! $ignore_diff_md5s && $hname eq $fastq_filename && $md5 ne $this_md5) {
            if ($file->is_processed('import')) {
                issue_warning("Lane ${lane_name}'s file $hname was updated - the lane may need to be reimported and remapped!");
                if ($do_major_updates) {
                    unless ($only_update_diff_md5s && $counts_same) {
                        # schedule the lane to be cleaned and reimported by
                        # HierarchyUpdate pipeline
                        unless ($file->is_processed('altered_fastq')) {
                            $file->is_processed('altered_fastq', 1);
                            $file->update || die "update failed";
                        }
                        unless ($lane->is_processed('altered_fastq')) {
                            $lane->is_processed('altered_fastq', 1);
                            $lane->update || die "altered_fastq update failed";
                        }
                    }
                    $file->md5($md5);
                    $file->update || die "update failed";
                    unless ($only_update_diff_md5s && $counts_same) {
                        issue_warning("Lane $lane_name was scheduled for being reset to allow reimport and mapping");
                    }
                    else {
                        issue_warning("Lane $lane_name had the md5 for $fastq_filename updated, but it was NOT scheduled for reimport and remapping");
                    }
                }
            }
            else {
                # the file was probably entered into the db as withdrawn and got
                # it's md5 changed before being unwithdrawn and imported
                $file->md5($md5);
                $file->update || die "update failed";
            }
        }
    }
    else {
        $file = $lane->add_file($path);
        unless ($file) {
            my $lname = $lane->name;
            my $lid = $lane->id;
            die "lane $lname ($lid) had a new file '$fastq_filename' but it was already in the db?\n";
        }
        $file->hierarchy_name($fastq_filename);
        $file->md5($md5);
        
        if ($pair_path) {
            my $lane_name = $lane->name;
            my ($type) = $path =~ /${lane_name}_(\d)/;
            die "Couldn't determine forward/reverse for fastq [$path] (checked for ${lane_name}_(\\d))\n" unless $type;
            $file->type($type);
        }
        else {
            $file->type(0);
        }
        
        if ($raw_reads =~ /^\d+$/) {
            $file->raw_reads($raw_reads);
            my $existing = $lane->raw_reads;
            $existing ||= 0;
            $lane->raw_reads($existing + $raw_reads);
        }
        if ($raw_bases =~ /^\d+$/) {
            $file->raw_bases($raw_bases);
            my $existing = $lane->raw_bases;
            $existing ||= 0;
            $lane->raw_bases($existing + $raw_bases);
        }
        
        # mean_q ?
        # read_len ?
        # a separate script that will actually download the fastqs will have to
        # update these details
        
        $file->update;
        
        # since we've just added a file, we'll need to import it, which means
        # we'll need to change the lane's import status
        if ($lane->is_processed('import')) {
            $lane->is_processed('import', 0);
            $lane->is_processed('mapped', 0);
            $lane->update || die "update failed";
        }
    }
}

# set lane->is_processed(swapped => 1) on all decsendant lanes
sub swap_descendants {
    my $obj = shift;
    
    my @descendants = @{$obj->descendants};
    
    foreach my $child (@descendants) {
        if ($child->isa('VRTrack::Lane')) {
            $child->is_processed(swapped => 1);
            $child->update;
            my $swapped = $child->is_processed('swapped');
            unless ($swapped) {
                my $name = $child->name;
                die "Failed to set swapped on lane $name\n";
            }
        }
    }
}

sub issue_warning {
    my $warning = shift;
    return if exists $issued_warnings{$warning};
    print STDERR "\b" if $spinner;
    print $warning."\n";
    $issued_warnings{$warning} = 1;
    print STDERR " " if $spinner;
}

sub main_loop {
    my $rh = shift;
    
    # provide some feedback since this is really slow
    if ($spinner) {
        $spin++;
        if ($spin == 1) {
            print STDERR "\b-" if $spinner;
        }
        elsif ($spin == 2) {
            print STDERR "\b\\";
        }
        elsif ($spin == 3) {
            print STDERR "\b|";
        }
        elsif ($spin == 4) {
            print STDERR "\b/";
        }
        elsif ($spin == 5) {
            print STDERR "\b-";
            $spin = 0;
        }
    }
    
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
    #[25] ANALYSIS_GROUP
    
    $rh->[0] || die "got a sequence.index entry with no fastq file!\n";
    
    $rh->[2] || die "fastq $rh->[0] had no run id!\n";
    
    if ($debug) {
        return unless $rh->[2] eq $debug;
        warn "Found a $debug line...\n";
    }
    
    my $lane = VRTrack::Lane->new_by_name($vrtrack, $rh->[2]);
    if ($debug) {
        if ($lane) {
            warn "already in db, with id ", $lane->id, "\n";
        }
        else {
            warn "not in db\n";
        }
    }
    $si_lanes{$rh->[2]} = 1;
    # don't try and represent brand new lanes that are already withdrawn:
    # they may have been withdrawn due to some kind of swap, in which case
    # we can't correctly represent their incorrect heierarchy in the db,
    # since the db can only contain a single (the correct) heierarchy
    $si_lines++;
    if ($rh->[20] && ! $lane) {
        issue_warning("$rh->[2] is brand new and withdrawn, skipping") if $verbose;
        return;
    }
    $ok_si_lines++;
    
    # we do per-file transactions, since we can die at various points within
    # this loop
    $vrtrack->transaction_start();
    
    warn "started transaction\n" if $debug;
    
    # create project if it doesn't already exist
    $rh->[9] || die "fastq $rh->[0] had no sample name!\n";
    my $sample_lookup = uc($rh->[9]);
    my $population_name = exists $sample_data{$sample_lookup} ? $sample_data{$sample_lookup}->{pop} : '';
    unless ($population_name) {
        # not in sample file, see if it's in the database thanks to
        # load_vrtrack_allocations.pl
        my $ind = VRTrack::Individual->new_by_name($vrtrack, $sample_lookup);
        if ($ind) {
            my $pop = $ind->population();
            if ($pop) {
                $population_name = $pop->name;
            }
        }
    }
    die "Couldn't map population name from sample '$rh->[9]'\n" unless $population_name;
    # sanity-check that it matches population in the sequence.index, which in
    # new 2010 sequence.index files should match
    unless ($population_name eq $rh->[10]) {
        warn "population name in the sample file ($population_name) does not match the name in the sequence.index file ($rh->[10]) for $rh->[2], skipping\n";
        $vrtrack->transaction_rollback;
        return;
    }
    $rh->[3] || die "fastq $rh->[0] had no study id!\n";
    
    # on disc, instead of naming a project after the actual project code, we
    # now name it after the population and analysis group, so that the same
    # sample in different projects that were made the same way are merged
    # together come release time. In the db we have a seperate project for each
    # population-analysisgroup-study so that the correct study info is not lost
    my $ag = $rh->[25] || 'generic';
    my $psuedo_project_name = $rh->[3];
    my $psuedo_project_hname = $psuedo_project_name;
    if ($merge_projects) {
        $psuedo_project_hname = $population_name.'_'.$ag;
        $psuedo_project_hname =~ s/\s/_/g;
        $psuedo_project_name = $psuedo_project_hname.'_'.$rh->[3];
    }
    warn "psuedo project name is $psuedo_project_name\n" if $debug;
    my $project = $object_cache{project}->{$psuedo_project_name} || VRTrack::Project->new_by_name($vrtrack, $psuedo_project_name);
    unless ($project) {
        $project = VRTrack::Project->create($vrtrack, $psuedo_project_name);
        $project->hierarchy_name($psuedo_project_hname);
        $project->update;
    }
    unless ($project->hierarchy_name eq $psuedo_project_hname) {
        $project->hierarchy_name($psuedo_project_hname);
        $project->update || die "failed to update project hierarchy name to new population_analysis-group format\n";
    }
    $object_cache{project}->{$psuedo_project_name} ||= $project;
    warn "got project with id ", $project->id, "\n" if $debug;
    
    # create study as necessary, which ties into project
    my $study = $object_cache{study}->{$rh->[3]} || $project->study();
    unless ($study) {
        $study = $project->study($rh->[3]);
        $study = $project->add_study($rh->[3]) unless $study;
        
        unless ($study->acc eq $rh->[3]) {
            $study->acc($rh->[3]);
            $study->update || die "update failed";
        }
        
        $project->update;
    }
    $object_cache{study}->{$rh->[3]} ||= $study;
    unless ($project->study_id) {
        $project->study_id($study->id);
        $project->update;
    }
    warn "got study with id ", $study->id, "\n" if $debug;
    
    # create sample if it doesn't already exist (it uses the same identifiers
    # as will the individual object)
    my $sample = $object_cache{sample}->{"$rh->[3].$rh->[9]"} || $project->get_sample_by_name($rh->[9]);
    unless ($sample) {
        $sample = $project->add_sample($rh->[9]);
        if (exists $sample_data{$sample_lookup}) {
            unless ($sample_data{$sample_lookup}->{acc} eq $rh->[8]) {
                die "$sample_to_population_map_file and $si_file disagree on the accession for $rh->[9] ([".$sample_data{$sample_lookup}->{acc}."] vs [$$rh[8]])\n";
            }
        }
        $sample->update;
        
        $project->update;
    }
    $object_cache{sample}->{"$rh->[3].$rh->[9]"} ||= $sample;
    warn "got sample with id ", $sample->id, "\n" if $debug;
    
    # create individual, population and species as necessary (which tie into
    # sample)
    my $individual = $object_cache{individual}->{"$rh->[3].$rh->[9]"} || $sample->individual();
    unless ($individual) {
        $individual = $sample->individual($rh->[9]);
        $individual = $sample->add_individual($rh->[9]) unless $individual;
        if (exists $sample_data{$sample_lookup}) {
            $individual->alias($sample_data{$sample_lookup}->{alias});
            $individual->sex($sample_data{$sample_lookup}->{sex});
        }
        $individual->acc($rh->[8]);
        
        my $population = $individual->population;
        unless ($population) {
            $population = $individual->population($population_name);
            $individual->add_population($population_name) unless $population;
        }
        
        my $species = $individual->species();
        unless ($species) {
            die "species not set in db for $rh->[9]!\n" unless exists $sample_data{$sample_lookup};
            $species = $individual->species($sample_data{$sample_lookup}->{spp});
            $individual->add_species($sample_data{$sample_lookup}->{spp}) unless $species;
        }
        
        $individual->update;
        $sample->update;
    }
    $object_cache{individual}->{"$rh->[3].$rh->[9]"} ||= $individual;
    warn "got individual with id ", $individual->id, "\n" if $debug;
    
    # creating library (next step) won't be possible if we don't have a
    # sequencing centre
    unless ($rh->[5]) {
        #issue_warning("no seq centre for $rh->[2], rolling back and skipping!");
        #$vrtrack->transaction_rollback();
        #return;
        
        # This happens for ~100 lanes in the April release, and DCC basically
        # says they need to go through regardless, so just call them 'UNKNOWN'
        $rh->[5] = 'UNKNOWN';
    }
    
    # create library if it doesn't already exist.
    # since in rare cases a library can have been legitimately sequenced by
    # multiple centres, the library name in the db concatenates the name with
    # centre name, with hierarchy_name set to the original library name.
    # A library can also be used in multiple projects (eg. first used in pilot,
    # then used for topups of a individual for the main project), so we also
    # concatenate project id. Ideally the database would follow reality: one
    # library obj for one sample obj that belongs to many project objs, but
    # the api does not allow for the same sample_id to belong to multiple
    # project_ids. So now we have multple library objs associated with multiple
    # sample objs, each associated with their own project_id, even though we're
    # describing the same physical sample and library.
    $rh->[14] || die "fastq $rh->[0] had no library name!\n";
    my $old_library_name = $rh->[14].'|'.$rh->[5];
    my $library_name = $rh->[14].'|'.$rh->[5].'|'.$rh->[3];
    my $library = $object_cache{library}->{"$rh->[3].$rh->[9].$library_name"} || $sample->get_library_by_name($library_name);
    my $tech_name = $platform_to_tech{uc($rh->[12])} || die "Could not map platform '$rh->[12]' to a technology\n";
    unless ($library) {
        # conversion of old pre-project_id-concat to new form
        $library = $sample->get_library_by_name($old_library_name);
        if ($library) {
            $library->name($library_name);
            $library->update || die "Could not update library name";
        }
        else {
            $library = VRTrack::Library->new_by_name($vrtrack, $library_name);
            if ($library) {
                # shouldn't happen; it means a different sample has this library
                # - a swap?
                my $sid = $library->sample_id;
                unless ($sid) {
                    issue_warning("Library $library_name belongs to an unknown sample in the database, but $rh->[9] in sequence.index! (fastq $rh->[0])");
                }
                my $conflict_sample = VRTrack::Sample->new($vrtrack, $sid);
                my $expected_sid = $sample->id;
                my $conflict_sample_name = $conflict_sample->name;
                if ($sid != $expected_sid) {
                    issue_warning("A different sample has this library? Library $library_name belongs to sample $conflict_sample_name (db id $sid) in the database, but $rh->[9] (corresponding to db id $expected_sid) in sequence.index!") unless $rh->[20];
                    if ($do_major_updates && ! $rh->[20]) {
                        swap_descendants($library);
                        $library->sample_id($sample->id);
                        $library->update || die "update failed";
                        issue_warning("Library ${library_name}'s sample was updated");
                        #issue_warning("Library ${library_name}'s sample was NOT updated because not quite sure what to do in this situation yet");
                    }
                }
                else {
                    die "did sample($expected_sid)->get_library_by_name($library_name) and got no library, yet VRTrack::Library->new_by_name(vrtrack, $library_name) gave a library with sample_id $sid!\n";
                }
            }
            else {
                $library = $sample->add_library($library_name);
                my $hname = $rh->[14];
                $hname =~ s/[^\w\d]/_/g;
                $library->hierarchy_name($hname);
                
                if ($rh->[17]) {
                    # (unpaired lanes don't have insert sizes)
                    # allow for removal of insert size in schema v12 
                    if ($schema_version >= 12) {
						$library->fragment_size_from($rh->[17]);
						$library->fragment_size_to($rh->[17]);
                    } else {
	                    $library->insert_size($rh->[17]);
					}
                }
                
                # sort out the associated objects seq_centre and seq_tech
                my $seq_centre = $library->seq_centre();
                unless ($seq_centre) {
                    $seq_centre = $library->seq_centre($rh->[5]);
                    $library->add_seq_centre($rh->[5]) unless $seq_centre;
                }
                
                my $seq_tech = $library->seq_tech();
                unless ($seq_tech) {
                    $seq_tech = $library->seq_tech($tech_name);
                    $library->add_seq_tech($tech_name) unless $seq_tech;
                }
                
                $library->update;
            }
        }
    }
    $object_cache{library}->{"$rh->[3].$rh->[9].$library_name"} ||= $library;
    warn "got library with id ", $library->id, "\n" if $debug;
    
    # create lane if it doesn't already exist
    # we don't ask the library for the lane, since there may have been a swap,
    # the lane now no longer in the library as recorded in our database
    unless ($lane) {
        $lane = $library->add_lane($rh->[2]);
        $new_lanes_count++;
        
        issue_warning("just created new lane $rh->[2] and will add file $rh->[0]");
        
        add_file($lane, $rh->[0], $rh->[1], $rh->[19], $rh->[23], $rh->[24]);
        
        $lane->update;
        $library->update;
    }
    else {
        # if lane does exist, check for changes. Allow only withdrawn status to
        # change and for new fastqs to be added; everything else only report and
        # require a special command line option to actually make the change
        add_file($lane, $rh->[0], $rh->[1], $rh->[19], $rh->[23], $rh->[24]);
        $lane->update;
        
        # now check for major inconsistencies (swaps); only report on these
        # unless the special option has been supplied
        if (! $rh->[20] || $update_withdrawn) {
            warn "checking for major inconsistencies...\n" if $debug;
            
            my $current_lib = VRTrack::Library->new($vrtrack, $lane->library_id());
            my $current_lib_name = $current_lib->name;
            my $current_lib_id = $current_lib->id;
            my $correct_lib_id = $library->id;
            warn "current_lib_name is $current_lib_name ($current_lib_id vs $correct_lib_id)\n" if $debug;
            unless ($current_lib_name eq $library_name && $current_lib_id == $correct_lib_id) {
                if ($library_name =~ /UNKNOWN/ && $current_lib_name =~ /NCBI/) {
                    issue_warning("Lane $rh->[2] belongs to library $current_lib_name ($current_lib_id) in the database, so allowing conversion to $library_name ($correct_lib_id)");
                    $lane->is_processed(swapped => 1);
                    $lane->library_id($correct_lib_id);
                    $lane->update || die "Failed to update lane $rh->[2]\n";
                    issue_warning("Lane $rh->[2]'s library was updated");
                }
                else {
                    issue_warning("Lane $rh->[2] belongs to library $current_lib_name ($current_lib_id) in the database, but $library_name ($correct_lib_id) in sequence.index!");
                    if ($do_major_updates) {
                        $lane->is_processed(swapped => 1);
                        $lane->library_id($correct_lib_id);
                        $lane->update || die "Failed to update lane $rh->[2]\n";
                        issue_warning("Lane $rh->[2]'s library was updated");
                    }
                }
            }
            
            my $current_sample = VRTrack::Sample->new($vrtrack, $library->sample_id());
            my $current_sample_name = $current_sample->name;
            my $current_sample_id = $current_sample->id;
            my $correct_sample_id = $sample->id;
            warn "current_sample_name is $current_sample_name ($current_sample_id vs $correct_sample_id)\n" if $debug;
            unless ($current_sample_name eq $rh->[9] && $current_sample_id == $correct_sample_id) {
                issue_warning("Library $library_name belongs to sample $current_sample_name ($current_sample_id) in the database, but $rh->[9] ($correct_sample_id) in sequence.index!");
                if ($do_major_updates) {
                    swap_descendants($library);
                    $library->sample_id($correct_sample_id);
                    $library->update || die "update failed";
                    issue_warning("Library ${library_name}'s sample was updated");
                }
            }
            
            my $current_project = VRTrack::Project->new($vrtrack, $sample->project_id());
            my $current_project_name = $current_project->name;
            my $current_project_id = $current_project->id;
            my $correct_project_id = $project->id;
            warn "current_project_name is $current_project_name ($current_project_id vs $correct_project_id)\n" if $debug;
            unless ($current_project_name eq $psuedo_project_name && $current_project_id == $correct_project_id) {
                issue_warning("Sample $rh->[9] belongs to project $current_project_name ($current_project_id) in the database, but $psuedo_project_name ($correct_project_id) in sequence.index!");
                if ($do_major_updates) {
                    swap_descendants($sample);
                    $sample->project_id($project->id);
                    $sample->update || die "update failed";
                    issue_warning("Sample $rh->[9]'s project was updated");
                }
            }
        }
    }
    
    # we want to know afterwards if a lane became completely unwithdrawn, or
    # at least partially withdrawn
    $withdrawn_statuses{$rh->[2]}->{$rh->[0]} = $rh->[20];
    # likewise with paired status
    $paired_statuses{$rh->[2]}->{$rh->[0]} = $rh->[19];
    
    # check for other consistency problems with the sub-objects of the main
    # objects
    if (! $rh->[20] || $update_withdrawn) {
        warn "checking for sub-object inconsistencies...\n" if $debug;
        
        my $population = $individual->population;
        if ($population->name ne $population_name) {
            my $db_name = $population->name;
            issue_warning("The population name for individual $rh->[9] is supposed to be $population_name, not $db_name");
            if ($do_major_updates) {
                $population = $individual->population($population_name);
                $individual->add_population($population_name) unless $population;
                $individual->update;
                issue_warning("Population name for individual $rh->[9] corrected");
            }
        }
        
        my $species = $individual->species();
        if (exists $sample_data{$sample_lookup} && $species->name ne $sample_data{$sample_lookup}->{spp}) {
            issue_warning("Species is supposed to be $sample_data{$sample_lookup}->{spp}, but isn't for individual $rh->[9]!");
            if ($do_major_updates) {
                $species = $individual->species($sample_data{$sample_lookup}->{spp});
                $individual->add_species($sample_data{$sample_lookup}->{spp}) unless $species;
                $individual->update;
                issue_warning("Species name for individual $rh->[9] corrected");
            }
        }
        
        my $seq_centre = $library->seq_centre();
        if ($seq_centre->name ne $rh->[5]) {
            my $db_name = $seq_centre->name;
            issue_warning("The sequencing centre for library $library_name is supposed to be $rh->[5], not $db_name");
            if ($do_major_updates) {
                swap_descendants($library);
                $seq_centre = $library->seq_centre($rh->[5]);
                $library->add_seq_centre($rh->[5]) unless $seq_centre;
                $library->update;
                issue_warning("Sequencing centre for library $library_name corrected");
            }
        }
        
        my $seq_tech = $library->seq_tech();
        if ($seq_tech->name ne $tech_name) {
            my $db_name = $seq_tech->name;
            issue_warning("The sequencing technology for library $library_name is supposed to be $tech_name, not $db_name");
            if ($do_major_updates) {
                swap_descendants($library);
                $seq_tech = $library->seq_tech($tech_name);
                $library->add_seq_tech($tech_name) unless $seq_tech;
                $library->update;
                issue_warning("Sequencing technology for library $library_name corrected");
            }
        }
    }
    
    $vrtrack->transaction_commit();
}
