#!/usr/bin/env perl
use strict;
use warnings;
no warnings 'uninitialized';

use Getopt::Long;
use VertRes::Utils::VRTrackFactory;
use VRTrack::Lane;
use VertRes::Utils::FileSystem;
use File::Path qw(remove_tree);
use Cwd 'abs_path';
use File::Find;
use Filesys::DiskUsage qw/du/;

# TODO: proper software engineering.

my ( $db, $root, $clean, $help, $verbose, $dry_run );

GetOptions(
    'd|db=s'    => \$db,
    'r|root=s'  => \$root,
    'y|dry_run' => \$dry_run,
    'v|verbose' => \$verbose,
    'h|help'    => \$help,
);

( $db && !$help ) or die <<USAGE;
Usage: $0 [options] <file of lanes to delete>  
  -d STR specify db name
  -r STR root directory for the analyses
  -y     dont delete anything
  -v     verbose output
  -h     this message
  
  This script takes in a list of lanes and deletes any mapping and SNP calling 
  where there are duplicates (same reference+version, same mapper+version).
  It keeps the latest SNP calling/mapping and double checks the lanes do contain duplicates.
USAGE

my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(
    database => $db,
    mode     => 'rw'
);
unless ($vrtrack) {
    die "Can't connect to tracking database\n";
}

my $total_mappings_deleted    = 0;
my $total_snp_calling_deleted = 0;
print join( "\t", ( "Progress:", 'Lane', "MappingsDel", "TotalMappingsDel", "snpCallingDel", "TotalsnpCallingDel", "\n" ) );

# get lanes from file of lane names or stdin
while (<>) {
    my $lanename = $_;
    chomp $lanename;

    next if ( length($lanename) < 4 );

    #Check lane actually exists or bail out
    my $vrlane = VRTrack::Lane->new_by_name( $vrtrack, $lanename );
    unless ($vrlane) {
        print "Can't get lane $lanename\n";
        next;
    }
    next unless ( -d $root );
    my $lane_suffix_dir = $vrtrack->hierarchy_path_of_lane_name( $vrlane->name );

    # If you dont check this exists then you end up deleting the root directory of the pipeline
    next if ( !defined($lane_suffix_dir) || ( length($lane_suffix_dir) < 10 ) );
    my $lanedir = $root . $lane_suffix_dir . '/';
    next unless ( -d $lanedir );
    next unless ( $vrlane->is_processed('mapped') );

    print "Working on lane:\t$lanename\n";

    # mapped => 4,
    # snp_called => 256,

    my @files_to_delete     = [];
    my $mappings_deleted    = 0;
    my $snp_calling_deleted = 0;

    # Get all mappings (excluding QC)
    my %ref_mapper_mapstat_ids;
    my @raw_mapstat_ids;
    my %mapstats_to_prefix;
    for my $mapping_obj ( @{ $vrlane->mappings() } ) {
        next if ( $mapping_obj->is_qc() == 1 );
        next unless ( defined( $mapping_obj->bases_mapped() ) );

        # add the mapstats ID to a hash were the key is the mapper ID (which has a unique mapper+version) and reference (unique version)
        push( @{ $ref_mapper_mapstat_ids{ $mapping_obj->mapper_id() . '_' . $mapping_obj->assembly_id() } }, $mapping_obj->id() );
        push( @raw_mapstat_ids,                                                                              $mapping_obj->id() );
        $mapstats_to_prefix{ $mapping_obj->id() } = $mapping_obj->prefix();
    }

    # Work out which mappings have SNP calling
    my %snp_called_map_stats_ids;
    if ( $vrlane->is_processed('snp_called') ) {
        for my $mapstats_id (@raw_mapstat_ids) {
            my $snp_done_file = $lanedir . ${mapstats_id} . '.pe.markdup.snp/' . '.snps_done';
            next unless ( -e $snp_done_file );
            $snp_called_map_stats_ids{$mapstats_id}++;
        }
    }

    # filter out mappings where there are no duplicates
    for my $mapper_ref ( keys %ref_mapper_mapstat_ids ) {
        if ( @{ $ref_mapper_mapstat_ids{$mapper_ref} } <= 1 ) {
            delete( $ref_mapper_mapstat_ids{$mapper_ref} );
        }
    }

    my @mapstats_to_keep;
    my @mapstats_to_delete;

    # for each mapper_ref get the list of mappings which have been snp called
    for my $mapper_ref ( keys %ref_mapper_mapstat_ids ) {
        my @snp_called_mapstats;
        for my $mapstats_id ( @{ $ref_mapper_mapstat_ids{$mapper_ref} } ) {
            if ( $snp_called_map_stats_ids{$mapstats_id} ) {
                push( @snp_called_mapstats, $mapstats_id );
            }
        }

        my @sorted_mapstats;
        if ( @snp_called_mapstats < 1 ) {

            # If none, then ignore the snp calling bit
            @sorted_mapstats = sort { $b <=> $a } @{ $ref_mapper_mapstat_ids{$mapper_ref} };
        }
        else {
            # if theres at least 1, then the highest mapstats ID  thats snp called is the one to keep.
            @sorted_mapstats = sort { $b <=> $a } @snp_called_mapstats;
        }

        my $mapstats_id_to_keep = shift(@sorted_mapstats);
        push( @mapstats_to_keep,   $mapstats_id_to_keep );
        push( @mapstats_to_delete, @sorted_mapstats );

        for my $m (@sorted_mapstats) {
            if ( @snp_called_mapstats < 1 ) {
                print "MappingOnly:\tKeep\t${mapstats_id_to_keep}\tDel\t$m\n";
            }
            else {
                print "SNPandMapping:\tKeep\t${mapstats_id_to_keep}\tDel\t$m\n";
            }
        }
    }

    # snp calling
    if ( $vrlane->is_processed('snp_called') ) {
        for my $mapstats_id (@mapstats_to_delete) {

            my $snp_dir = $lanedir . ${mapstats_id} . '.pe.markdup.snp/';

            for my $file ($snp_dir) {
                push( @files_to_delete, $file );
            }
            $snp_calling_deleted++;
        }
    }

    # mapping
    for my $mapstats_id (@mapstats_to_delete) {
        my $mapping_prefix = $mapstats_to_prefix{$mapstats_id};
        next unless ( defined($mapping_prefix) || ( length($mapping_prefix) < 1 ) );

        for my $file (
            "split_pe_${mapstats_id}", "statistics_${mapstats_id}.pe.raw.sorted.bam",
            "merge_pe_${mapstats_id}", "mark_duplicates_pe_${mapstats_id}",
          )
        {
            foreach my $suffix (qw(o e pl jids)) {
                push( @files_to_delete, $lanedir . $mapping_prefix . $file . '.' . $suffix );
            }
        }

        for my $file (
            "${mapstats_id}.pe.raw.sorted.bam",         "${mapstats_id}.pe.raw.sorted.bam.bai",
            "${mapstats_id}.pe.markdup.bam",            "${mapstats_id}.pe.markdup.bam.bc",
            "${mapstats_id}.pe.markdup.bam.bai",        "${mapstats_id}.pe.raw.sorted.bam_graphs",
            "${mapstats_id}.pe.raw.sorted.bam.cover",   "${mapstats_id}.pe.raw.sorted.bam.bc",
            "${mapstats_id}.pe.raw.sorted.bam.bas",     "${mapstats_id}.pe.raw.sorted.bam.flagstat",
            "${mapstats_id}.pe.raw.sorted.bam.tmp",     "${mapstats_id}.pe.raw.sorted.bam.bad",
            "${mapstats_id}.pe.raw.sorted.bam.checked", "${mapping_prefix}job_status",
            "${mapping_prefix}split.jids",              ".mapping_complete_pe_${mapstats_id}",
            ".split_complete_pe_${mapstats_id}",        ".${mapstats_id}.pe.raw.sorted.bam.checked",
            "split_pe_${mapstats_id}"
          )
        {
            push( @files_to_delete, $lanedir . $file );
        }
        $mappings_deleted++;
    }

    my $fsu = VertRes::Utils::FileSystem->new( reconnect_db => 1 );
    for my $file (@files_to_delete) {

        if ( -d $file ) {
            print "Delete directory:\t" . $file . "\n";
            remove_tree($file) unless ($dry_run);
        }
        elsif ( -l $file ) {
            print "Delete link:\t" . $file . "\n";
            unlink($file) unless ($dry_run);
        }
        elsif ( -e $file ) {
            print "Delete file:\t" . $file . "\n";
            unlink($file) unless ($dry_run);
        }

        if ( !$dry_run ) {
            $fsu->file_exists( $file, recurse => 1, wipe_out => 1 );
        }

    }

    for my $mapping_obj ( @{ $vrlane->mappings() } ) {
        next if ( $mapping_obj->is_qc() == 1 );
        next unless ( defined( $mapping_obj->bases_mapped() ) );

        #update pipeline database to remove mapstats_id row
        for my $mapstats_id (@mapstats_to_delete) {
            if ( $mapping_obj->id() == $mapstats_id ) {

                print "Deleting mapstats DB entry:\t$mapstats_id\n";
                next if ($dry_run);
                $vrtrack->transaction_start();
                $mapping_obj->is_latest(0);
                $mapping_obj->update();
                $vrtrack->transaction_commit();
            }
        }

    }

    $total_mappings_deleted    += $mappings_deleted;
    $total_snp_calling_deleted += $snp_calling_deleted;

    print join( "\t",
        ( "Progress:", $vrlane->name, $mappings_deleted, $total_mappings_deleted, $snp_calling_deleted, $total_snp_calling_deleted, "\n" ) );

}

