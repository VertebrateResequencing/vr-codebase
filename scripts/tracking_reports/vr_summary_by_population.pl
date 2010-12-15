#!/usr/bin/env perl
use strict;
use warnings;
no warnings 'uninitialized';
use Getopt::Long;
use VertRes::Utils::VRTrackFactory;
use Date::Manip;

my ($projfile, $spp, $all_samples, $beforedate, $help);

GetOptions(
    'p|projects=s'  =>  \$projfile,
    's|spp=s'       =>  \$spp,
    'a|all'         =>  \$all_samples,
    'd|date=s'      =>  \$beforedate,
    'h|help'	    =>  \$help,
    );

my %db_for_spp = ( 'mouse'  => 'mouse_reseq_track',
		    'g1k'   => 'g1k_track',
	        );

my $db = $db_for_spp{$spp};

$beforedate ||= "3000-01-01";  # far future before date...
my $before = ParseDate($beforedate);

($db && !$help) or die <<USAGE;
    Usage: $0   
                --spp       <species, i.e. g1k or mouse>
               [--projects  <project name or file of SequenceScape project names>]
               [--all       <G1K only. Report all individuals including non-Sanger allocations.  Defaults to Sanger-only.>]
                --help      <this message>

Summarises VRTrack projects at the population level.

USAGE

warn "Species: $spp\n";
warn "Database: $db\n";
my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $db,
                                                          mode => 'rw');

unless ($vrtrack){
    die "Can't connect to tracking database\n";
}

my @projects;
if (-s $projfile){
    open (my $PROJ, "$projfile") or die "Can't open $projfile: $!\n";
    while (<$PROJ>){
	chomp;
        my ($project) = split /\t/, $_;
        push @projects, $project;
    }
    close $PROJ;
}
elsif ($projfile){
    push @projects, $projfile;
}
else {
    foreach (@{$vrtrack->projects}){
        push @projects, $_->name;   # wasteful, yes.
    }
}

print join "\t",('project','pop','individ.','raw individ. > 4X','pass individ. > 4X','total lanes','new lanes','qc pend lanes','passed lanes','failed lanes','total Gbp','total Gbp needed','passed Gbp','passed Gbp needed','new Gbp','pending Gbp','failed Gbp','Gbp submitted','Gbp to submit');
print "\n";
foreach my $pname (@projects){
    my $project = $vrtrack->get_project_by_name($pname);
    unless ($project){
        warn "Can't get project $pname.  Skipping\n";
        next;
    }
    my %sumdata;
    foreach my $sample(@{$project->samples}){
        if ($spp eq 'g1k' &! $all_samples){
            next unless $sample->is_sanger_sample;
        }
        my $pop = $sample->individual->population->name;
        my $ind = $sample->individual->name;
        unless (exists $sumdata{$pop}{$ind}){
            $sumdata{$pop}{$ind} = {};
        }
        my $ind_data = $sumdata{$pop}{$ind};
        foreach my $library (@{$sample->libraries}){
            $ind_data->{libs}++;
            foreach my $lane (@{$library->lanes}){
                my $dateflag = Date_Cmp($lane->run_date,$before);
                if ($dateflag >= 1){ 
                    # $run_date is not earlier than $before
                    next;
                }
                $ind_data->{lanes}++;
                $ind_data->{qc}{$lane->qc_status}++;
                $ind_data->{qc_bp}{$lane->qc_status} += $lane->raw_bases;
                my $bp = $lane->raw_bases;
                $ind_data->{bp_tot} += $bp;
                if (defined $lane->submission_id){
                    $ind_data->{bp_sub} += $bp;
                }
                elsif ($lane->qc_status eq 'passed'){
                    $ind_data->{bp_to_sub} += $bp;
                }
            }
        }
    }
    # stats & output
    foreach my $pop (keys %sumdata){
        my $indcount = 0;
        my $ind_over_4gb = 0;
        my $ind_over_4gb_pass = 0;
        my $total_bp = 0;
        my $total_bp_no_qc = 0;
        my $total_bp_pending = 0;
        my $total_bp_passed = 0;
        my $total_bp_failed = 0;
        my $bp_to_seq = 0;
        my $bp_to_seq_pass = 0;
        my $subbed_bp = 0;
        my $bp_to_sub = 0;
        my $total_lanes = 0;
        my $no_qc_lanes = 0;
        my $pending_lanes = 0;
        my $passed_lanes = 0;
        my $failed_lanes = 0;
        foreach my $ind (keys %{$sumdata{$pop}}){
            $indcount++;
            my $ind_data = $sumdata{$pop}{$ind};
            if ($ind_data->{bp_tot} > 12e9){    # > 4x sequenced
                $ind_over_4gb++;
            }
            else {
                $bp_to_seq += 12e9-$ind_data->{bp_tot};
            }
            if ($ind_data->{qc_bp}{'passed'} > 12e9){    # > 4x sequenced
                $ind_over_4gb_pass++;
            }
            else {
                $bp_to_seq_pass += 12e9-$ind_data->{qc_bp}{'passed'};
            }
            $total_bp += $ind_data->{bp_tot};
            $total_bp_no_qc += $ind_data->{qc_bp}{'no_qc'};
            $total_bp_pending += $ind_data->{qc_bp}{'pending'};
            $total_bp_passed += $ind_data->{qc_bp}{'passed'};
            $total_bp_failed += $ind_data->{qc_bp}{'failed'};
            $subbed_bp += $ind_data->{bp_sub};
            $bp_to_sub += $ind_data->{bp_to_sub};
            $total_lanes += $ind_data->{lanes};
            $no_qc_lanes += $ind_data->{qc}{no_qc};
            $pending_lanes += $ind_data->{qc}{pending};
            $passed_lanes += $ind_data->{qc}{passed};
            $failed_lanes += $ind_data->{qc}{failed};
        }

        $total_bp = sprintf('%.2f',$total_bp/1e9);
        $total_bp_no_qc = sprintf('%.2f',$total_bp_no_qc/1e9);
        $total_bp_pending = sprintf('%.2f',$total_bp_pending/1e9);
        $total_bp_passed = sprintf('%.2f',$total_bp_passed/1e9);
        $total_bp_failed = sprintf('%.2f',$total_bp_failed/1e9);
        $bp_to_seq = sprintf('%.2f',$bp_to_seq/1e9);
        $bp_to_seq_pass = sprintf('%.2f',$bp_to_seq_pass/1e9);
        $subbed_bp = sprintf('%.2f',$subbed_bp/1e9);
        $bp_to_sub = sprintf('%.2f',$bp_to_sub/1e9);

        print join "\t",($pname,$pop,$indcount,$ind_over_4gb, $ind_over_4gb_pass,$total_lanes,$no_qc_lanes,$pending_lanes,$passed_lanes,$failed_lanes,$total_bp,$bp_to_seq,$total_bp_passed,$bp_to_seq_pass,$total_bp_no_qc,$total_bp_pending,$total_bp_failed,$subbed_bp,$bp_to_sub);
        print "\n";
    }
}

