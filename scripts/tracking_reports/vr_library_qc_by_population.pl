#!/usr/bin/env perl
use strict;
use warnings;
no warnings 'uninitialized';
use Getopt::Long;
use VertRes::Utils::VRTrackFactory;

my ($projfile, $spp, $all_samples, $help);

GetOptions(
    'p|projects=s'  =>  \$projfile,
    's|spp=s'       =>  \$spp,
    'a|all'         =>  \$all_samples,
    'h|help'	    =>  \$help,
    );

my %db_for_spp = ( 'mouse'  => 'mouse_reseq_track',
		    'g1k'   => 'g1k_track',
	        );

my $db = $db_for_spp{$spp};

($projfile && $db && !$help) or die <<USAGE;
    Usage: $0   
                --spp       <species, i.e. g1k or mouse>
               [--projects  <project name or file of SequenceScape project names>]
               [--all       <G1K only. Report all individuals including non-Sanger allocations.  Defaults to Sanger-only.>]
                --help      <this message>

Summarises VRTrack library qc_status for each population

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

print join "\t",('project','pop','total libs','no_qc','pending','passed','failed');
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
            if (@{$library->lanes}){
                $ind_data->{libs}++;
                $ind_data->{$library->qc_status}++;
            }
        }
    }
    # stats & output
    foreach my $pop (keys %sumdata){
        my $total_libs = 0;
        my $no_qc_libs = 0;
        my $pending_libs = 0;
        my $passed_libs = 0;
        my $failed_libs = 0;
        foreach my $ind (keys %{$sumdata{$pop}}){
            my $ind_data = $sumdata{$pop}{$ind};
            $total_libs += $ind_data->{libs};
            $no_qc_libs += $ind_data->{no_qc};
            $pending_libs += $ind_data->{pending};
            $passed_libs += $ind_data->{passed};
            $failed_libs += $ind_data->{failed};
        }

        print join "\t",($pname,$pop,$total_libs,$no_qc_libs,$pending_libs,$passed_libs,$failed_libs);
        print "\n";
    }
}

