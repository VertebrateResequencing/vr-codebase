#!/usr/bin/env perl
use strict;
use warnings;
no warnings 'uninitialized';
use Getopt::Long;
use VertRes::Utils::VRTrackFactory;
use VRTrack::Multiplex_pool;

my ($projfile, $db, $coverage, $help);

GetOptions(
    'p|projects=s'  =>  \$projfile,
    'd|db=s'        =>  \$db,
    'c|coverage=f'  =>  \$coverage,
    'h|help'	    =>  \$help,
    );

($db && !$help) or die <<USAGE;
    Usage: $0   
                --db        <database, e.g. vrtrack_uk10k_cohort>
               [--projects  <project name or file of SequenceScape project names>]
               [--coverage  <only report samples with less than this X coverage>]
                --help      <this message>

Lists each individual with its coverage, and any qc_passed libraries.

USAGE

warn "Database: $db\n";
my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $db,
                                                          mode => 'r');

unless ($vrtrack){
    die "Can't connect to tracking database\n";
}

# build list of project names to get
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

print join "\t",('project','pop','individ.','open reqs','pass reqs','fail reqs','total Gbp','GT conf Gbp','QC pass Gbp','map pass Gbp','pend Gbp','GT unchecked Gbp', 'GT unconf Gbp','GT unknown Gbp','GT wrong Gbp','pass cov','map pass cov','pend cov','pass libs','pend libs','fail libs');
print "\n";
foreach my $pname (@projects){
    my $project = $vrtrack->get_project_by_name($pname);
    unless ($project){
        warn "Can't get project $pname.  Skipping\n";
        next;
    }
    my %sumdata;
    my %pop_for_ind;
    foreach my $sample(@{$project->samples}){
        my $ind = $sample->individual->name;
        my $pop = $sample->individual->population->name;
        $pop_for_ind{$ind} = $pop;
        unless ($sumdata{$ind}){
            $sumdata{$ind}{total_bp} = 0;           # total sequenced basepairs from all lanes
            # sequenced basepairs by genotype status
            $sumdata{$ind}{gt_bp} = {   confirmed => 0,
                                        unchecked => 0,
                                        unconfirmed => 0,
                                        unknown => 0,
                                        wrong => 0,
                                    };             

            # sequencing requests
            # ignore 'cancelled' and 'hold'
            $sumdata{$ind}{open_reqs} = 0;      # requests started, pending, unknown
            $sumdata{$ind}{pass_reqs} = 0;    
            $sumdata{$ind}{fail_reqs} = 0;   

            $sumdata{$ind}{total_pass_bp} = 0;      # total basepairs from passed lanes
            $sumdata{$ind}{map_pass_bp} = 0;        # extrapolated from the mapping QC %
            $sumdata{$ind}{total_pend_bp} = 0;     # no_qc and pending lanes
            $sumdata{$ind}{pass_libs} = [];
            $sumdata{$ind}{fail_libs} = [];
            $sumdata{$ind}{pend_libs} = [];
        }

        foreach my $library (@{$sample->libraries}){
            # might need to revisit these later - if we have an open request that we actually
            # have data for, then it isn't open any more?
            my @seq_reqs;
            push @seq_reqs, @{$library->seq_requests};
            # and now for multiplex pools:
            # note that we'll count the seq_request more than once on a multiplex
            foreach my $mpp(@{$library->library_multiplex_pools}){
                my $mplex = VRTrack::Multiplex_pool->new($vrtrack, $mpp->multiplex_pool_id);
                push @seq_reqs, @{$mplex->seq_requests};
            }

            foreach my $sr (@seq_reqs){
                if ($sr->seq_status eq 'passed'){
                    $sumdata{$ind}{pass_reqs}++;
                }
                elsif ($sr->seq_status eq 'failed'){    
                    $sumdata{$ind}{fail_reqs}++;
                }
                elsif ( $sr->seq_status eq 'unknown' ||
                        $sr->seq_status eq 'pending'||
                        $sr->seq_status eq 'started'){
                    $sumdata{$ind}{open_reqs}++;
                }
            }


            # only show open and uncancelled passed & pending libs.
            # failed libs can be closed or cancelled, doesn't really matter
            if ($library->qc_status eq 'passed'){
                if ($library->open() && $library->prep_status ne 'cancelled'){
                    push @{$sumdata{$ind}{pass_libs}}, $library->name;
                }
            }
            elsif ($library->qc_status eq 'failed'){
                push @{$sumdata{$ind}{fail_libs}}, $library->name;
            }
            else {
                # this includes library not yet made, or in progress, as well as those
                # completed with a test lane and awaiting our qc.
                if ($library->open() && $library->prep_status ne 'cancelled'){
                    push @{$sumdata{$ind}{pend_libs}}, $library->name;
                }
            }

            foreach my $lane (@{$library->lanes}){
                # 2010-05-06 jws changed to include all genotype statuses.
                #next unless $lane->genotype_status eq 'confirmed';

                my $bp = $lane->raw_bases;
                $sumdata{$ind}{total_bp} += $bp;
                unless (defined $sumdata{$ind}{gt_bp}{$lane->genotype_status}){
                    warn "Unexpected lane genotype status : ",$lane->genotype_status,"\n";
                }
                else {
                    $sumdata{$ind}{gt_bp}{$lane->genotype_status} += $bp;
                }
                if ($lane->genotype_status eq 'confirmed'){
                
                    if ($lane->qc_status eq 'passed'){
                        $sumdata{$ind}{total_pass_bp} += $bp;
                        my $mapstats = $lane->latest_mapping;
                        if ($mapstats && $mapstats->raw_bases){
                           $sumdata{$ind}{map_pass_bp} += $bp * ($mapstats->rmdup_bases_mapped/$mapstats->raw_bases);
                        }
                        else {
                            $sumdata{$ind}{map_pass_bp} += $bp * 0.9;  # not sure what else to do here?
                        }
                    }
                    elsif ($lane->qc_status eq 'no_qc' or $lane->qc_status eq 'pending' or $lane->qc_status eq "investigate"){
                        # shouldn't be any no_qc lanes with genotype confirmed, as this
                        # is part of genotype checking, but leave it in anyway
                        $sumdata{$ind}{total_pend_bp} += $bp
                    }
                }
                
            }
        }

    }

    
    # output
    foreach my $ind (keys %sumdata){
        if ($coverage){
            next unless $sumdata{$ind}{map_pass_bp} < $coverage*3e9;
        }
        
        my $pass_reqs = $sumdata{$ind}{pass_reqs};
        my $fail_reqs = $sumdata{$ind}{fail_reqs};
        my $open_reqs = $sumdata{$ind}{open_reqs};

        my $total_bp = $sumdata{$ind}{total_bp};
        my $total_pass_bp = $sumdata{$ind}{total_pass_bp};
        my $map_pass_bp = $sumdata{$ind}{map_pass_bp};
        my $total_pend_bp = $sumdata{$ind}{total_pend_bp};
        my $pop = $pop_for_ind{$ind};
        my $pass_libs = $sumdata{$ind}{pass_libs};
        my $fail_libs = $sumdata{$ind}{fail_libs};
        my $pend_libs = $sumdata{$ind}{pend_libs};
        my $total_gtconf_bp = sprintf('%.2f',$sumdata{$ind}{gt_bp}{confirmed}/1e9);
        my $total_gtuncheck_bp = sprintf('%.2f',$sumdata{$ind}{gt_bp}{unchecked}/1e9);
        my $total_gtunconf_bp = sprintf('%.2f',$sumdata{$ind}{gt_bp}{unconfirmed}/1e9);
        my $total_gtunk_bp = sprintf('%.2f',$sumdata{$ind}{gt_bp}{unknown}/1e9);
        my $total_gtwrong_bp = sprintf('%.2f',$sumdata{$ind}{gt_bp}{wrong}/1e9);

        my $total_pass_cov = sprintf('%.1f',$total_pass_bp/3e9);
        my $map_pass_cov = sprintf('%.1f',$map_pass_bp/3e9);
        my $total_pend_cov = sprintf('%.1f',$total_pend_bp/3e9);

        $total_bp = sprintf('%.2f',$total_bp/1e9);
        $total_pass_bp = sprintf('%.2f',$total_pass_bp/1e9);
        $map_pass_bp = sprintf('%.2f',$map_pass_bp/1e9);
        $total_pend_bp = sprintf('%.2f',$total_pend_bp/1e9);

        print join "\t",($pname,$pop,$ind, $open_reqs, $pass_reqs,$fail_reqs, $total_bp,$total_gtconf_bp,$total_pass_bp,$map_pass_bp,$total_pend_bp,$total_gtuncheck_bp,$total_gtunconf_bp,$total_gtunk_bp,$total_gtwrong_bp,$total_pass_cov, $map_pass_cov, $total_pend_cov, (join ",",@$pass_libs),(join ",",@$pend_libs),(join ",",@$fail_libs));
        print "\n";
    }
}

