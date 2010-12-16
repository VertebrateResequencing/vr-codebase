#!/usr/bin/env perl
use strict;
use warnings;
no warnings 'uninitialized';
# GET LOCAL (DEV) VRTRACK, NOT $G1K/MODULES
use lib ".";

use Getopt::Long;
use File::Find;
use VRTrack::VRTrack;
use VRTrack::Lane;

my ($help, $era_dir);

GetOptions(
    'h|help'	    =>  \$help,
    'a|era_dir=s'   =>  \$era_dir,
    );

(!$help && -d $era_dir) or die <<USAGE;
    Usage: $0   
                --help <this message>
                --era_dir <directory that all the ERA xmls are under>

Loads submission info into G1K tracking db.

era_dir is the directory that the ERA xml tarball has been unpacked into.  The
script 'finds' into this directory to locate the appropriate xmls for our
submissions.

e.g. to run:
cd \$G1K/SRA/ERA_accessions
rm -rf ERA00*

wget ftp://ftp.era-xml.ebi.ac.uk/meta/xml/xml.era.tar.gz
tar -zxvf xml.era.tar.gz

$0 --era_dir \$G1K/SRA/ERA_accessions

Note that this script only updates unaccessioned lanes that have a submission id.  It does not check existing accessions, or check for lanes that have no submission id.

USAGE

my $vrtrack = VRTrack::VRTrack->new({ host => $ENV{VRTRACK_HOST},
                                    port => $ENV{VRTRACK_PORT},
                                    user => $ENV{VRTRACK_RW_USER},
                                    password => $ENV{VRTRACK_PASSWORD},
                                    database => 'g1k_track',
                                   });

unless ($vrtrack){
    die "Can't connect to tracking database\n";
}

# get the submission aliases for lanes without accessions
my $sql = qq[select distinct sub.name from submission sub, latest_lane l where l.submission_id = sub.submission_id and l.acc is null];
my %submissions;
my $sth = $vrtrack->{_dbh}->prepare($sql);

if ($sth->execute()){
    foreach(@{$sth->fetchall_arrayref()}){
        $submissions{$_->[0]} = "";
    }
}
else{
    die(sprintf('Cannot retrieve submissions: %s', $DBI::errstr));
}

unless (keys %submissions){
    warn "All lanes are already accessioned\n";
    exit;
}


# Find all the submission xmls for these submission aliases.  Should be exactly 1:1.
# 1) Find all the submission.xmls
# 2) grep each for each of the submission aliases and die if find more than one match
open( my $FIND, "-|", "find", $era_dir, "-type", "f", "-name","*.submission.xml");
while (<$FIND>) {
    chomp;
    foreach my $sub_alias(keys %submissions){
        my $subxml = `grep -l 'alias="$sub_alias"' $_`;
        if ($subxml){
            if ($submissions{$sub_alias}){
                die "Found $subxml but already found ",$submissions{$sub_alias}," match for $sub_alias\n";
            }
            else {
                $submissions{$sub_alias} = $subxml;
            }
        }
    }
}

close ($FIND);

foreach my $sub_alias(keys %submissions){
    my $subxml = $submissions{$sub_alias};
    unless ($subxml){
        warn "submission alias $sub_alias not found in ERA files\n";
        next;
    }

    # get the submission accession for this alias
    # Hack the xml parsing with regex.
    open (my $SUBXML, $subxml) or die "Can't open $subxml: $!\n";
    my $sub_acc;
    while (<$SUBXML>){
        next unless /accession="([^"]+)"/i;
        $sub_acc = $1;
        last;
    }
    close ($SUBXML);
    unless ($sub_acc){
        warn "Can't get accession for $sub_alias in $subxml\n";
        # don't "next" as it's still useful to get any run accessions we can find
    }

    my $runxml = $subxml;
    $runxml =~ s/submission/run/;

    # open the run.xml file for this submission and pull out the runs, aliases, and accessions.
    open (my $RUNXML, $runxml) or die "Can't open $runxml: $!\n";
    while (<$RUNXML>){
        next unless /<RUN /;
        $_ =~ /alias="([^"]+)"/i;        
        my $lane_name = $1;
        $lane_name =~ s/_s_/_/; # just in case any submissions have the old _s_ names

        $_ =~ /accession="([^"]+)"/i;        
        my $lane_acc = $1;

        unless ($lane_name && $lane_acc){
            warn "Can't get name and accession for runline:\n$_\n";
            next;
        }

        my $lane = VRTrack::Lane->new_by_hierarchy_name($vrtrack,$lane_name);
        unless ($lane){
           print  "Can't find lane $lane_name in database\n";
           next;
        }
        if ($lane->submission){
            if ($lane->submission->name ne $sub_alias){
                # major problem somewhere
                die "Existing submission for $lane_name is ",$lane->submission->name, " not $sub_alias";
            }
        }
        else {
            warn "$lane_name has no associated submission.  Attempting to add.\n";
        }

        my $submission = $lane->submission($sub_alias);
        unless ($submission){
            $submission = $lane->add_submission($sub_alias);
            $sub_alias =~ /sc-(200\d{5})/;   # get date from name
            my $sub_date = $1;
            $submission->acc($sub_acc);
            $submission->date($sub_date) if $sub_date;
            $submission->update;
        }
        # update existing submission acc
        unless ($submission->acc){
            $submission->acc($sub_acc);
            $submission->update;
        }

        if ($lane->acc){
            if ($lane->acc ne $lane_acc){
                die "Lane has acc of ",$lane->acc," rather than $lane_acc\n";
            }
        }
        $lane->acc($lane_acc);

        $lane->update;
    }
    close ($RUNXML);

}



