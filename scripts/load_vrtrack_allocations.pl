#!/usr/bin/env perl
use strict;
use warnings;
no warnings 'uninitialized';

use Getopt::Long;
use VertRes::Utils::VRTrackFactory;
use VRTrack::Study;
use VRTrack::Individual;
use VRTrack::Seq_centre;
use VRTrack::Allocations;

my ($allocfile, $db, $help);

GetOptions(
    'a|alloc=s'     =>  \$allocfile,
    'd|db=s'        =>  \$db,
    'h|help'	    =>  \$help,
    );


(-f $allocfile && $db && !$help) or die <<USAGE;
    Usage: $0   
                --alloc     <file of allocation info>
                --db        <database to update, e.g. g1k_track>
                --help      <this message>

Loads G1K vrtrack database with new individuals and allocations of those
individuals to sequencing centres within projects.

Allocation file is munged from Assya's spreadsheet, and has fields:

0       Population
1       Coriell plate
2       sample id
3       SRA accession
4       family
5       sex
6       relationship
7       rel data
8       rel type
9       cell line
10      karyotype 
11      full proj center
12      full proj platform
13      full proj coverage
14      full proj >= 4x
15      P1 alloc centre
16      P1 Platform
17      P1 coverage
18      seq >= 4X?
19      LC29
20      P2 alloc centre
21      P2 platform
22      P2 coverage
23      P3 center
...      ...

Of course, it will probably have a different format, so you'll need to fix the
code...

USAGE

print "Database: $db\n";
my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $db,
                                                          mode => 'rw');

unless ($vrtrack){
    die "Can't connect to tracking database\n";
}

my $vrallocs = VRTrack::Allocations->new($vrtrack) or die "Can't get VRTrack::Allocations";

my %study_for_pop = (   'pilot1'    => 'SRP000031',
                        'pilot2'    => 'SRP000032',
                        'pilot3'    => 'SRP000033',
                        'CEU'	=> 'SRP000547',
                        'YRI'	=> 'SRP000542',
                        'CHB'	=> 'SRP000546',
                        'JPT'	=> 'SRP000544',
                        'TSI'	=> 'SRP000540',
                        'GBR'	=> 'SRP001294',
                        'FIN'	=> 'SRP000808',
                        'IBS'	=> 'SRP001514',
                        'CHS'	=> 'SRP001293',
                        'CDX'	=> 'SRP001515',
                        'KHV'	=> 'SRP001516',
                        'CHD'	=> 'SRP001517',
                        'LWK'	=> 'SRP000543',
                        'GWD'	=> 'SRP001518',
                        'GHN'	=> 'SRP001519',
                        'MAB'	=> 'SRP001520',
                        'ASW'	=> 'SRP000805',
                        'AJM'	=> 'SRP001521',
                        'ACB'	=> 'SRP001522',
                        'MXL'	=> 'SRP000803',
                        'CLM'	=> 'SRP001523',
                        'PEL'	=> 'SRP001524',
                        'PUR'	=> 'SRP001525',
                        'GIH'	=> 'SRP000806',
                        'MKK'	=> 'SRP000807',
                    );

my %code_for_centre = (
                        'ab'	    => 'ABI',
                        'bgi'	    => 'BGI',
                        'baylor'    => 'BCM',
                        'broad'	    => 'BI',
                        'illum'	    => 'ILLUMINA',
                        'illumina'  => 'ILLUMINA',
                        'maxplanck' => 'MPIMG',
                        'roche'	    => '454MSC',
                        'sanger'    => 'SC',
                        'washu'	    => 'WUGSC',
                     );

open (my $ALLOC, $allocfile) or die "Can't open $allocfile:$!\n";
while (<$ALLOC>){
    chomp;

    my ($pop, $samp, $acc, $sex);
    my @fields = map {$_ =~ s/\s+$//;$_} split "\t",$_;
    if ($fields[0] =~ /\((\w+)\)/){  # capture pop code from parens
        $pop = $1;
    }
    else {
        print "No population found in ",$fields[0],". skipping\n";
        next;
    }
    $samp = $fields[2];
    $acc = $fields[3];
    $sex = $fields[5];
    if ($sex eq 'male'){
        $sex = "M";
    }
    elsif ($sex eq 'female'){
        $sex = "F";
    }
    else {
        print "sex $sex not recognised.  Skipping\n";
        next;
    }
    # strip out whitespace and any trailing coverage number
    # e.g. go from "Illum 2X" to 'illum'
    foreach my $num(11,15,20,23){
        $fields[$num] =~ s/ .*//;
        $fields[$num] = lc($fields[$num]);
    }
    my %centre_for_project;
    $centre_for_project{main} = $fields[11] if $fields[11];
    $centre_for_project{pilot1} = $fields[15] if $fields[15];
    $centre_for_project{pilot2} = $fields[20] if $fields[20];
    $centre_for_project{pilot3} = $fields[23] if $fields[23];

    # Has there been any allocation at all?
    # If so, get the vr individual
    # If not, go to next record
    my $vind;
    if (scalar keys %centre_for_project){
        $vind = get_vr_individual($samp,$sex,'Homo sapiens',$pop, $acc);
    }
    else {
        print "No allocations for $pop $samp.  Skipping\n";
        next;
    }

    # for each of the projects with an allocation:
    # - get study code
    # - get study and individual from db, or create if not there
    # - add allocation
    # Note that this is not a switching statement - there can be all of the
    # projects specified on a single line
    foreach my $proj(keys %centre_for_project){
        my $centrename = $centre_for_project{$proj};
        my $centre = $code_for_centre{$centrename};
        unless ($centre){
            print "Centre $centrename not recognised.  Skipping\n";
            next;
        }
        my $study;
        if ($proj eq 'main'){
            $study = $study_for_pop{$pop};
        }
        else {
            $study = $study_for_pop{$proj}; # pilots have one study each
        }

        unless ($study){
            print "No study for $proj $pop.  Skipping\n";
            next;
        }
        my $vstudy = get_vr_study($study);
        my $vcentre = get_vr_centre($centre);
        add_allocation($vstudy, $vind, $vcentre);
        #print "$samp $proj $study $centre\n";

    }

}

close $ALLOC;

# retrieve a vrtrack study, or create one
sub get_vr_study {
    my $studyacc = shift;

    my $vstudy = VRTrack::Study->new_by_acc($vrtrack,$studyacc);
    unless ($vstudy){
        print "Adding study: $studyacc.\n";
        $vstudy = VRTrack::Study->create($vrtrack, $studyacc);
    }

    unless ($vstudy){
        die "ERROR can't create study : $studyacc\n";
    }

    return $vstudy;
}


# retrieve a vrtrack individual, or create one
sub get_vr_individual {
    my ($name,$sex,$spp,$pop,$acc) = @_;
    my $vind = VRTrack::Individual->new_by_name($vrtrack, $name);
    if ($vind){
        # check details
        my $vacc = $vind->acc;
        if ($vacc){
            unless ($vacc eq $acc){
                print "MISMATCH: $name $acc does not match $vacc in DB!\n";
            }
        }
        else {
            $vind->acc($acc);
        }
        my $vsex = $vind->sex;
        if ($vsex){
            unless ($vsex eq $sex){
                print "MISMATCH: $name $sex does not match $vsex in DB!\n";
            }
        }
        else {
            $vind->sex($sex);
        }
        $vind->update;
    }
    else {  # create new individual
        print "Adding individual: $name.\n";
        # need to set sex, alias, population, species.
        $vind = VRTrack::Individual->create($vrtrack, $name);
        $vind->sex($sex);
        $vind->alias($name);
        $vind->acc($acc) if $acc;
        my $vpop = $vind->population($pop);
        unless($vpop){
            print "\tNew population $pop\n";
            $vpop = $vind->add_population($pop);
        }
        my $vspp = $vind->species($spp);
        unless($vspp){
            print "\tNew spp ",$spp,"\n";
            $vspp = $vind->add_species($spp);
        }
        $vind->update;
    }

    unless ($vind){
        die "ERROR can't create individual : $name\n";
    }
    return $vind;
}


# retrieve a vrtrack centre, or create one
sub get_vr_centre {
    my $centrename = shift;
    my $vcentre = VRTrack::Seq_centre->new_by_name($vrtrack, $centrename);
    unless ($vcentre){
        print "Adding centre: $centrename.\n";
        $vcentre = VRTrack::Seq_centre->create($vrtrack, $centrename);
    }
    unless ($vcentre){
        die "ERROR can't create centre : $centrename\n";
    }
    return $vcentre;
}


sub add_allocation {
    my ($study, $ind, $centre) = @_;
    my $study_id = $study->id;
    my $centre_id = $centre->id;
    my $ind_id = $ind->id;

    if ($vrallocs->is_allocation_in_database($study_id, $ind_id, $centre_id)){
        print "Allocation ",join ", ",($study->acc,$ind->name,$centre->name)," already in db\n";
        # already in db
    }
    else {
        print "Adding allocation ",join ", ",($study->acc,$ind->name,$centre->name);
        print "\n";
        unless ($vrallocs->add_allocation($study_id, $ind_id, $centre_id)){
            die "Can't add allocation for $study_id, $ind_id, $centre_id\n";
        }
    }
}


