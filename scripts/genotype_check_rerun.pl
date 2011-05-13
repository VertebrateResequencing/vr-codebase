#!/usr/bin/env perl

#
# genotype_check_rerun.pl --root /path/to/TRACKING --db dbname <file-of-lane-names>
#
# For example:
# genotype_check_rerun.pl --root /lustre/scratch106/projects/kuusamo/TRACKING 
# --db vrtrack_uk10k_cohort lanes.txt
#
# Sets up the lanes that need to have the genotype checking performed again. It deletes the gtype and 
# gtypex files and resets the processed flag on the lane table to 1.
#
# Author: John Maslen <jm23@sanger.ac.uk>
#

use strict;
use warnings;
no warnings 'uninitialized';

use Getopt::Long;
use VertRes::Utils::VRTrackFactory;
use VRTrack::Lane;
use File::Path;

my ($db, $root, $help);
my @filetypes = ('.gtype', '.gtypex');

GetOptions(
    'd|db=s'        =>  \$db,
    'r|root=s'		=>  \$root,
    'h|help'	    =>  \$help,
    );

($db && $root && !$help) or die <<USAGE;
    Usage: $0 <file of lanes that need genotype reanalysis performed>  
                --db        <specify db name>
                --root      <root directory for the analyses>
                --help      <this message>

Restarts the genotype checking for the lanes specified in the file provided.

USAGE
$db || die "a known --db must be specified\n";
if ('/' ne substr $root,-1,1) {
	$root = $root.'/';
}
my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $db,
                                                          mode => 'rw');
unless ($vrtrack){
    die "Can't connect to tracking database\n";
}

# get lanes from file of lane names or stdin
while (<>){
    my $lanename = $_;
    my @delfiles;
    chomp $lanename;
    #Check lane actually exists or bail out
    my $lane = VRTrack::Lane->new_by_name($vrtrack,$lanename);
    unless ($lane){
       print "Can't get lane $lanename\n";
       next;
    }
    #Delete files first
    #Get full path to qc-sample directory
    my $lanedir = $root.$vrtrack->hierarchy_path_of_lane_name($lane->name);
    my @mapping_ids = sort { $b<=>$a } @{$lane->mapping_ids()};
    my $qcdir = $lanedir.'/qc-sample.'.$mapping_ids[0].'/'; 
    #Get file names and delete if they exist
    for my $filetype ( @filetypes ) {
    	my $filepath = $qcdir.$lanename.$filetype;
    	if ( -e $filepath ) {
    		unlink $filepath || die "Unable to delete files for $lanename, error = $!\n"; 
    		print "Deleted: ", $filepath, "\n";
    	}
    	else {
       		print "$filepath does not appear to exist - check $lanename details.\n";
       	}
    }
    #update database
    my $sql = qq[update lane set processed = 1 where name = ?];
	my $sth = $vrtrack->{_dbh}->prepare($sql);
	if ($sth->execute(($lanename))) {
		print "Processed flag has been reset to 1 for lane ", $lanename, "\n";
	}
}