#!/usr/bin/env perl

#
# genotype_check_rerun.pl --root /path/to/TRACKING --db dbname <file-of-lane-names>
# (optional --verbose flag can be set to output actions of the script)
#
# For example:
# genotype_check_rerun.pl --root /lustre/scratch106/projects/kuusamo/TRACKING 
# --db vrtrack_uk10k_cohort --verbose lanes.txt
#
# Sets up the lanes that need to have the genotype checking performed again. It deletes the gtype and 
# gtypex files and resets the lane processed flag for qc to 0 ($lane->is_processed(qc => 0)).
#
# Author: John Maslen <jm23@sanger.ac.uk>
#

use strict;
use warnings;
no warnings 'uninitialized';

use Getopt::Long;
use VertRes::Utils::VRTrackFactory;
use VRTrack::Lane;

my ($db, $root, $help, $verbose);
my @filetypes = ('.gtype', '.gtypex','.gtypey','.glf','.bcf');

GetOptions(
    'd|db=s'        =>  \$db,
    'r|root=s'		=>  \$root,
    'v|verbose'     =>  \$verbose,
    'h|help'	    =>  \$help,
    );

($db && $root && !$help) or die <<USAGE;
    Usage: $0 <file of lanes that need genotype reanalysis performed>  
                --db        <specify db name>
                --root      <root directory for the analyses>
                --verbose   [be extra verbose about what it is doing]
                --help      <this message>

Restarts the genotype checking for the lanes specified in the file provided.

USAGE

if ('/' ne substr $root,-1,1) {
	$root = $root.'/';
}
die "Root directory $root does not exist\n" unless (-e $root);

my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $db,
                                                          mode => 'rw');
unless ($vrtrack){
    die "Can't connect to tracking database\n";
}

# get lanes from file of lane names or stdin
while (<>){
    my $lanename = $_;
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
    #if $qcdir does not exist, use /qc-sample/
    unless (-d $qcdir) {
    	$qcdir = $lanedir.'/qc-sample/'; 
    }
    
    #Get file names and delete if they exist
    for my $filetype ( @filetypes ) {
    	my $filepath = $qcdir.$lanename.$filetype;
    	warn "WARNING: $filepath does not exist\n" unless ( -e $filepath );
    	if ( -e $filepath ) {
    		unlink $filepath || die "Unable to delete $filepath, error = $!\n"; 
    		print "Deleted: ", $filepath, "\n" if ($verbose);
    	}
    }
    # also remove error and output files
    for my $outfile ('_glf.o','_glf.e'){
    	my $filepath = $qcdir.'_'.$lanename.$outfile;
    	if ( -e $filepath ) {
            unlink $filepath || die "Unable to delete $filepath, error = $!\n"; 
            print "Deleted: ", $filepath, "\n" if ($verbose);
    	}
    }

    #update database
    $lane->is_processed(qc => 0);
    $lane->update;
    print "Processed flag for qc has been reset to 0 for lane ", $lanename, "\n" if ($verbose);
}
