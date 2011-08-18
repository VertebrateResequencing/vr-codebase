#!/usr/bin/env perl
# 
# Author:       tk2
# Maintainer:   tk2
# Created:      Fri Sep 10 14:07:53 BST 2010 @588 /Internet Time/
# Updated:      Fri Sep 10 14:08:03 BST 2010 @588 /Internet Time/

use Carp;
use strict;
use warnings;

use File::Path;
use Cwd;

use VertRes::Utils::VRTrackFactory;
use VRTrack::VRTrack;
use VRTrack::Lane;

use Getopt::Long;
use Carp;

my ($database, $help, $lanes, $newRoot, $dry);

GetOptions
(
    'd|database=s'  =>\$database,
    'l|lanes=s'       =>  \$lanes,
    'dry|dry-run'     =>  \$dry,
    'r|root=s'        => \$newRoot,
    'h|help'	    =>  \$help,
);

( $database && -f $lanes && -d $newRoot && !$help) or die <<USAGE;
Constrcut a hierarchy from a list of lanes and a VRTrack db
maintaining the sym links to the lanes according to the storage path

    Usage: $0   
                --database <database name>
                --lanes	    <lanes to consider>
                --dry-run   <print out what changes will be made>
                --root      <new root of hierarchy>
                --help      <this message>
USAGE

my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $database, mode => 'r');

croak "Cant find database: $database" unless $vrtrack;

open( my $ifh, $lanes ) or die $!;
while( <$ifh> )
{
    chomp;
    
    my $lname = $_;
    my $lane = VRTrack::Lane->new_by_name($vrtrack, $lname);
    croak qq[cant find lane: $lane] unless $lane;
    
    my $lane_hier = $vrtrack->hierarchy_path_of_lane_name($lname);
    
    my $storagePath = $lane->storage_path();

    my $new_lane_path = qq[$newRoot/$lane_hier];
    if( -d $new_lane_path || -l $new_lane_path )
    {
        print qq[Already found lane path: $new_lane_path\n];
        next;
    }
    
    if( $storagePath )
    {
        #remove the lane from the path
        print qq[$lane_hier\n];
        my @t = split(/\//, $lane_hier );
        pop( @t );
        my $lib_path = join( '/', @t );
        
        my $path = qq[$newRoot/$lib_path];
        print qq[PATH: $path\n];
        if( ! -d $path )
        {
            if( $dry )
            {
                print qq[Make Dir: $path\n];
            }
            else
            {
                if( File::Path::mkpath( $path ) )
                {
                    print qq[Created lib path: $path\n];
                }else{croak qq[failed to make lib path: $path];}
            }
        }
        
        #now create the sym link to the storage path
        chdir( $path ) or croak qq[Failed to change to $path] unless $dry;
        
        if( $dry )
        {
            print qq[Making symlink: $storagePath -> $lname\n];
        }
        else
        {
            symlink( $storagePath, $lname ) or croak qq[Failed to symlink lane: $storagePath $lname];
        }
    }
    else
    {
        #make the hierarchy directory
        if( ! -d $new_lane_path )
        {
            if( $dry )
            {
                print qq[Making Dir: $new_lane_path\n];
            }
            else
            {
                File::Path::mkpath($new_lane_path) or croak qq[Failed to make directory: $new_lane_path\n];
            }
        }
    }
}
close( $ifh );

print qq[Verifying lanes in new root....\n];
#run a check on the new hierarchy to see if the lanes are all there
open( $ifh, $lanes ) or die $!;
while( <$ifh> )
{
    chomp;
    chomp;
    
    my $lname = $_;
    my $lane = VRTrack::Lane->new_by_name($vrtrack, $lname);
    croak qq[cant find lane: $lane] unless $lane;
    
    my $lane_hier = $vrtrack->hierarchy_path_of_lane_name($lname);
    
    my $path = qq[$newRoot/$lane_hier];
    if( ! -d $path && ! -l $path )
    {
        croak qq[Cant find lane $path\n];
    }
}
close( $ifh );

print qq[All done\n];
