#!/usr/local/bin/perl

use strict;
use warnings;
no warnings 'uninitialized';
use Getopt::Long;

$|++;

my ($srf_file, $help);

GetOptions(
    's|srf=s'	=>  \$srf_file,
    'h|help'    =>  \$help,
    );

(-s $srf_file && !$help) or die <<USAGE;
    Usage: $0   
                --srf  <file of srf filenames.  Files will be pulled out of fuse>
		--help <this message>

Script to automate aspera copy of a set of srfs out of fuse.

USAGE

open (my $SRF, $srf_file) or die "Can't open $srf_file: $!\n";
my $ASCP = '~jws/apps/aspera/connect/bin/ascp';
my $password="aGsu7aiT";

while (<$SRF>){
    chomp;
    my $file = "/fuse/mpsafs/all/$_";
    unless (-f $file){
	warn "No file at $file\n";
	next;
    }
    my $command = "|$ASCP -QT -m 50M -l 400M $file era-drop-2\@fasp.era.ebi.ac.uk:";
    open(P,$command) || confess("Failed to open ".$command);
    print P "$password\n";
    close P;
}

