package G1K::Cleanup;

use Exporter;
@ISA = qw(Exporter);

@EXPORT = qw(cleanup);

use G1K::G1K;
use Utility;
use File::Basename;
use Cwd;

use strict; 

# g1k.pl cln
#
# Clean up e.g. v1 directory, deleting unneeded stuff in preparation
# for an rsync (backup).
#
# Can have 1 to 4 args after the "cln"; will recurse downwards through
# directory structure.

sub cleanup {
    my (@args) = @_;
    my (undef,$dir) = G1K::G1K::srcAndTgtDirs(@args);
    report("$dir: cleaning up");
    cleanupDir($dir);
    report();
}

sub cleanupDir {
    my ($d) = @_;
    my $b = basename($d);
    if ($b =~ /^NA.....$/) {
	cleanupIndividualDir($d);
    } elsif ($b =~ /^lib/) {
	cleanupLibraryDir($d);
    } elsif ($b =~ /^lane/) {
	cleanupLaneDir($d);
    } elsif ($b eq 'sample') {
	cleanupSampleDir($d);
    }
    foreach my $dd (glob("$d/*")) {
	cleanupDir($dd) if (-d $dd);
    }
}

# These might do something one day...

sub cleanupIndividualDir {
    my ($d) = @_;
}

sub cleanupLibraryDir {
    my ($d) = @_;
}

sub cleanupLaneDir {
    my ($d) = @_;
    # First, do the cleaning that we do anyway before a "g1k.pl map".
    G1K::RunMaq::maqMapClean($d);
    # Delete ref files if they're hard links, to save space in
    # the backup (symlinks should be OK). They're easily regenerated.
    foreach my $b (qw(ref.bfa ref.fa)) {
	my $f = "$d/$b";
	if (-f $f) {
	    unlink($f);
	}
    }
    # See if we have recalibrated.fastq.gz. If so, see if a "map maq"
    # run on this dir would use the recal file from the repository instead.
    # If so, we don't need this one. But if they're different sizes,
    # and the repository one is newer, warn and don't delete (something
    # weird is going on).
    my $rf = "$d/recalibrated.fastq.gz";
    if (-s $rf) {
	my $srcD = correspondingDataDir($d);
	# $crf should be the recal file in the repository if that's
	# present and up to date
	my $crf = G1K::RunMaq::chooseRecalibratedFastq($srcD,$d);
	if ($crf && ($crf ne $rf)) {
	    if (-s $crf != -s $rf) {
		report("WARNING: sizes differ: $rf $crf",1);
	    } else {
		`rm -f $d/recalibrated.fastq.gz*`;
	    }
	}
    }
}

sub correspondingDataDir {
    my ($d) = @_;
    my @comps = split(/\//,$d);
    $comps[-5] = 'DATA';
    return join("/",@comps);
}

# If qualmap file has been created, delete the sample dir altogether.

sub cleanupSampleDir {
    my ($d) = @_;
    my $dd = dirname($d);
    if (-s "$dd/qualmapBayesian.txt") {
	`rm -rf $d`;
    }
}
