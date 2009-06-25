package G1K::Calibrate;

use Exporter;
@ISA = qw(Exporter);

@EXPORT = qw(sampleFastq calibrate);

use G1K::G1K;
use Utility;
use Bsub;
use File::Basename;
use Cwd;

use strict; 

my $BIN = $ENV{G1K} . "/bin";
my $MAP2QMAP = "$BIN/map2qMapFile.pl";
my $SNP126 = "$ENV{G1K}/misc/dbSNP/snp129.snp";

# g1k.pl sam command
sub sampleFastq {
    my ($anaType,$proj,$indiv,$lib,$lane) = @_;
    # Find the source (DATA/...) and target (v1/...) directories.
    my ($srcD,$tgtD) = srcAndTgtDirs(@_);
    # Bail out if this is a bad lane.
    return unless (G1K::RunMaq::inLaneList($tgtD));
    # See if we've already got a quality map file. If so, we don't
    # do sampling, even if there's no sample dir. This allows us
    # to delete sample dirs when we're done.
    return if (-s "$tgtD/qualmapBayesian.txt");

    # Check we can locate the fastq file(s) to recalibrate from.
    my ($src1, $src2) = laneFastqLink($srcD);
    return unless ($src);
    # We'll try and make $tgt from $src.
    my $tgt = "$tgtD/sample/30m.fastq";
    my $ms = makeStatus([$src],[$tgt]);
    return unless ($ms > 0 && stamp($tgt));
    my $fh = openToRead($src);
    # nrSum will be undef if no fastqcheck; so count how long
    # the readlength is manually
    my ($nrSum,$rl) = fastqReadCount($src);
    #warn "NRSUM: $nrSum,$rl";
    my $nrWanted; # number of reads we want
    if ($nrSum) {
	# $nrWanted is 30M/readLength, but if that's more than nrSum,
	# we can only have nrSum. NOTE: 30M/readLength is chosen so that
	# we have an upper limit on the time taken by g1k.pl mca,
	# regardless of read length. But it might be better to have a
	# fixed number of reads, as that would tend to equalise the
	# count per cell used for recalibration.
	$nrWanted = min($nrSum,30000000/$rl);
	my $pc = $nrWanted*100/$nrSum; # percent we'll sample:
	report(sprintf("%s: creating with %d%%",$tgt,$pc));
    } else {
	# we'll take the first nrWanted = 30M/readLength
	$nrSum = $nrWanted = 30000000/$rl;
	report(sprintf("%s: creating with first %d good reads",$tgt,$nrSum));
    }
    my $gh = openToWrite($tgt);
    my $nrOutput = 0;
    my $nrOrigWanted = $nrWanted;
    while ($nrSum >= 0) {
	# At any point in this loop, $nrSum is the number of reads left
	# in the file, and nrWanted is the number we still want, so
	# this definition of $want samples uniformly. If $nrSum was
	# initially undefined, $nrSum <= $nrWanted throughout and
	# $want will always be 1.
	my $want = (rand()*$nrSum <= $nrWanted);
	if ($want) {
	    # Read the fastq entry (4 lines)
	    my @quad;
	    foreach my $i (0 .. 3) {
		# 2008-09-27 jws strip down header lines so maq doesn't choke:
		my $line = <$fh>;
		if ($i == 0 || $i == 2){
		    $line =~ s/ .*//;	# delete anything after first space
		}
		$quad[$i] = $line;
	    }
	    # Read is good if the first 36 bases contains neither N
	    # nor a poly-A. But for paired-end we should probably use
	    # both halves. Occasionally this leads to problems when
	    # the first run was good but the second bad.
	    my $sq = substr($quad[1],0,36);
	    unless ($sq =~ /N/ || $sq =~ /AAAAAAAAAAAAAAA/) {
		print $gh @quad;
		$nrWanted --;
		$nrOutput ++;
	    }
	} else {
	    # Read and discard 4 lines.
	    foreach my $i (0 .. 3) {
		<$fh>;
	    }
	}
	$nrSum --;
	last if eof;
    }
    unstamp($tgt);
    report();
    report(sprintf("%d good reads wanted, %d found and output",$nrOrigWanted,$nrOutput),1);
}

# Return the number of reads in a fastq file or files and the read (pair) length.
sub fastqReadCount {
    my @files = @_;
    my ($sum,$rl);
    foreach my $f (@files) {
	$f = readlink($f) if (-l $f);
	my $fc; # find the fastqcheck file corresponding to a fastq
	foreach my $suf (qw(fastq fastq.gz fq.gz)) {
	    if ($f =~ /$suf$/) {
		($fc = $f) =~ s/$suf$/fastqcheck/;
		last if (-s $fc);
	    }
	}
	if (-f $fc) {
	    # If we got it, read off the number of reads directly:
	    # first field of first line
	    my $fh = openToRead($fc);
	    my $line = <$fh>;
	    my ($nr) = split(" ",$line);
	    $sum += $nr;
	}
	else {
	    # Otherwise count the number of lines directly
	    my $fh = openToRead($f);
	    while (<$fh>){
		$sum++;
	    }
	    $sum = int($sum/4);	# fastq, so reads = lines/4
	}
	# Get the read length from the first read in the # file itself
	my $fh = openToRead($f);
	<$fh>;  # skip header line
	chomp (my $line = <$fh>);
	$rl = length($line);
    }
    return ($sum,$rl);
}

# g1k.pl cal ...
# Use the output of g1k.pl mca (i.e. sample/aln.map) to derive a quality
# map file qualmapBayesian.txt, then apply that to the uncalibrated fastq
# file to derive recalibrated.fastq.gz.

sub calibrate {
    my ($anaType,$proj,$indiv,$lib,$lane) = @_;
    # Find source (DATA) and target (v1) dirs as usual, and return if this
    # is a bad lane.
    my ($srcD,$dir) = srcAndTgtDirs(@_);
    return unless (G1K::RunMaq::inLaneList($dir));
    # Input files needed
    my $src = "$dir/sample/aln.map";
    my $qmf = "$dir/qualmapBayesian.txt"; # or just qualmap.txt when fix done...
    my $ref = "$dir/ref.bfa";
    # We want to create this:
    my $tgt = "$dir/recalibrated.fastq.gz";
    # If we've already got the qualmap and there's a recalibrated fastq
    # file in the repository that's newer than the qualmap, assume we
    # don't need to do anything.
    if (-s $qmf && ! -s $tgt) {
	my $existing = G1K::RunMaq::chooseRecalibratedFastq($srcD,$dir);
	if (newer($existing,$qmf)) {
	    report("$dir: not creating recalibrated.fastq.gz because $existing is newer than qualmapBayesian.txt",1);
	    return;
	}
    }
    # Link to the right reference if required. It ought to be in the
    # sample directory.
    if (! -s $ref && -s "$dir/sample/ref.bfa") {
	link("$dir/sample/ref.bfa",$ref);
    }
    # If $qmf exists, we only require that; otherwise, we require $src and $ref.
    my $pReq = -s $qmf ? [$qmf] : [$src,$ref];
    my $ms = makeStatus($pReq,[$tgt]);
    if ($ms > 0 && stamp($tgt)) {
	my @cmds;
	# If $qmf doesn't exist, we have to create it.
	unless (-s $qmf) {
	    # Check we mapped enough reads to make it possible to create a
	    # reasonably accurate map.
	    if (enoughReadsMapped($src,1)) {
		my (undef,$isSingle) = G1K::G1K::getReadLength($srcD);
		my $single = $isSingle ? '-s' : '';
		# "-p $SNP126" masks out dbSNP positions.
		# Create qualmapLH.txt, then fix the recalibration with makeQualitiesBayesian.pl.
		# We ought to integrate the fix into $MAP2QMAP really.
		my $qmlh = "$dir/qualmapLH.txt";
		my $mqb = "$ENV{G1K}/bin/makeQualitiesBayesian.pl";
		unless (-s $qmlh) {
		    push @cmds,"$MAP2QMAP $single -p $SNP126 $src $ref > $qmlh";
		}
		push @cmds,"$mqb $qmlh > $qmf";
	    }
	}
	# If we either plan to make $qmf or already have it, we'll then create $tgt.
	if (@cmds || -s $qmf) {
	    my $fqLink = laneFastqLink($srcD);
	    die("$srcD: cannot find .fastq, .fastq.gz or .fq.gz link") unless ($fqLink);
	    # Apply $qmf to the uncalibrated data to create calibrated (and then remove the touchfile)
	    my $aqm = "$ENV{G1K}/bin/applyQualityMap.pl";
	    push @cmds,"$aqm $fqLink $qmf $tgt";
	    # Submit the job...
	    my $id = bsub($tgt,join('; ',@cmds));
	    report("submitted job $id for $tgt",1);
	} else {
	    # If we don't do anything, e.g. because not enough reads mapped, just unstamp.
	    unstamp($tgt);
	}
    }
}

sub enoughReadsMapped {
    my ($map,$rep) = @_;
    # Look for mapped-read counts in the aln.stat file corresponding to $map.
    (my $stat = $map) =~ s/map$/stat/;
    if (-s $stat) {
	my $n = 0;
	my $fh = openToRead($stat);
	while (<$fh>) {
	    if (/mapped [SP]E reads: (\d+) /) {
		$n += $1;
	    }
	}
	# Assume we need at least 100,000 reads for recalibration to be accurate.
	if ($n >= 100000) {
	    return 1;
	} else {
	    report("$map: only $n reads mapped, too few to recalibrate",1) if ($rep);
	    return 0;
	}
    } else {
	report("$map: aln.stat file missing or empty, not recalibrating",1) if ($rep);
	return 0;
    }
}

# @args will be the arguments to g1k.pl after (and not including)
# the action (command). For example if @args is
# (v1 Trio-CEU NA12878 libSC-NA12878-1 lane123_4), then
#   $srcD = $ENV{G1K}/DATA/Trio-CEU/NA12878/libSC-NA12878-1/lane123_4
#   $tgtD = $ENV{G1K}/v1/Trio-CEU/NA12878/libSC-NA12878-1/lane123_4
#   $badD = $ENV{G1K}/v1_BAD/Trio-CEU/NA12878/libSC-NA12878-1/lane123_4
# The first is where we find the pointer into the repository; the
# second is the analysis directory; and the third is where we move
# the analysis directory to if this lane is deemed to be bad data.

# jws: 2008-12-05 this is a copy from G1K::G1K because recal data is in
# RECALIBRATION not in DATA.  The DCC are taking over recalibration, so we want
# to keep it separate for now.
sub srcAndTgtDirs {
    my @args = @_;
    
    my $srcD = join("/",$ENV{G1K},'RECALIBRATION',@args[1..5]);
    my $tgtD = join("/",$ENV{G1K},@args[0..5]);
    my $badD = join("/",$ENV{G1K},$args[0] . '_BAD',@args[1..5]);
    $srcD =~ s/\/+$//;
    $tgtD =~ s/\/+$//;
    return ($srcD,$tgtD,$badD);
}


# the fastq to be recalibrated is symlinked into the RECALIBRATION directory
sub laneFastqLink {
    my ($d) = @_;
    my @links;
    foreach my $suf (qw(fastq fastq.gz)) {
	my @files = grep{ -l $_} glob("$d/*.$suf"); # find symlinks
	push @links, @files;
    }
    if (@links){
	if (scalar @links > 2){
	    die "More than two fastq file links in $d";
	}
	else {
	    return @links;
	}
    }
    return undef;
}

1;
