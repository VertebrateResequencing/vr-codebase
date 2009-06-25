package G1K::RunMaq;

use Exporter;
@ISA = qw(Exporter);

@EXPORT = qw(runMaqMap runMaqMapCal runMaqMerge runMaqRmdup runMaqPost runMaqValidate);

use G1K::G1K;
use Utility;
use Bsub;
use File::Basename;
use Cwd;

use strict; 

my $BIN = "$ENV{G1K}/bin";

my $MAQ = "$BIN/maq";
my $MAQ_NT = "$BIN/maq-nt";
my $MAQ_POST = "$BIN/maq_post.pl";
my $MAQ_SANGER = "$BIN/maq_sanger.pl";

# INPUTS:  $refGenome
#          $fqLink (will be readlinked to $fqFile)
# OUTPUTS: $tgtD/aln.map

# "map" action: run maq map on all data for a lane
sub runMaqMap {
    runMaqMap1('map',@_);
}

# "mca" action: run maq map on a sample, for calibration
sub runMaqMapCal {
    runMaqMap1('mca',@_);
}

sub runMaqMap1 {
    my $mapOrMca = shift @_;
    # Arguments to g1k.pl, from v1 onwards:
    my ($anaType,$proj,$indiv,$lib,$lane) = @_;
    # Not sure why this name has to be qualified; I'm exporting it...
    my $root = G1K::G1K::rootDir();
    my ($srcD,$tgtD,$badD) = G1K::G1K::srcAndTgtDirs(@_);
    # Check this lane is a valid one. If not, move tgtD to badD if
    # it exists, and return.
    if (!inLaneList($tgtD)) {
	report("$tgtD: in badLanes.txt, or not in goodLanes.txt",1);
	moveToBad($tgtD,$badD) if (-d $tgtD);
	return;
    }
    # Decide what file to map.
    my $fqFile;
    if ($mapOrMca eq 'map') {
	# "g1k.pl map" runs on recalibrated data. Fairly soon if not
	# now, this should be in the repository. But to allow
	# map to be run immediately after recalibrating, the best
	# thing is probably to set $fqFile as here if it exists,
	# and otherwise to look for *.recal.fastq.gz in the repos
	# (via $srcD).
	$fqFile = chooseRecalibratedFastq($srcD,$tgtD);
    } else {
	# "g1k.pl mca". Stuff happens in $tgtD/sample.
	$tgtD = "$tgtD/sample";
	$fqFile = "$tgtD/30m.fastq";
    }
    # (Gender-dependent) reference genome to use; will be linked
    # into place in $tgtD if does not already exist.
    my $refGenome = getRefGenome($root,$indiv,$tgtD);
    # File we're aiming to create.
    my $tgtF = "$tgtD/aln.map";
    # Do we need to create it, and can we? If so, we stamp it
    # to prevent another process doing so, and off we go.
    # The touch file ($tgtF.touch) will be deleted by the last
    # to finish of the jobs that maq_sanger.pl sets in motion.
    # Its continued existence is a sign that things have gone wrong.
    my $ms = makeStatus([$fqFile,$refGenome],[$tgtF],0,1);
    if ($ms > 0 && stamp($tgtF)) {
	# We're good to go. First make sure any old maq map debris
        # has been deleted.
	maqMapClean($tgtD);
        # Find the read length and whether
	# this lane is single or paired.
	my ($readLength,$isSingle) = G1K::G1K::getReadLength($srcD);
	# Dig out max insert length if defined.
	my $maxIL = maxInsertLength(dirname($tgtD));

	# 2008-10-05 jws - change default to 1000 if maxIL not set
	my $switches = defined $maxIL ? "-a $maxIL " : "-a 1000 ";
	# Tell maq the read length and whether single/paired.
	if ($isSingle) {
	    $switches .= " -s -3 $readLength";
	} else {
	    my $rl2 = 2*$readLength;
	    $switches .= " -l $readLength";
	}
	# Warn maq if fastq file is compressed.
	$switches .= ' -z' if ($fqFile =~ /gz$/);
	# We'll call maq_sanger.pl with the specified arguments. It's
	# a black box from the point of view of this code (though maybe
        # its contents should be integrated here. Perhaps talk to Heng).
	my $cmd = "$MAQ_SANGER $switches $tgtD $refGenome $fqFile";
	print "$cmd\n";
	# cmd0 will contain the maq_sanger.pl command itself, in
	# case we need it later for debugging.
	my $cmdFile0 = "$tgtD/cmd0";
	my $gh = openToWrite($cmdFile0);
	print $gh "$cmd\n";
	$gh->close();
	report("$tgtD: running maq_sanger.pl",1);
	# Run $cmd, which itself generates commands; save those
	# commandsin $tgtD/cmd; and run it, which should submit
	# some jobs to do the real work.
	my $cmdF = "$tgtD/cmd";
	`$cmd >$cmdF; sh $cmdF`;
    }
}

sub maqMapClean {
    my ($tgtD) = @_;
    return unless (-s "$tgtD/aln.map");
    my ($touch) = glob("$tgtD/*.touch");
    return if ($touch); # don't clean a directory that's being worked on
    my @torm = qw(aln-rmdup.log aln-rmdup.map *.bfq.* *.cat.* *.map.* *.post.* map read[12] ref.*fa run_script.sh unmap);
    my $cwd = cwd();
    chdir($tgtD);
    my $ae = "aln.err";
    unless (-s $ae) {
	`cat *.map.err/*.err > $ae 2>/dev/null`;
    }
    foreach my $pat (@torm) {
	foreach my $f (glob($pat)) {
	    if (-d $f) {
		`rm -rf $f`;
	    } else {
		unlink($f);
	    }
	}
    }
    chdir($cwd);
}

# Choose whichever is newer (and exists) of recalibrated fastq file in srcD and
# tgtD. The former is for long-term storage, but a more recent version might have
# arrived in the latter.

sub chooseRecalibratedFastq {
    my ($srcD,$tgtD) = @_;
    my $cand1 = recalFileInRepository($srcD);
    my $cand2 = "$tgtD/recalibrated.fastq.gz";
    return newer($cand2,$cand1) ? $cand2 : $cand1;
}

sub recalFileInRepository {
    my ($srcD) = @_;
    my $f = G1K::G1K::laneFastqLink($srcD);
    if ($f) {
	$f = readlink($f);
	my ($most,$suf) = ($f =~ /^(.*)\.(f[^\/]+)$/);
	$f = "$most.recal.$suf";
	$f .= ".gz" unless ($f =~ /gz$/);
    }
    return -s $f ? $f : undef;
}

sub moveToBad {
    my ($tgtD,$badD) = @_;
    mkdir_p(dirname($badD));
    `rm -rf $badD`; # just in case
    report("moving existing $tgtD to $badD",1);
    rename($tgtD,$badD);
    # Propose to remove higher-level analysis files, i.e.
    # print out commands for doing so, but be cautious and
    # don't actually do it.
    my $super = $tgtD;
    foreach my $i (1,2) {
	$super = dirname($super);
	foreach my $pfx (qw(aln cmd maq_post)) {
	    my $pat = sprintf("%s/%s*",$super,$pfx);
	    if (glob($pat)) {
		print "# rm -f $pat\n";
	    }
	}
    }
}

# Get the reference genome for the specified individual (NAnnnnn)
# and hard-link to it from $tgtD.
sub getRefGenome {
    my ($root,$indiv,$tgtD) = @_;
    my $ans;
    foreach my $suf (qw(fa bfa)) {
	$ans = "$tgtD/ref.$suf";
	unless (-s $ans) {
	    my $tgt = sprintf("%s/ref/human_b36_%s.%s",$root,getGender($root,$indiv),$suf);
	    $tgt = absoluteReadLink($tgt);
	    if (-s $tgt) {
		mkdir_p(dirname($ans));
		link($tgt,$ans);
	    } else { # if ($suf eq 'bfa') {
		die("$tgt: missing or empty");
	    }
	}
    }
    return $ans; # i.e. the bfa one
}

# If .../libN/libParams.txt exists and contains a value for maxInsertLength, return it.
# (I don't think I currently create libParams.txt...)
sub maxInsertLength {
    my ($d) = @_;
    my $f = "$d/libParams.txt";
    if (-s $f) {
	my $fh = openToRead($f);
	while (<$fh>) {
	    chomp;
	    my ($key,$val) = split;
	    if ($key eq 'maxInsertLength') {
		return $val;
	    }
	}
    }
    return undef;
}

# Returns whether the specified (lane) directory is good.
# $op allows us to specified lanes are good for some purposes
# but not others. "Good" means occurrence in goodLanes.txt
# if that exists, else non-occurrence in badLanes.txt if that
# exists, else good. Currently we only used badLanes.txt.

sub inLaneList {
    my ($dir,$op) = @_;
    # Lane name is final component.
    my $lane = basename($dir);
    my $goodFile = "$ENV{G1K}/misc/goodLanes.txt";
    my $badFile = "$ENV{G1K}/misc/badLanes.txt";
    if (-s $goodFile) {
	my $fh = openToRead($goodFile);
	my $tag = "no_$op";
	while (<$fh>) {
	    chomp;
	    my ($lane1,$tag1) = split;
	    return 1 if ($lane1 eq $lane && $tag1 ne $tag);
	}
	return 0; # goodLanes.txt exist and we're not in it
    } elsif (-s $badFile) {
	my $fh = openToRead($badFile);
	while (<$fh>) {
	    chomp;
	    # We allow the items in badLanes.txt to be
	    # (Perl-style) regular expressions, so we can
	    # mark a whole run as bad by saying e.g. lane123_.
	    return 0 if ($lane =~ /^$_$/);
	}
	return 1; # OK -- not bad
    } else {
	return 1; # No goodLanes or badLanes => anything is OK
    }
}

sub getGender {
    my ($root,$indiv) = @_;
    my $gFile = "$root/ref/genders.txt";
    my $fh = openToRead($gFile); # this had better exist!
    while (<$fh>) {
	chomp;
	my ($ind,$gen) = split;
	return $gen if ($ind eq $indiv);
    }
    die("No gender in $gFile for $indiv");
}

# Called by "g1k.pl mrg" action, which should specify as far
# as individual or library, but not lane, because we'll
# merge to this level from the level below.
sub runMaqMerge {
    my ($anaType,$proj,$indiv,$lib) = @_;
    my $root = G1K::G1K::rootDir();

    # jws 2008-10-31
    # srcD comes from DATA in srcAndTgtDirs - i.e. it wants to include anything
    # in data.  For merges, though, I think it is better to include everything
    # in the _target_ tree (although this means you have to be sure it's
    # complete) So, I'm changing srcD to tgtD here.
    
    my ($srcD,$tgtD) = G1K::G1K::srcAndTgtDirs(@_);
    $srcD = $tgtD;

    unless (-d $tgtD) {
	report("$tgtD: does not exist",1);
	return;
    }
    my $cmdF = "$tgtD/cmdMerge";
    my $tgtF = "$tgtD/aln.map";
    my @inFiles;
    # If we're merging to a library from its lanes, we'll
    # read aln.map files, but if from a library to individuals,
    # we expect g1k.pl rmd to have already been run on the
    # libraries, so we'll read from aln-rmdup.map.
    my $inFileBase = defined $lib ? 'aln.map' : 'aln-rmdup.map';
    # Not actually needed here, but created for completeness:
    getRefGenome($root,$indiv,$tgtD);
    # For each subdirectory of $srcD (not $tgtD, because we
    # want toget all the "good" subdirs, not just those that 
    # have been mapped)...
    foreach my $srcSubDir (glob("$srcD/*")) {
	if (-d $srcSubDir) {
	    my $b = basename($srcSubDir);
	    # Find the corresponding subdir in $tgtD.
	    my $subDir = sprintf("%s/%s",$tgtD,$b);
	    my $inFile = "$subDir/$inFileBase";
	    # Only include this dir if it's "good"...
	    if (inLaneList($subDir,'mrg')) {
		# Note we haven't checked for $inFile existing!
		push @inFiles,$inFile;
	    }
	}
    }
    if (@inFiles == 0) {
	report("$tgtD: no subdirectories have $inFileBase files",1);
	return;
    }
    # Check all the inFiles are ready and tgtF can be written.
    my $ms = makeStatus(\@inFiles,[$tgtF]);
    return unless ($ms > 0);
    if (@inFiles == 1) {
	# If we only have one inFile, we don't need to merge, we
	# can just link to it.
	my $bd = basename(dirname($inFiles[0]));
	report("$tgtD: linking aln.map directly to $bd/$inFileBase",1);
	link($inFiles[0],$tgtF);
    } elsif (stamp($tgtF)) {
	# Otherwise, if we can stamp $tgtF, we submit maq mapmerge.
	my $gh = openToWrite($cmdF);
	print $gh "$MAQ mapmerge $tgtD/aln.map @inFiles\n";
	print $gh "rm $tgtF.touch\n";
	$gh->close();
	my $id = bsub($cmdF,"sh $cmdF");
	report("$tgtD: submitted job $id to create aln.map",1);
    }
}

# g1k.pl rmd. Should only be run at library level.
sub runMaqRmdup {
    my ($anaType,$proj,$indiv,$lib) = @_;
    # Find target dir for these arguments.
    my (undef,$tgtD) = G1K::G1K::srcAndTgtDirs(@_);
    # Can and should we derive outF from inF?
    my $inF = "$tgtD/aln.map";
    my $outF = "$tgtD/aln-rmdup.map";
    my $ms = makeStatus([$inF],[$outF]);
    if ($ms > 0 && stamp($outF)) {
	# If so, bsub maq rmdup ...
	my $id = bsub($outF,"$MAQ rmdup $outF $inF; rm $outF.touch\n");
	report("$tgtD: submitted job $id to create aln-rmdup.map",1);
    } else {
	# If already done, see how the statistics look.
	checkRmdupStatus($inF,$outF,@_);
    }
}

sub checkRmdupStatus {
    my ($inF,$outF,@args) = @_;
    my $inSz = -s $inF;
    my $outSz = -s $outF;
    return unless ($inSz && $outSz);
    my $berr = "$outF.berr";
    return unless (-s $berr);
    my $fh = openFromPipe("tail -1 $berr");
    chomp (my $line = <$fh>);
    if ($line =~ /^\[maq_rmdup_core\] (\d+) \/ (\d+) =/) {
	my ($del,$tot) = ($1,$2);
	unless ($tot) {
	    #print "WARNING: last line $line\n  in $berr\n";
	    $tot = 1;
	}
	my $prop = $del/$tot;
	my $szProp = 1-$outSz/$inSz;
	if ($prop > 3 || $szProp < 0.55*$prop || $szProp > 0.85*$prop) {
	    printf("%s has %7.4f%% removed, by file size %7.4f%%\n",
		   join(" ",@args),100*$prop,100*$szProp);
	}
	my $tf = "$outF.touch";
	if (-s $tf && time()-modTime($outF) >= 86400) {
	    printf("Removing %s\n",$tf);
	    unlink($tf);
	}
    } else {
	print "Bad last line in $berr: $line\n";
    }
}

sub looksReady {
    my ($f) = @_;
    return 0 unless (-s $f);
    my @a = stat($f);
    my $age = time()-$a[9];
    return ($age >= 300);
}

sub runMaqPost {
    my ($anaType,$proj,$indiv,$lib) = @_;
    my $root = G1K::G1K::rootDir();
    my (undef,$tgtD) = G1K::G1K::srcAndTgtDirs(@_);
    unless (-d $tgtD) {
	report("$tgtD: does not exist",1);
	return;
    }
    my $cwd = cwd();
    chdir($tgtD);
    my $refGenome = getRefGenome($root,$indiv,$tgtD);
    my $inF = "aln.map";
    my $miscD = "$root/misc";
    # h2s is optional...
    my $h2s = "$miscD/hapmap2/${indiv}-hapmap2.snp.gz";
    my $h3s = "$miscD/hapmap3/${indiv}-hapmap3.snp.gz";
    my $s126 = "$miscD/dbSNP/snp126.snp";
    my $tgtF = "maq_post/cns.flt";
    my @needed = ($refGenome,$inF,$s126);
    my ($h2sSw,$h3sSw);
    if (-s $h2s) {
	push @needed,$h2s;
	$h2sSw = "-2 $h2s";
    }
    push @needed,$h3s;
    $h3sSw = "-3 $h3s";
    my $ms = makeStatus(\@needed,[$tgtF]);
    if ($ms > 0 && stamp($tgtF)) {
	report("$tgtD: submitting jobs",1);
	# indexing the alignment:
	bsub("aln.map.index","$MAQ_NT index -i aln.map");
	# everything else:
	my $cmd = "$MAQ_POST $h2sSw $h3sSw -D 30 -S $s126 maq_post $refGenome $inF";
	print "$cmd\n";
	`$cmd | sh`;
    }
    chdir($cwd);
}

# Validate g1k.pl map outputs. Not essential for further processing, but
# nice for peace of mind.
sub runMaqValidate {
    my ($anaType,$proj,$indiv,$lib,$lane) = @_;
    my (undef,$d) = G1K::G1K::srcAndTgtDirs(@_);
    runMaqValidate1($d,'aln.map');
    runMaqValidate1($d,'aln-rmdup.map') unless ($lane);
}

sub runMaqValidate1 {
    my ($d,$b) = @_;
    my $src = "$d/$b";
    my $tgt = "$d/$b.valid";
    if (makeStatus([$src],[$tgt]) > 0 && stamp($tgt)) {
	my $id = bsub($tgt,"$MAQ mapvalidate $src > $tgt; rm $tgt.touch");
	report("$d: submitted job $id to create $b.valid",1);
    }
}

1;
