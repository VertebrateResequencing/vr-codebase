package Utility;

### DMC's toolkit of incredibly useful everyday functions ###

use Exporter;
@ISA = qw(Exporter);

@EXPORT = qw(
  absoluteFileName absoluteReadLink alias allBinaries allFilesReady
  appendBeforeGZ arch await estimatedFinishTime fileToHash fileToList
  filesFromSmallest hashToFile identicalFiles isComplete listToFile max
  min mkdir_p moveAside newer nonZeroBinaries openFromPaste openFromPipe
  openToAppend openToPipe openToRead openToWrite readAliases
  replaceBeforeGZ report sameList stamp touch unstamp readFastaFile
  writeFastaFile avoidInitialDigit printGFF readQuad makeStatus modTime
);

use FileHandle;
use Sys::Hostname;
use File::Basename;
use Cwd;
use Carp;
use File::Spec; # for rel2abs
use Sys::Hostname;

use strict;

############## Wrappers around FileHandle methods #######################

# openToRead() is like new FileHandle(), but it correctly handles
# gzipped files, and return *STDIN if there's no file.

sub openToRead {
    my ($f) = @_;
    if (! $f) {
	#die("Empty filename to read from\n");
	return *STDIN; # no arg means standard input
    } elsif ($f =~ /\.gz$/ || ($f =~ /\.gz\./ && $f !~ /\.gz\.touch$/)) {
	# allow .gz.sv1 etc ...
	return new FileHandle("gunzip -c $f |");
    } elsif (-s "$f.gz" && ! -s $f) {
	openToRead("$f.gz");
    } else {
	openSafely($f,"r","read from");
    }
}

# openFromPipe($cmd) is like new FileHandle("$cmd |")

sub openFromPipe {
    my ($exp) = @_;
    $exp =~ s/\s+$//;
    $exp .= "|" unless (substr($exp,-1) eq "|");
    return openSafely($exp,undef,"pipe from");
}

# openToPipe($cmd) is like new FileHandle("| $cmd")

sub openToPipe {
    my ($exp,@files) = @_;
    $exp =~ s/^\s+//;
    $exp =~ s/^\|//; # delete pipe at front
    return new FileHandle("| $exp");
}

# openToWrite is like new FileHandle(">...") except that it
# correctly handles gzipped files, and returns *STDOUT if
# there's no file.

sub openToWrite {
    my ($g) = @_;
    return *STDOUT unless ($g);
    my $dir = dirname($g);
    mkdir_p($dir);
    if ($g =~ /gz$/) {
	return new FileHandle("| gzip -c > $g");
    } else {
	return openSafely($g,"w","write to");
    }
}

# openToAppend is like new FileHandle(...,'a')

sub openToAppend {
    my ($f) = @_;
    openSafely($f,"a","append to");
}

# If the open fails, report the fact properly...
sub openSafely {
    my ($f,$mode,$verb) = @_;
    if ($f) {
	if ($mode) {
	    return new FileHandle($f,$mode) || carp ("Cannot $verb $f");
	} else {
	    return new FileHandle($f) || carp("Cannot $verb $f");
	}
    } else {
	return *STDOUT;
    }
}

# Paste the @files together and open the result for reading.

sub openFromPaste {
    my (@files) = @_;
    foreach my $f (@files) {
	croak("openFromPaste: cannot find $f\n") unless (-f $f);
    }
    return openFromPipe("paste @files");
}

########### Read/write files between hashes, lists and aliases #############

# Read the contents of $file into $pHash, using field $keyField
# (counting from 0) of each line as the key and field $valField as the
# value.

sub fileToHash {
    my ($file,$pHash,$keyField,$valField) = @_;
    $keyField = 0 if ($keyField eq "");
    $valField = 1 if ($valField eq "");
    my $maxField = max($valField,$keyField);
    my $fh = openToRead($file);
    while (<$fh>) {
	chomp;
	s/^\s+//;
	if (/^[^#]/) { # skip lines starting with hash
	    my @arr = split;
	    if ($#arr >= $maxField) {
		# If there are enough fields in @arr, set $pHash->{key} = val
		$pHash->{$arr[$keyField]} = $arr[$valField];
	    }
	}
    }
}

# Write a hash to a file: each line is "$key $value"

sub hashToFile {
    my ($pHash,$file,$mode) = @_;
    my $gh = openSafely($file,$mode,"write to");
    foreach my $key (sort keys %$pHash) {
	printf($gh "%s %s\n",$key,$pHash->{$key});
    }
    $gh->close();
}

# Read $file to list $pList; like fileToHash, but $keyField
# is for the index of the list.

sub fileToList {
    my ($file,$pList,$keyField,$valField) = @_;
    $keyField = 0 if ($keyField eq "");
    $valField = 1 if ($valField eq "");
    my $maxField = $valField;
    $maxField = $keyField if ($keyField > $maxField);
    my $fh = openToRead($file);
    while (<$fh>) {
	chomp;
	if (/^[^#]/) {
	    my @arr = split;
	    if ($#arr >= $maxField && $arr[$keyField] =~ /^\d+$/) {
		$pList->[$arr[$keyField]] = $arr[$valField];
	    }
	}
    }
}

# Write a list ot a file, one (non-empty) item per line.

sub listToFile {
    my ($pList,$file,$mode) = @_;
    my $gh = openSafely($file,$mode,"write to");
    foreach my $key (0 .. $#$pList) {
	if ($pList->[$key] ne "") {
	    printf($gh "%d %s\n",$key,$pList->[$key]);
	}
    }
    $gh->close();
}

# Used in SGRP project.

sub readAliases {
    my ($f) = @_;
    return undef unless (-s $f);
    my %alias;
    my $fh = openToRead($f);
    my ($multi,$var);
    while (<$fh>) {
	chomp;
	unless (/^\s*\#/) {
	    next unless (/\S/);
	    if (/^\s/) {
		s/^\s+//;
		my @vals = split;
		push @{$alias{$var}},@vals;
		$multi = 1;
	    } else {
		my @vals;
		($var,@vals) = split;
		if (defined $var && @vals) {
		    my %h;
		    my $isHash = 1;
		    foreach my $val (@vals) {
			if ($val =~ /^(.*)=(.*)$/) {
			    $h{$1} = $2;
			} else {
			    $isHash = 0;
			    last;
			}
		    }
		    if ($isHash) {
			$alias{$var} = \%h;
		    } else {
			$alias{$var} = \@vals;
		    }
		    $multi = 1 if ($isHash || @vals > 1);
		}
	    }
	}
    }
    unless ($multi) {
	foreach my $var (keys %alias) {
	    $alias{$var} = $alias{$var}[0];
	}
    }
    return \%alias;
}

sub alias {
    my ($var,$pAli) = @_;
    return defined $pAli->{$var} ? $pAli->{$var} : $var;
}

#### Touchfile code ####

# Create an empty file (not used much if at all)

sub touch {
    my ($f) = @_;
    my $fh = new FileHandle($f,"w");
    if ($fh) {
	$fh->close();
    } else {
	croak("Can't touch $f\n");
    }
}

# Wait for a file to exist, checking every $secs seconds.

sub await {
    my ($file,$secs,$repDir) = @_;
    my $repFile = $file;
    $repFile = "$repDir/$file" if ($repDir);
    $secs = 60 unless ($secs > 0);
    my $nWaits = 0;
    while (-f "$file.touch" || ! -f $file) {
	if ($nWaits >= 10) {
	    my $gh = openToAppend($ENV{'HOME'} . "/await.log");
	    printf($gh "Awaiting %s at %s in command %s\n",
		   $repFile,scalar(localtime),
		   join(" ",$0,@ARGV));
	    $nWaits = 0;
	}
	print "Awaiting $repFile ...\n";
	sleep($secs);
	$nWaits ++;
    }
}

# Return 1 iff all files exist and don't still have an $f.touch file.

sub allFilesReady {
    my (@files) = @_;
    foreach my $f (@files) {
	return 0 if (-f "$f.touch" || ! -s $f);
    }
    return 1;
}

# Try to "stamp" (claim) file $f. If either $f or $f.touch
# exists, fail. Else write our identity (hostname and process)
# to $f.touch, wait 3 seconds and make sure we can read back
# the same identity (this minimises race conditions but
# can't altogether eliminate them...but this never seems
# to fail!)

sub stamp {
    my ($f,$app) = @_;
    my $t = "$f.touch";
    # Fail if it's already there:
    return undef if (-f $t || ((-f $f || -d $f) && !$app));
    # id is my hostname and process ID; should be unique...
    my $id = sprintf("%s.%d",hostname(),$$);
    my $gh = openToWrite($t);
    print $gh "$id\n";
    $gh->close();
    # Wait for any competitors
    sleep(3);
    return undef unless (-s $t);
    my $fh = openToRead($t);
    chomp ($_ = <$fh>);
    # Only succeed if we didn't get overwritten
    return ($id eq $_);
}

# Remove $f.touch, croaking if it's not there, unless $nullOK.

sub unstamp {
    my ($f,$nullOK) = @_;
    my $t = "$f.touch";
    croak("missing $t") unless ($nullOK || -f $t);
    unlink($t);
    return 1;
}

# Return 1 if $f exists and $f.touch doesn't.

sub isComplete {
    my ($f,$rep) = @_;
    if (! -f $f) {
	report("$f: missing",1) if ($rep);
	return 0;
    }
    if (-f "$f.touch") {
	report("$f.touch: still exists",1) if ($rep);
	return 0;
    }
    return 1;
}

############ Estimate finish time for iterative process #############

sub estimatedFinishTime {
    # Started at time() = $t0; done $k out of $n iterations
    my ($t0,$k,$n) = @_;
    my $t1 = time();
    my $t2 = $t1; # wild guess for now
    if ($k > 0) {
	$t2 = $t1 + ($t1-$t0)*($n-$k)/$k;
    }
    return scalar(localtime($t2));
}

############# Progress reporter #####################

my $tReport; # static variable ...

# Call at the beginning and end of an operation like this, to time it:
#
#   report($msg) => print date+time and $msg, then "..."
#   report()     => print "done" and how long since report($msg)
#
# Call to report an event at a specific time but with no "done" to follow:
#
#   report($msg,1)

sub report {
    my ($msg,$nl,$gh) = @_;
    $gh = *STDOUT unless ($gh);
    if ($msg && $msg ne 'failed') {
	$tReport = time();
	printf($gh "%s: %s",scalar(localtime($tReport)),replaceEnvironmentVariables($msg));
	if ($nl) {
	    print $gh "\n";
	} else {
	    print $gh " ... ";
	}
    } else {
	my $d = time()-$tReport;
	my ($dStr,$units);
	if ($d >= 3600) {
	    $dStr = sprintf("%d:%2d:%2d",int($d/3600),int(($d-3600*int($d/3600))/60),$d-60*int($d/60));
	    $units = "hours";
	} elsif ($d >= 60) {
	    $dStr = sprintf("%d:%2d",int($d/60),$d-60*int($d/60));
	    $units = "minutes";
	} elsif ($d > 0) {
	    $dStr = $d;
	    $units = "seconds";
	} else {
	    $dStr = "no";
	    $units = "time";
	}
	$dStr =~ s/ /0/g;
	$msg = "done" unless ($msg);
	print $gh "$msg in $dStr $units\n";
    }
}

# Replace value of $G1K by "$G1K" and same for $SGR (and potentially
# other variables...). This shortens messages for easy reading.

sub replaceEnvironmentVariables {
    my ($msg) = @_;
    foreach my $v (qw(G1K SGR)) {
	my $d = $ENV{$v};
	if ($d) {
	    my $i = index($msg,$d);
	    if ($i >= 0) {
		substr($msg,$i,length($d)) = '$' . $v;
	    }
	}
    }
    return $msg;
}

############# min, max and list-comparison functions ###############

# Minimum of two numbers; undef counts as +infinity.

sub min {
    my ($x,$y) = @_;
    if (defined $y && $x > $y) {
	return $y;
    } elsif (defined $x) {
	return $x;
    } else {
	return $y;
    }
}
   
# Maximum of two numbers; undef counts as -infinity.

sub max {
    my ($x,$y) = @_;
    if (defined $y && $y > $x) {
	return $y;
    } elsif (defined $x) {
	return $x;
    } else {

	return $y;
    }
}

# Succeed if the two lists are the same length and element-wise eq.

sub sameList {
    my ($pL1,$pL2) = @_;
    return 0 unless ($#$pL1 == $#$pL2);
    foreach my $k (0 .. $#$pL1) {
	return 0 unless ($pL1->[$k] eq $pL2->[$k]);
    }
    return 1;
}

############# File and directory manipulation #############

# Move $f to the first of $f.sv1, $f.sv2, ... that
# doesn't exist yet. (cf my "sv -m" shell script)

sub moveAside {
    my ($f) = @_;
    return unless (-f $f || -d $f);
    return if ($f =~ /\.sv\d+$/); # already an .sv file
    my $n=1;
    while (-s "$f.sv$n" || -d "$f.sv$n") {
	$n++;
    }
    rename($f,"$f.sv$n");
    return "$f.sv$n";
}    

# Return whether $f1 is newer than $f2. A non-existent file
# is infinitely old.

sub newer {
    my ($f1,$f2) = @_;
    return 0 unless (-f $f1);
    return 1 unless (-f $f2);
    #croak("missing comparator in newer: $f2") unless (-f $f2);
    return 0 if ($f1 eq $f2); # equal is not newer!
    return (modTime($f1) > modTime($f2));
}

# "mkdir -p": mkdir, but create directories en route when necessary.

sub mkdir_p {
    my ($dir) = @_;
    $dir = absoluteFileName($dir);
    unless (-d $dir) {
	my $dn = dirname($dir);
	mkdir_p($dn);
	mkdir($dir);
    }
}

# Return the absolute name of a (relative) file.

sub absoluteFileName {
    my ($f) = @_;
    return File::Spec->rel2abs(File::Spec->canonpath($f));
}

# readlink until we get to a non-link, then absoluteFileName.

sub absoluteReadLink {
    my ($f) = @_;
    if (-l $f) {
	return absoluteReadLink(readlinkAbs($f));
    } else {
	return absoluteFileName($f);
    }
}

sub readlinkAbs {
    my ($f) = @_;
    return $f unless (-l $f);
    my $r = readlink($f);
    unless (substr($r,0,1) eq '/') {
	$r = sprintf("%s/%s",dirname($f),$r);
    }
    return $r;
}

# Return 1 if $f and $g both exist and are identical in contents.

sub identicalFiles {
    my ($f,$g) = @_;
    return 0 unless (-f $f && -f $g);
    return 0 if (-s $f != -s $g); # quicker to check sizes, maybe...
    my $fh = openFromPipe("cmp $f $g 2>&1");
    my $msg = <$fh>;
    return !$msg;
}

############# All patterns of 0 and 1 of length n ###############

sub nonZeroBinaries {
    my ($n) = @_;
    my @ans = allBinaries($n);
    shift @ans; # to get rid of the 0000 one
    return @ans;
}

sub allBinaries {
    my ($n) = @_;
    if ($n == 0) {
	return ("");
    } else {
	my @tails = allBinaries($n-1);
	my @ans;
	foreach my $k (0,1) {
	    foreach my $tail (@tails) {
		push @ans,"$k$tail";
	    }
	}
	return @ans;
    }
}
		
###############

# Return the value of "arch" (or uname -m, more generally)

sub arch {
    my $fh = openFromPipe("uname -m"); # arch not defined on all machines
    $_ = <$fh>;
    chomp;
    return $_;
}

# Add $suf to a file name, before any final .gz but otherwise
# at the end.

sub appendBeforeGZ {
    my ($f,$suf) = @_;
    if ($f =~ /\.gz$/) {
	$f =~ s/\.gz$/${suf}.gz/;
    } else {
	$f .= $suf;
    }
    return $f;
}

# Change $from to $to at the end of $f, but allow .gz to be
# right at the end.

sub replaceBeforeGZ {
    my ($f,$from,$to) = @_;
    if ($f =~ /${from}\.gz$/) {
       $f =~ s/${from}\.gz$/${to}.gz/;
    } elsif ($f =~ /${from}$/) {
	$f =~ s/${from}$/$to/;
    } else {
	croak("Cannot find $from at end of $f");
    }
    return $f;
}

# Return the files in increasing order of size.

sub filesFromSmallest {
    my (@files) = @_;
    my %sz;
    foreach my $f (@files) {
	$sz{$f} = -s $f;
    }
    return sort {$sz{$a} <=> $sz{$b}} @files;
}

# Return a hash of names and sequences from a fasta file.
# The name of a sequence is the characters after the ">" but
# before any space.

sub readFastaFile {
    my ($f) = @_;
    my $fh = openToRead($f);
    my (%seq,$buf,$chr);
    while (<$fh>) {
	chomp;
	if (/^>(\S+)/) {
	    $seq{$chr} = $buf if ($buf);
	    $chr = $1;
	    $buf = "";
	} else {
	    $buf .= $_;
	}
    }
    $seq{$chr} = $buf if ($buf);
    return \%seq;
}

# Given a hash from names to sequences, write 
# the specified fastq file with line length 60.

sub writeFastaFile {
    my ($pSeq,$g) = @_;
    my $gh = openToWrite($g);
    foreach my $chr (sort keys %$pSeq) {
	my $seq = $pSeq->{$chr};
	print $gh ">$chr\n";
	my $len = length($seq);
	for (my $i=0; $i<$len; $i+=60) {
	    print $gh substr($seq,$i,60) . "\n";
	}
    }
    $gh->close();
}

# SGRP only.

sub avoidInitialDigit {
    my ($str) = @_;
    $str =~ s/^(\d+)/S\1/;
    my $lstr = $str;
    $lstr =~ tr/A-Z/a-z/;
    return ($str,$lstr);
}

# SGRP only.

sub printGFF {
    my ($gh,@arr) = @_;
    if (@arr > 9) {
	confess("Overlong GFF line: @arr");
    }
    unless (defined $gh) {
	confess("Undefined output handle for @arr");
    }
    foreach my $i (0 .. 8) {
	if (!defined $arr[$i]) {
	    $arr[$i] = '.';
	} elsif ($arr[$i] eq 'F') {
	    $arr[$i] = '+';
	} elsif ($arr[$i] eq 'C') {
	    $arr[$i] = '-';
	}
    }
    # Ensure 2nd field doesn't start with a digit
    #$arr[1] = avoidInitialDigit($arr[1]);
    print $gh join("\t",@arr) . "\n";
}

# Intended for Maq-style Fastq files, which consist of
# a header line starting in "@", one or more lines of DNA letters,
# a header line starting in "+", then the same number of quality
# characters as we had DNA letters. This is an unambiguous spec
# because "+" is not a DNA symbol.
#
# This should work for normal (one line per sequence or quality
# string) Fastq files too.

sub readQuad {
    my ($fh) = @_;
    my @a;
    # Set nucs header line
    my $line = <$fh>;
    return undef unless ($line);
    die("Expected line starting in @, got: $line") unless (substr($line,0,1) eq '@');
    chomp ($a[0] = $line);
    report("reading quad: $a[0]");
    # Read nucs
    while (1) {
	chomp (my $line = <$fh>);
	if (substr($line,0,1) eq '+') {
	    $a[2] = $line;
	    last;
	} else {
	    $a[1] .= $line;
	}
    }
    # Read same number of quals as we have nucs
    my $wanted = length($a[1]);
    while ($wanted > 0) {
	chomp (my $line = <$fh>);
	$wanted -= length($line);
	$a[3] .= $line;
    }
    report();
    return \@a;
}

# makeStatus returns the "status" of a ("make") operation that
# would create the files in @$pOutFiles from the inputs in
# @$pInFiles. The status values returned are:
#
#   2 -> need to make; all outfiles missing
#   1 -> need to make; some outfiles missing or out of date
#   0 -> no need to make; all infiles present [and non-empty], all outfiles present and in date,
#        or one has .touch so already being made
#  -1 -> cannot make; some infiles missing [or empty]
#
# Options: z => infiles don't have to be non-empty, just present
#          n => DO NOT remove any existing outfiles when status == 1
#          a => min age in seconds of infiles (default 0)

sub makeStatus {
    my ($pInFiles,$pOutFiles,$pOpt) = @_;
    my $inTime;
    my ($zok,$rmo,$recent,$rep) = (undef,1,undef,1);
    if ($pOpt) {
	$zok = $pOpt->{z};
	$rmo = !$pOpt->{n};
	$recent = time()-$pOpt->{a};
	$rep = !$pOpt->{q};
    }
    foreach my $f (@$pInFiles) {
	my $af = $f;
	$af = absoluteFileName($f);
	if (-f "$f.touch") { # $f still being made
	    report("$af: still being made",1) if ($rep);
	    return -1;
	}
	$f = readlinkAbs($f);
	my $t = modTime($f,$zok);
	unless ($t > 0) {
	    report("$af: not found",1) if ($rep);
	    return -1;
	}
	if (defined $recent && $t > $recent) {
	    report("$af: too recent",1);
	    return -1;
	}
	$inTime = max($inTime,$t);
    }
    my $outTime;
    my @presOut;
    foreach my $f (@$pOutFiles) {
	return 0 if (-f "$f.touch"); # $f already being made
	$f = readlinkAbs($f);
	my $t = modTime($f,1);
	if ($t > 0) {
	    $outTime = min($outTime,$t);
	    push @presOut,$f if ($rmo);
	}
    }
    if (!defined $outTime) {
	return 2;
    } elsif ($outTime < $inTime) {
	unlink @presOut if (@presOut);
	return 1;
    } else {
	return 0;
    }
}

# Return the time (in seconds since 1970) at which $f was modified, or
# zero if it's missing or (unless $zok) empty.

sub modTime {
    my ($f,$zok) = @_;
    return 0 unless ($zok ? -f $f : -s $f);
    my @a = stat($f);
    return $a[9];
}    

# This is probably only needed for the old farm ...

sub END {
    my $tmpd = sprintf("/tmp/%s.%s",$ENV{USER},$ENV{LSB_JOBID});
    if (-d $tmpd) {
	sleep(15); # to give separate lsrcp processes etc time to finish
	`rm -rf $tmpd`;
    }
}

1;
