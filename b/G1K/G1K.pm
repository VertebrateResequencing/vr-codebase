package G1K::G1K;

use Exporter;
@ISA = qw(Exporter);

@EXPORT = qw(g1k g1kCommand rootDir srcAndTgtDirs getReadLength laneFastqLink);

use G1K::RunMaq;
use G1K::Calibrate;
use G1K::GCDepth;
use G1K::Cleanup;

use Utility;
use Bsub;
use File::Basename;
use Getopt::Std;

use Carp;

use strict;

autoflush STDOUT 1;

# Special comment block, printed out when command is invoked without
# arguments. The sample commands should give an instance of each
# available (or commonly-used) action, in the right order.

##
##Example commands:
##
##  g1k.pl sam v1 LowCov-CEU NA12156 SLX libSC_NA12156_1 lane445_7
##  g1k.pl mca v1 LowCov-CEU NA12156 SLX libSC_NA12156_1 lane445_7
##  g1k.pl cal v1 LowCov-CEU NA12156 SLX libSC_NA12156_1 lane445_7
##
##  g1k.pl map v1 LowCov-CEU NA12156 SLX libSC_NA12156_1 lane445_7
##  g1k.pl mrg v1 LowCov-CEU NA12156 SLX libSC_NA12156_1
##  g1k.pl rmd v1 LowCov-CEU NA12156 SLX libSC_NA12156_1
##  g1k.pl mrg v1 LowCov-CEU NA12156 SLX
##  g1k.pl pst v1 LowCov-CEU NA12156 SLX
##  g1k.pl gcd v1 LowCov-CEU NA12156 SLX
##
##  g1k.pl cln v1

# These are the commands (actions) that come immediately after g1k.pl,
# and the functions they call. There is a special command "all" too,
# defined in g1k code below (ugh).

my %CMD = (# Recalibration. sam = collect a (random) sample of 30m
           # reads and put them in <laneDir>/sample/30m.fastq.
           # mca = run maq as part of calibration.
           # cal = actually calibrate: generate quality map and apply
           # it to the uncalibrated fastq file to create
           # <laneDir>/recalibrated.fastq.gz. See Calibrate.pm.
           sam => \&sampleFastq,
	   mca => \&runMaqMapCal,
           cal => \&calibrate,

	   # Main mapping/analysis sequence. See RunMaq.pm.
	   map => \&runMaqMap,   # Run maq map on a lane of data
	   mrg => \&runMaqMerge, # Merge lanes to libs or libs to indivs.
	   rmd => \&runMaqRmdup, # Remove duplicates (at lib level)
	   pst => \&runMaqPost,  # Post-processing: call SNPs etc
           # This is a good example to start from when someone asks you
           # to implement a new g1k.pl action...
           gcd => \&gcDepth,    # GC-dependent depth analysis (GCDepth.pm)
	   # Cleaning up, e.g. before a backup
	   # jws: 2008-12-05: this is inadvisable right now.
	   #cln => \&cleanup);

# Fixed numbers of arguments.
# jws: 2008-12-05 changed to include the tech level, e.g. SLX
my %NFIXED = (cpy => 6,
	      pst => 4,
	      gcd => 4,
              cal => 6
	      ); # and default is -1 for G1K

# Note rootDir => $ENV{G1K} below. So your G1K environment variable
# must be set to the root of the analysis hierarchy: /lustre/sf4/1kgenomes/G1K 
# on sflogin, /lustre/thougen/thougen on cbi4 (if you ever run there).
sub g1k {
    my ($pOpt,@args) = @_;
    # G1K-specific information to be passed to expandAndSubmit, which
    # is general.
    my %info = (nFixed => \%NFIXED,
		nFixedDefault => -1,
		command => \%CMD,
		commandSet => {},
		decide => \&decideBsubParametersG1K,
		rootDir => $ENV{G1K},
		shellCommand => 'g1k.pl',
		instructionModule => findInstructionModule(),
		expansionRoot => "$ENV{G1K}/DATA",
		callback => \&g1k);
    # To allow directory paths to be pasted in, allow components
    # to be separated by "/" as an alternative to " ".
    @args = splitAtSlashes(@args);
    # Define the 'all' action. It should be "g1k.pl all v1"
    # followed by up to five arguments.
    if ($args[0] eq 'all' && @args >= 2) {
	# Action
	my $ac = $args[1];
	# Extend with dots (equivalent to shell "*" i.e. anything)
	while (@args < 7) {
	    push @args,'.';
	}
	expandAndSubmit($pOpt,\%info,'map',@args[1..6]);
	expandAndSubmit($pOpt,\%info,'mrg',@args[1..5]);
	expandAndSubmit($pOpt,\%info,'rmd',@args[1..5]);
	expandAndSubmit($pOpt,\%info,'mrg',@args[1..4]);
	expandAndSubmit($pOpt,\%info,'pst',@args[1..4]);
    } else {
	expandAndSubmit($pOpt,\%info,@args);
    }
}

# Where to look for the special (##) comment block (see above):
sub findInstructionModule {
    my $im;
    foreach my $d (@INC) {
	$im = "$d/G1K/G1K.pm";
	return $im if (-s $im);
    }
    return undef;
}

sub splitAtSlashes {
    my (@args) = @_;
    my @ans;
    foreach my $arg (@args) {
	if (substr($arg,0,1) eq '/') {
	    # absolute: not safe to split
	    push @ans,$arg;
	} else {
	    $arg =~ s/\/$//;
	    push @ans,split("/",$arg);
	}
    }
    return @ans;
}

# Given a command and its arguments, return parameters for bsub:
# the prefix of the command to give to it, the queue, and the
# requirements.
sub decideBsubParametersG1K {
    my ($cmd,@args) = @_;
    my $req;
    if (substr($ENV{G1K},0,8) eq '/lustre/') {
	$req = 'select[lustre]';
    }
    my $fixedSwitches = "-r $ENV{G1K}";
    my $fixedCmdPfx = "g1k.pl $fixedSwitches $cmd";
    my $queue = 'normal';
    return ($fixedCmdPfx,$queue,$req);
}

# Callback
sub rootDir {
    return $ENV{G1K};
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

sub srcAndTgtDirs {
    my @args = @_;
    
    my $srcD = join("/",$ENV{G1K},'DATA',@args[1..5]);
    my $tgtD = join("/",$ENV{G1K},@args[0..5]);
    my $badD = join("/",$ENV{G1K},$args[0] . '_BAD',@args[1..5]);
    $srcD =~ s/\/+$//;
    $tgtD =~ s/\/+$//;
    return ($srcD,$tgtD,$badD);
}

sub getReadLength {
    my ($dir) = @_;
    # $rlFile should be created by updateDataDirectory.pl.
    my ($rlFile) = glob("$dir/*readlengths");
    if (-s $rlFile) {
	my $fh = openToRead($rlFile);
	chomp ($_ = <$fh>);
	my ($len1,$len2) = split;
	return ($len1,$len2==0);
    }
    # No readlengths file -- have to look in fastq[check].
    # The suffix of fastq files varies between .fastq, .fastq.gz
    # and .fq.gz.
    my $pairLen;
    foreach my $fqFile (glob("$dir/*fastq*"),glob("$dir/*fq*")) {
	$fqFile = readlink($fqFile) if (-l $fqFile);
	# $fqc should be the fastqcheck file. This always ought
	# to exist for Sanger lanes, and Tom Skelly can create them
	# for other sites' lanes on request.
	my $fqc = $fqFile . 'check';
	$fqc =~ s/\.gzcheck$/check/;
	if (-s $fqc) {
	    # Pair length is 6th field on 1st line.
	    my $fh = openToRead($fqc);
	    chomp ($_ = <$fh>);
	    my @a = split();
	    $pairLen = int($a[5]);
	    last;
	} elsif (-s $fqFile) {
	    # Open the fastq file itself, and assume all reads
	    # are the same length.
	    my $fh = openToRead($fqFile);
	    <$fh>; # header
	    chomp ($_ = <$fh>); # nucleotide sequence
	    $pairLen = length($_);
	    last;
	}
    }
    unless (defined $pairLen) {
	die("Cannot find readlengths or fastq files in $dir");
    }
    if ($pairLen < 60 || $pairLen%2 > 0) {
	# Assume it's a single-end run if "pair" length <= 60 or is odd.
        # This assumption is NOT future proof, which is why we
	# prefer to have readlengths files!
	report("$dir: assuming length $pairLen is for single end",1);
	return ($pairLen,1);
    } else {
	report("$dir: assuming length $pairLen is for paired end",1);
	return ($pairLen/2,0);
    }
}


sub laneFastqLink {
    my ($d) = @_;
    foreach my $suf (qw(fastq fastq.gz fq.gz)) {
	my ($link) = glob("$d/*.$suf");
	if ($link && -l $link) {
	    return $link;
	}
    }
    return undef;
}

1;
