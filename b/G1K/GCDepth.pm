package G1K::GCDepth;

use Exporter;
@ISA = qw(Exporter);

@EXPORT = qw(gcDepth);

use G1K::G1K;
use Utility;
use Bsub;
use File::Basename;
use Cwd;

use strict; 

# The specific contents of this file come from Aylwyn (as6), but as it's
# relatively simple, it's a good template for building new g1k.pl actions.
# Here's a check list:
#
# * Choose a 3-letter action name, e.g. foo
# * Add example to the ## comments in G1K.pm
# * Add a tie to a function in %CMD in G1K.pm, e.g. foo => \&doFooFunc
# * Add an entry in %NFIXED if the number of arguments (after g1k.pl foo) is fixed
# * Add "use G1K::DoFoo" to top of G1K.pm
# * Create G1K/DoFoo.pm, e.g. from this file
# * Define doFooFunc and export it.

# Where to find useful commands like "maq".

my $BIN = $ENV{G1K} . "/bin";

# All functions tied in %CMD get arguments like these; typical values are
#    v1 Trio-CEU NA12878 libSC_NA12878_1 lane354_1
# When the user specified patterns on the command line, the function will be
# called multiple times with different arguments.

sub gcDepth {
    my ($anaType,$proj,$indiv,$lib,$lane) = @_;
    # Get the source directory (under $G1K/DATA) and the target
    # (under $G1K/$anatype/...). We don't always need the source.
    my (undef,$dir) = G1K::G1K::srcAndTgtDirs(@_);
    # Check this lane isn't "bad"; return if it is.
    return unless (G1K::RunMaq::inLaneList($dir));
    my $src = "$dir/aln.map";
    my $tgt = "$dir/depth/depth.txt";
    my $root = G1K::G1K::rootDir();
    G1K::RunMaq::getRefGenome($root,$indiv,$dir); # ensure it's there
    my $ms = makeStatus([$src],[$tgt]);
    if ($ms > 0 && stamp($tgt)) {
	report("$tgt: creating");
	my $cwd = cwd();
	mkdir("$dir/depth");
	chdir("$dir/depth");
	my $binSize = $dir =~ /LowCov/ ? 500 : 50;
	`$BIN/maq mapview ../aln.map | $BIN/mapdepth ../ref.fa -b=$binSize > bindepth-$binSize`;
	`$BIN/gcdepth.py bindepth-$binSize -i $indiv > $tgt`;
	chdir($cwd);
	unstamp($tgt);
	report();
    }
}

1;
