package G1K::CopyContents;

use Exporter;
@ISA = qw(Exporter);

@EXPORT = qw(copyContents);

use G1K::G1K;
use Utility;
use File::Basename;
use Cwd;

use strict;

# This code is to copy maq outputs from Heng's directory
# into the analysis hierarchy. I don't think we will ever need it now.

my $HENG_DIR = '/lustre/scratch1/lh3/reseq*'; # can be a pattern ...

sub copyContents {
    my ($anaType,$proj,$indiv,$lib,$run) = @_;
    my $srcBase;
    if ($run =~ /^lane(\d+)_(\d)$/) {
	$srcBase = sprintf("%s_s_%s",$1,$2);
    } elsif ($run =~ /lane(.*)$/) {
	my $runId = $1;
	if ($lib =~ /^libBI/) {
	    $srcBase = 'BI*' . substr($runId,3);
	    $srcBase =~ s/_(.)$/.\1/;
	} else {
	    $srcBase = $runId;
	}
    } else {
	report("$run: not in form lane*",1);
	return;
    }
    my @srcDirs = glob("$HENG_DIR/*/$srcBase");
    my (undef,$tgtDir) = G1K::G1K::srcAndTgtDirs(@_); # why qualified??
    if (@srcDirs == 1) {
	unless (-d $tgtDir) {
	    report("$tgtDir: creating from $srcDirs[0]");
	    copyContents1($srcDirs[0],$tgtDir);
	    report();
	}
    } elsif (@srcDirs == 0) {
	if (-d $tgtDir) {
	    copyContents1($tgtDir,$tgtDir);
	} else {
	    report("$HENG_DIR/*/$srcBase: not found",1);
	}
    } else {
	report(sprintf("%s/*/%s: %d matches",$HENG_DIR,$srcBase,scalar(@srcDirs)));
    }
}

# If srcD eq tgtD, remove all files other than the @keep ones; else
# copy from srcD to tgtD.
sub copyContents1 {
    my ($srcD,$tgtD) = @_;
    mkdir_p($tgtD);
    my $mE = "$tgtD/aln.map.err";
    unless (-s $mE) {
	my @meList = glob("$srcD/*map.err/*");
	`cat $srcD/*map.err/* > $mE` if (@meList);
    }
    my @keep = qw(aln.map aln-rmdup.log aln-rmdup.mapcheck aln-rmdup.mapstat unmap.txt.gz);
    if ($srcD eq $tgtD) {
	my %toKeep = map {$_ => 1} ('aln.map.err',@keep);
	my @torm;
	my $cwd = cwd();
	chdir($srcD);
	foreach my $f (glob("*")) {
	    unless ($toKeep{$f} || $f =~ /(bout|berr)$/) {
		push @torm,$f;
	    }
	}
	if (@torm) {
	    `rm -rf @torm`;
	    report(sprintf("%s: removed %d files",$tgtD,scalar(@torm)),1);
	}
	chdir($cwd);
    } else {
	my $sameFS = sameFilesystem($srcD,$tgtD);
	foreach my $b (@keep) {
	    copyOrLink($b,$srcD,$tgtD,$sameFS);
	}
    }
}

# Return 2 if same dir, 1 if same filesystem but not same dir, else 0.

sub sameFilesystem {
    my ($dir1,$dir2) = @_;
    return 2 if ($dir1 eq $dir2);
    my $fh1 = openFromPipe("df $dir1");
    <$fh1>;
    $_ = <$fh1>;
    my ($sys1) = split;
    my $fh2 = openFromPipe("df $dir2");
    <$fh2>;
    $_ = <$fh2>;
    my ($sys2) = split;
    return ($sys1 eq $sys2) ? 1 : 0;
}

sub copyOrLink {
    my ($b,$srcD,$tgtD,$sameFS) = @_;
    my $srcF = "$srcD/$b";
    return unless (-f $srcF);
    my $tgtF = "$tgtD/$b";
    return if (-s $tgtF == -s $srcF);
    if ($sameFS) {
	link($srcF,$tgtF);
    } else {
	`cp -p $srcF $tgtF`;
    }
}

1;
