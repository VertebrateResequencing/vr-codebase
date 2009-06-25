#!/usr/bin/perl -w

# Version: 16Jul2008

use strict;
use warnings;
use Getopt::Std;
use Cwd qw/abs_path getcwd/;

my %opts = (P=>'maq', p=>'', q=>undef);
getopts('P:sp:q:', \%opts);
die("Usage: map2qMapFile.pl [-P $opts{P}] [-s] [-p null] <input.map> <ref.bfa>\n") if (@ARGV < 2);
my $maq = gwhich($opts{P}) || die("Cannot find maq executable.\n");
die("The input alignment file name must be ended with '.map'") unless ($ARGV[0] =~ /\.map$/);
$ARGV[0] =~ s/\.map$//;
#&print_qMapFile(1, "$ARGV[0].1.mapcheck"); exit; # for testing only
$opts{p} = "-P $opts{p}" if ($opts{p});
$opts{q} = "-q $opts{q}" if defined ($opts{q});
if (defined $opts{s}) {
  warn("Run mapcheck...\n");
  system("$maq mapcheck -cS1 $opts{p} $opts{q} $ARGV[1] $ARGV[0].map > $ARGV[0].mapcheck");
  warn("Convert to qMapFile...\n");
  &print_qMapFile(0, "$ARGV[0].mapcheck");
} else {
  warn("Break pair...\n");
  system("$maq breakpair $ARGV[0]");
  warn("Run mapcheck...\n");
  system("$maq mapcheck -cS1 $opts{p} $opts{q} $ARGV[1] $ARGV[0].1.map > $ARGV[0].1.mapcheck");
  system("$maq mapcheck -cS1 $opts{p} $opts{q} $ARGV[1] $ARGV[0].2.map > $ARGV[0].2.mapcheck");
  warn("Convert to qMapFile...\n");
  &print_qMapFile(1, "$ARGV[0].1.mapcheck");
  &print_qMapFile(2, "$ARGV[0].2.mapcheck");
}

sub print_qMapFile {
  my ($which, $fn) = @_;
  my ($stage, $cycle) = (0, 0);
  my $fh;
  open($fh, $fn) || die;
  while (<$fh>) {
	if (/AC\s+AG\s+AT\s+CA/) {
	  $stage = 1;
	} if ($stage == 1 && /^\s*[\s\d.]+:\s+[\d\s]+:\s+(\d+.*)/) {
	  my @t = split(/\s+/, $1);
	  my $max_q = (@t-1) / 2 - 1;
	  ++$cycle;
	  for (0..$max_q) {
		my ($c1, $c2) = @t[$_, 2+$max_q+$_];
		my $err = ($c2+0.1/($max_q+1)) / ($c1+0.1);
		print join("\t", $which, $cycle, $_, int(-4.343*log($err)+0.5), $c1, $c2), "\n";
	  }
	}
  }
  close($fh);
}

# the following codes are copied from treefam::generic

sub dirname
{
	my $prog = shift;
	my $cwd = getcwd;
	return $cwd if ($prog !~ /\//);
	$prog =~ s/\/[^\s\/]+$//g;
	return $prog;
}
sub which
{
	my $file = shift;
	my $path = (@_)? shift : $ENV{PATH};
	return if (!defined($path));
	foreach my $x (split(":", $path)) {
		$x =~ s/\/$//;
		return "$x/$file" if (-x "$x/$file" && -f "$x/$file");
	}
	return;
}
sub gwhich
{
	my $progname = shift;
	my $addtional_path = shift if (@_);
	my $dirname = &dirname($0);
	my $tmp;

	chomp($dirname);
	if (-x $progname && -f $progname) {
		return abs_path($progname);
	} elsif (defined($addtional_path) && ($tmp = &which($progname, $addtional_path))) {
		return $tmp; # lh3: Does it work? I will come back to this later
	} elsif (defined($dirname) && (-x "$dirname/$progname" && -f "$dirname/$progname")) {
		return abs_path("$dirname/$progname");
	} elsif (($tmp = &which($progname))) { # on the $PATH
		return $tmp;
	} else {
		warn("[gwhich] fail to find executable $progname anywhere.");
		return;
	}
}
