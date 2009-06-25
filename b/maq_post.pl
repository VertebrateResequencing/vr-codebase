#!/usr/bin/perl -w

# Author: lh3

use strict;
use warnings;
use Cwd qw/getcwd/;
use Getopt::Std;
use File::Spec;

my $version = '0.1.2';
my %opts = (2=>'', 3=>'', S=>'', Q=>60, n=>20, d=>3, D=>256, q=>10);
getopts('2:3:S:Q:n:d:D:q:', \%opts);
die(qq(
Program: maq_post.pl (run maq jobs on Sanger Inst's computing farm)
Contact: lh3
Version: $version
Usage:   maq_post.pl [options] <run.dir> <ref.bfa> <aln.map>

Options: -2 FILE   HapMap-1+2 SNPs [null]
         -3 FILE   HapMap-3 SNPs [null]
         -S FILE   dbSNP in .snp format [null]
         -Q INT    minimum mapping quality [$opts{Q}]
         -d INT    minimum read depth [$opts{d}]
         -n INT    minimum neighbouring quality [$opts{n}]
         -D INT    maximum read depth [$opts{D}]
         -q INT    minimum consensus quality [$opts{q}]

)) if (@ARGV < 3);
my $maq = gwhich('maq') || die("[maq_post] fail to find 'maq'");
my $maq_pl = gwhich('maq.pl') || die("[maq_post] fail to find 'maq.pl'");
my $maq_eval_pl = gwhich('maq_eval.pl') || die("[maq_post] fail to find 'maq_eval.pl'");
my $asub = gwhich('asub') || die("[maq_post] fail to find 'asub'");
my $outdir = shift(@ARGV);
my $ref = File::Spec->rel2abs(shift(@ARGV));
die("[maq_post] fail to stat reference file '$ref'\n") unless (-e $ref);
my $aln = File::Spec->rel2abs(shift(@ARGV));
die("[maq_post] fail to stat alignment file '$aln'\n") unless (-e $aln);
$_ = $outdir;
s/.*\///;
my $jn = "$_.$$";
my $cwd = getcwd;

# check files

for ('2', '3', 'S') {
  die("[maq_post] fail to find '$opts{$_}'.\n") if ($opts{$_} && !-e $opts{$_});
  $opts{$_} = File::Spec->rel2abs($opts{$_}) if ($opts{$_});
}

print qq(mkdir -p $outdir; cd $outdir;\n);

# make genotype file

if ($opts{2} && $opts{3}) {
  print qq(test ! -e geno.raw.snp && awk 'BEGIN{while((getline<"$opts{2}")>0)l[\$1","\$2]=\$4}l[\$1","\$2]{print \$1,\$2,\$3,\$4,l[\$1","\$2]}' $opts{3}|tr " " "\\t" > geno.raw.snp;
awk '\$4==\$5' geno.raw.snp > geno.snp;\n);
} elsif ($opts{2}) {
  print qq(cp $opts{2} geno.snp;\n);
} elsif ($opts{3}) {
  print qq(cp $opts{3} geno.snp;\n);
}

# compose the command lines

my $filt = qq/-Q $opts{Q} -n $opts{n} -d $opts{d} -D $opts{D}/;

print qq(
$asub -j "$jn.core" <<EOF
  test ! -e cns.cns      && $maq assemble cns.cns $ref $aln 2> cns.cns.log
  test ! -e cns.indelsoa && $maq indelsoa $ref $aln > cns.indelsoa
  test ! -e cns.indelpe  && $maq indelpe  $ref $aln > cns.indelpe
  test ! -e aln.mapcheck && $maq mapcheck $ref $aln > aln.mapcheck
  test ! -e aln.mapstat  && $maq mapstat  $aln > aln.mapstat
EOF
$asub -j "$jn.geno" -w "done($jn.core)" <<EOF
  test -f 'geno.snp' && test ! -e cns.geno && $maq subpos cns.cns geno.snp > cns.geno
EOF
$asub -j "$jn.snpreg" -w "done($jn.core)" <<EOF
  test ! -e cns.snpreg && $maq snpreg $filt cns.cns 2> cns.snpreg
EOF
$asub -j "$jn.cns2snp" -w "done($jn.core)" <<EOF
  test ! -e cns.snp && $maq cns2snp cns.cns > cns.snp
EOF
$asub -j "$jn.filt" -w "done($jn.cns2snp)" <<EOF
  $maq_pl SNPfilter $filt -f cns.indelsoa -F cns.indelpe cns.snp > cns.flt
EOF
$asub -j "$jn.evalgeno" -w "done($jn.filt) && done($jn.geno)" <<EOF
  test -f 'geno.snp' && $maq_eval_pl geno -q $opts{q} $filt cns.geno cns.flt > cns.gev
EOF
test -e '$opts{S}' && test -e 'geno.snp' && bsub -o /dev/null -e /dev/null -J "$jn.evalsub" -w "done($jn.filt)" -R"select[mem>3000] rusage[mem=3000]" -W 3000000 <<EOF
  $maq_eval_pl sub -e cns.cns.log -p cns $opts{S} geno.snp cns.flt
EOF
bsub -w "done($jn.evalsub) && done($jn.evalgeno)" -o /dev/null -e /dev/null 'rm -f cns.flt.touch'
);
print "cd $cwd;\n";

exit;

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
		return File::Spec->rel2abs($progname);
	} elsif (defined($addtional_path) && ($tmp = &which($progname, $addtional_path))) {
		return $tmp; # lh3: Does it work? I will come back to this later
	} elsif (defined($dirname) && (-x "$dirname/$progname" && -f "$dirname/$progname")) {
		return File::Spec->rel2abs("$dirname/$progname");
	} elsif (($tmp = &which($progname))) { # on the $PATH
		return $tmp;
	} else {
		warn("[gwhich] fail to find executable $progname anywhere.");
		return;
	}
}
