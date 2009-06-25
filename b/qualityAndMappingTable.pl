#!/usr/bin/env perl

use Utility;
use strict;

# Create quality-and-mapping-quantity spreadsheet for Richard.

autoflush STDOUT 1;

print "Project,Population,Individual,Library,Lane,PercentMapped,InputQuality,RecalQuality,TotalSeq,EstMappableSeq,AdjustedDepth1,AdjustedDepth2\n";
foreach my $d (glob("$ENV{G1K}/v1/*-???/NA?????/lib*/lane*")) {
    my $stat = "$d/sample/aln.stat";
    my $qm = "$d/qualmapBayesian.txt";
    # checks on recalibrated are to ensure qm file is OK...
    if (-s $stat && -s $qm && -s "$d/recalibrated.fastq.gz" >= 100000 && ! -s "$d/recalibrated.fastq.gz.touch") {
	printQualityAndMapping($d,$stat,$qm);
    }
}

sub printQualityAndMapping {
    my ($d,$stat,$qm) = @_;
    my ($lane,$lib,$indiv,$proj) = reverse(split(/\//,$d));
    $proj =~ s/-/,/;
    $proj =~ s/Trio/Pilot2/;
    $proj =~ s/LowCov/Pilot1/;
    $lane =~ s/^lane//;
    $lib =~ s/^lib//;
    my $perc = percentMapped($stat);
    my ($qIn,$qOut) = meanQualities($qm);
    my @vals = ($proj,$indiv,$lib,$lane,$perc,sprintf("%3.1f",$qIn),sprintf("%3.1f",$qOut));
    (my $srcD = $d) =~ s/\/v1\//\/DATA\//;
    my ($fq) = glob("$srcD/*f*[qz]");
    if (-l $fq) {
	my $reposBase = readlink($fq);
	$reposBase =~ s/\.gz$//;
	$reposBase =~ s/fq$/fastq/;
	my $fqc = $reposBase .  "check";
	my $totalSeq;
	if (-s $fqc) {
	    my $fh = openToRead($fqc);
	    my $line = <$fh>;
	    my (undef,undef,$bp) = split(" ",$line);
	    $totalSeq = $bp/1000000;
	} else {
	    my $fh = openToRead($fq);
	    <$fh>;
	    chomp (my $seq = <$fh>);
	    my $len = length($seq);
	    $fh->close();
	    my $cmd = $fq =~ /gz$/ ? 'zcat' : 'cat';
	    chomp (my $lines = `$cmd $fq | wc -l`);
	    $lines =~ s/[^0-9]//g;
	    $totalSeq = $lines*$len/(4*1000000);
	}
	my $ems = $totalSeq*$perc/100;
	my $ad1 = $ems/2800;
	my $ad2 = $ad1*$qOut/25;
	foreach my $v ($totalSeq,$ems,$ad1,$ad2) {
	    push @vals,sprintf("%5.3f",$v);
	}
	print join(",",@vals) . "\n";
    }
}

sub percentMapped {
    my ($f) = @_;
    my $fh = openToRead($f);
    while (<$fh>) {
	chomp;
	if (/mapped [SP]E reads:/) {
	    s/.* //;
	    if (/%\)$/) {
		s/..$//;
		return sprintf("%4.2f",$_) if ($_ > 0);
	    }
	}
    }
    return undef;
}

sub meanQualities {
    my ($f) = @_;
    my $fh = openToRead($f);
    my ($sum2,$sum1,$sum);
    while (<$fh>) {
	chomp;
	my ($r,$p,$q1,$q2,$n) = split;
	$sum1 += $q1*$n;
	$sum2 += $q2*$n;
	$sum += $n;
    }
    return ($sum1/$sum,$sum2/$sum);
}
