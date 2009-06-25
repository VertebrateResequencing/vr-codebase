#!/usr/bin/perl -w

# Author: lh3
# Last Modified: 2007-12-06

use strict;
use warnings;
use Getopt::Std;

my $version = "0.1.4";

&usage if (@ARGV < 1);
my $command = shift(@ARGV);
my %func = (sub=>\&evalsub, geno=>\&evalgeno, indelpe=>\&evalindelpe, indelsoa=>\&evalindelsoa);
die("Unknown command \"$command\".\n") if (!defined($func{$command}));
&{$func{$command}}();
exit(0);

#
# indelsoa command
#

sub evalindelsoa {
  my %opts = (w=>5, s=>3, m=>1);
  getopts('w:s:m:', \%opts);
  die("
Usage:   maq_eval.pl indelsoa [-w $opts{w}] [-s $opts{s}] <true.indel> <predict.indel>\n
Options: -w INT     window size around the true indel [5]
         -s INT     score of SOA indels [3]\n\n") if (@ARGV < 2);
  my ($fn_true, $fn_pred) = @ARGV;
  my (%ins, %del, $fh);
  my $n_all = &read_true_indel($fn_true, \%ins, \%del, $opts{w});
  my @count = (0, 0);
  open($fh, $fn_pred) || die;
  while (<$fh>) {
	my @t = split;
	next if ($t[4] + $t[5] - $t[3] < $opts{s} || $t[3] > $opts{m});
	my $t1 = $t[1] + $t[2];
	if ($del{$t[0],$t[1]} || $ins{$t[0],$t[1]} || $del{$t[0],$t1} || $ins{$t[0],$t1}) {
	  ++$count[0];
	} else {
	  ++$count[1];
	}
  }
  close($fh);
  my $fn = ($n_all == 0)? -1.0 : (1 - $count[0] / $n_all) * 100.0;
  my $fp = ($count[0] + $count[1] == 0)? -1.0 : (1 - $count[0] / ($count[0] + $count[1])) * 100.0;
  printf qq(
IN: $count[0], OUT: $count[1]
FN: %.2f%%
FP: %.2f%%
), $fn, $fp;
}

#
# indelpe command
#

sub evalindelpe {
  my %opts = (w=>5);
  getopts('w:', \%opts);
  die("
Usage:   maq_eval.pl indelpe [-w $opts{w}] <true.indel> <predict.indel>\n
Options: -w INT     window size around the true indel [5]\n\n") if (@ARGV < 2);
  my ($fn_true, $fn_pred) = @ARGV;
  my (%ins, %del, $fh);
  my $n_all = &read_true_indel($fn_true, \%ins, \%del, $opts{w});

  my %fhash = ("*"=>0, "+"=>1, "-"=>2);
  my @count = (0, 0, 0, 0, 0, 0);
  open($fh, $fn_pred) || die;
  while (<$fh>) {
	my @t = split;
	next if ($t[2] eq ".");
	my $flag = $fhash{$t[2]};
	if ($t[4] < 0) {
	  if ($del{$t[0],$t[1]}) {
		++$count[$flag*2];
	  } else {
		++$count[$flag*2+1];
	  }
	} else {
	  if ($ins{$t[0],$t[1]}) {
		++$count[$flag*2];
	  } else {
		++$count[$flag*2+1];
	  }
	}
  }
  close($fh);
  my $n_called = $count[0] + $count[1] + $count[2] + $count[3];
  my $fn = ($n_all == 0)? -1.0 : ($n_all - $count[0] - $count[1]) / $n_all * 100.0;
  my $fp = ($n_called == 0)? -1.0 : ($count[1] + $count[3]) / $n_called * 100.0;
  print "\n";
  print "Type '*' - IN: $count[0], OUT: $count[1]\n";
  print "Type '+' - IN: $count[2], OUT: $count[3]\n";
  print "Type '-' - IN: $count[4], OUT: $count[5]\n\n";
  printf "*/+ FN: %.2f%%\n", $fn;
  printf "*/+ FP: %.2f%%\n", $fp;
}

sub read_true_indel {
  my ($fn, $ins, $del, $window) = @_;
  $window ||= 5;
  my $fh;
  my ($n_ins, $n_del) = (0, 0);
  open($fh, $fn) || die;
  while (<$fh>) {
	my @t = split;
	if ($t[2] eq "-") {
	  ++$n_ins;
	  for (my $i = $t[1] - $window; $i <= $t[1] + $window; ++$i) {
		$ins->{$t[0],$i} = 1;
	  }
	}
	if ($t[3] eq "-") {
	  ++$n_del;
	  for (my $i = $t[1] - $window; $i <= $t[1] + $window; ++$i) {
		$del->{$t[0],$i} = 1;
	  }
	}
  }
  close($fh);
  print "true insertions: $n_ins\n";
  print "true deletions:  $n_del\n";
  return $n_ins + $n_del;
}

#
# sub command
#

sub evalsub {
  die("Usage: maq_eval.pl sub [-e <.err>] [-g] -p <prefix> <dbSNP.snp> <genotyping.snp> <.snp>\n") if (@ARGV < 4);
  my (%opts, $fh);
  getopts('p:e:gs', \%opts);
  my $is_summary = (defined($opts{s}))? 1 : 0;
  my $is_gnuplot = (defined($opts{g}))? 1 : 0;
  my $err_out = ($opts{e})? $opts{e} : '';
  my $prefix = ($opts{p})? $opts{p} : "es-$$";
  my (%ref_hash, %good_hash, %test_hash);
  &read_snp_file($ARGV[0], \%ref_hash);
  &read_snp_file($ARGV[1], \%good_hash);
  &read_snp_file($ARGV[2], \%test_hash, 1);

  my (@n_hit, @n_in, @n_iden, @n_q, @n_het, @cover);
  my (@n_hom10, @n_het10, @n_homhit10, @n_hethit10, @n_cover10);
  for (my $i = 0; $i != 256; ++$i) { $n_hit[$i] = $n_in[$i] = $n_iden[$i] = $n_q[$i] = $n_het[$i] = $cover[$i] = 0; }
  for (my $i = 0; $i <= 9; ++$i) { $n_hom10[$i] = $n_het10[$i] = $n_homhit10[$i] = $n_hethit10[$i] = $n_cover10[$i] = 0; }

  # calculate how many candidate SNPs are in dbSNPs
  foreach my $pos (keys %test_hash) {
	$test_hash{$pos} =~ /(\S),(\d+)$/;
	my $q10 = int($2/10);
	my $is_het;
	$is_het = ($1 ne 'A' && $1 ne 'C' && $1 ne 'G' && $1 ne 'T')? 1 : 0;
	++$n_q[$2];
	if ($is_het) {
	  ++$n_het[$2]; ++$n_het10[$q10];
	} else { ++$n_hom10[$q10]; }
	if (defined($ref_hash{$pos})) {
	  ++$n_hit[$2];
	  if ($is_het) { ++$n_hethit10[$q10]; }
	  else { ++$n_homhit10[$q10]; }
	}
  }
  # calculate how many ExoSeq SNPs are covered/correct.
  foreach my $pos (keys %good_hash) {
	if (defined($test_hash{$pos})) {
	  $test_hash{$pos} =~ /(\S),(\d+)/;
	  ++$n_in[$2];
	  ++$n_iden[$2] if ($1 eq $good_hash{$pos});
	}
  }
  # accumulate
  for (my $i = 98; $i >= 0; --$i) {
	$n_hit[$i] += $n_hit[$i+1];
	$n_in[$i] += $n_in[$i+1];
	$n_iden[$i] += $n_iden[$i+1];
	$n_q[$i] += $n_q[$i+1];
	$n_het[$i] += $n_het[$i+1];
  }
  # calculate the percetage
  for (my $i = 0; $i <= 99; ++$i) {
	if ($n_q[$i] == 0) {
	  $n_hit[$i] = 100;
	  $n_het[$i] = 100;
	} else {
	  $n_hit[$i] /= $n_q[$i] * .01;
	  $n_het[$i] /= $n_q[$i] * .01;
	}
	if (scalar(keys %good_hash) == 0) {
	  $n_in[$i] = 100;
	  $n_iden[$i] = 100;
	} else {
	  $n_in[$i] /= scalar(keys %good_hash) * .01;
	  $n_iden[$i] /= scalar(keys %good_hash) * .01;
	}
  }
  # read ".err"
  if ($err_out) {
	my $tot_len;
	open($fh, $err_out) || die;
	while (<$fh>) {
	  if (/S0 reference length: (\d+)/) {
		$tot_len = $1;
	  } elsif (/S0 number of gaps in the reference: (\d+)/) {
		$tot_len -= $1;
	  } elsif (/^S1\s+(\d+)\s+(\d+)\s+(\d+)/) {
		$cover[$1] = $3 / $tot_len * 100.0;
		$n_cover10[int($1/10)] += $2;
	  }
	}
	close($fh);
  }
  # print out to ".txt"
  open($fh, ">$prefix.txt");
  for (my $i = 99; $i >= 0; --$i) {
	my $tmp = ($n_in[$i] == 0)? 100 : $n_iden[$i]/$n_in[$i]*100;
	printf $fh ("$i\t$n_q[$i]\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", $n_hit[$i], $n_in[$i], $tmp,
				$n_het[$i], $cover[$i]);
  }
  close($fh);
  # print out to ".fri"
  open($fh, ">$prefix.fri");
  printf $fh ("Q %12s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%12s\n", "coverage", "n_hom", "homHit", "homFrac",
			  "n_het", "hetHit", "hetFrac", "total", "totHit", "totFrac", "totCum", "covCum");
  my ($totSum, $covSum) = (0, 0);
  for (my $i = 9; $i >= 0; --$i) {
	$totSum += $n_hom10[$i]+$n_het10[$i];
	$covSum += $n_cover10[$i];
	my $tmp1 = $n_hom10[$i] > 0? 100.0*$n_homhit10[$i]/$n_hom10[$i] : 100;
	my $tmp2 = $n_het10[$i] > 0? 100.0*$n_hethit10[$i]/$n_het10[$i] : 100;
	my $tmp3 = $n_hom10[$i]+$n_het10[$i] > 0? ($n_homhit10[$i]+$n_hethit10[$i])/($n_hom10[$i]+$n_het10[$i])*100.0 : 100;
	printf $fh ("%dx%12d%8d%8d%7.2f%%%8d%8d%7.2f%%%8d%8d%7.2f%%%8d%12d\n", $i, $n_cover10[$i], $n_hom10[$i], $n_homhit10[$i],
				$tmp1, $n_het10[$i], $n_hethit10[$i], $tmp2,
				$n_hom10[$i]+$n_het10[$i], $n_homhit10[$i]+$n_hethit10[$i],
				$tmp3, $totSum, $covSum);
  }
  close($fh);
  # plot
  if ($is_gnuplot) {
	my $fn = "$prefix.txt";
	my $y2skip = int($n_q[0]/10+0.5);
	open($fh, "| gnuplot") || die;
	print $fh qq(
set t po eps so co;
set xlab "PHRED quality";
set ylab "Percentage (Cumulated)";
set y2lab "Cumulated Number of SNPs";
set out "$prefix.eps";
set size 1.2,1.2;
set key left bottom;
set yran [0:100];
set ytics 10;
set y2ran [0:$n_q[0]];
set y2tics $y2skip;
set grid;
plot "$fn" u 1:3 t "in_dbSNP%" w l lw 7, "$fn" u 1:4 t "GenoTyping-cover%" w l lw 7, "$fn" u 1:6 t "het%" w l, "$fn" u 1:7 t "coverage%" w l, "$fn" u 1:2 t "n_SNPs" w l ax x1y2, "$fn" u 1:5 t "GenoTyping-exact%" w l;
exit;
\n);
	close($fh);
  }
}

sub evalgeno {
  my %opts = (Q=>40, d=>3, D=>254, q=>20);
  my %ishom = (A=>1, C=>1, G=>1, T=>1, a=>1, c=>1, g=>1, t=>1);
  my %test_hash;
  getopts("Q:q:D:d:", \%opts);
  die("
Usage:   maq_eval.pl geno [options] <in.geno> <in.snp>\n
Options: -q INT      minimum consensus quality [$opts{q}]
         -Q INT      minimum mapping quality [$opts{Q}]
         -d INT      minimum read depth [$opts{d}]
         -D INT      maximum read depth [$opts{D}]
\n") if (@ARGV < 2);
  &read_snp_file($ARGV[1], \%test_hash);
  my ($fh, $n_total, $n_seg, @count);
  $n_total = $n_seg = 0;
  open($fh, $ARGV[0]) || die;
  for my $x (0..2) { for (0..3) { $count[$x][$_] = 0; } }
  while (<$fh>) {
	my @t = split;
	++$n_total;
	my $row = ($t[2] eq $t[10])? 0 : (($ishom{$t[10]})? 1 : 2);
	my $is_filt = ($t[5] >= $opts{d} && $t[5] <= $opts{D} && $t[6] >= 0.25
				   && $t[7] >= $opts{Q} && $t[4] >= $opts{q})? 0 : 1;
	$is_filt = 1 if ($t[2] ne $t[3] && !$test_hash{$t[0],$t[1]}); # correct $is_filt
	my $is_missing = ($row && $t[2] eq $t[3])? 1 : 0;
	my $col = $is_filt? 0 : ($is_missing? 1 : (($t[3] eq $t[10])? 3 : 2));
	++$n_seg if ($row);
	$count[$row][$col]++;
  }
  close($fh);
  $count[0][0] = ">=$count[0][0]";
  printf "\n# genotyped sites: $n_total\n";
  printf "# segregating sites: $n_seg\n\n";
  printf "%11s%11s%11s%11s%11s\n", "", "filtered", "missing", "wrong", "identical";
  my @label = ("true_mono", "true_hom", "true_het");
  for my $x (0..2) {
	printf "%11s", $label[$x];
	for my $y (0..3) { printf "%11s", $count[$x][$y]; }
	print "\n";
  }
  print "\n";
}

sub read_snp_file {
  my ($fn, $hash2, $is_qual) = @_;
  my $fh;
  open($fh, $fn) || die;
  while (<$fh>) {
	my @t = split;
	next if ($t[2] eq '-' || $t[3] eq '-'); # skip indel sites
	$t[4] = 99 if ($t[4] =~ /^\d+$/ && $t[4] > 99);
	$hash2->{$t[0],$t[1]} = ($is_qual)? "$t[3],$t[4]" : $t[3];
  }
  close($fh);
}

sub usage {
  die qq(
Program: maq_eval.pl (evaluate Maq SNPs/indels)
Version: $version
Contact: Heng Li <lh3\@sanger.ac.uk>\n
Usage:   maq_eval.pl <command> <arguments>\n
Options: sub        evaluate substitutions
         geno       evaluate genotyping results
         indelpe    evaluate indelpe results
         indelsoa   evaluate indelsoa results
\n);
}
