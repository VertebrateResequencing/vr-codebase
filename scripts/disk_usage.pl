#!/usr/bin/env perl
# print vertres disk usage
use strict;
use warnings;

my $inputfile = $ARGV[0] ? $ARGV[0] : '/lustre/scratch105/conf/vertres_disks.conf';
open FILE, $inputfile or die "$!: inputfile";

print "Type       Mounted on               Size     Used    Avail     Use%   OST Max%\n";
my ($Size_t, $Used_t, $Avail_t, $Use_pct_t);

foreach my $d (<FILE>) {
  chomp $d;
  $_=`cd $d;df -PTh .| grep -v '^Filesystem'`;

  #ons04-g1k:/g1k_data03/g1k-03 nfs   24T   18T  6.7T  73% /warehouse/g1k-03
  my (undef, $Type, $Size, $Used, $Avail, $Use_pct, $Mounted_on) = split;

  my $ost_max='';  
  if ($d =~ /lustre/) {
    $ost_max=`lfs df $d | sort -nr -k 5 | grep -v 'filesystem summary' | head -1 | awk '{print \$5}'`;
    chomp $ost_max;
  }

  printf ("%-10s %-20s %8s %8s %8s %8s %5s\n", $Type, $Mounted_on, $Size, $Used, $Avail, $Use_pct, $ost_max);

  # Totals; note df always rounds up, hence any discrepancies in addition
  $_=`cd $d;df -PT .| grep -v '^Filesystem'`;
  (undef, $Type, $Size, $Used, $Avail, $Use_pct, $Mounted_on) = split;
  $Size_t += $Size;
  $Used_t += $Used;
  $Avail_t += $Avail;
}

$Use_pct_t = ($Used_t/$Size_t) * 100;
printf ("%-10s %-20s %8s %8s %8s %7.0f%%\n", "", "Total", displayK($Size_t), displayK($Used_t), displayK($Avail_t),  $Use_pct_t);

sub displayK {
  my $num = shift;
  my $sfx = 'K';

  if ($num > 1000) {
    $num /= 1024;
    $sfx = 'M';
  }
  if ($num > 1000) {
    $num /= 1024;
    $sfx = 'G';
  }
  if ($num > 1000) {
    $num /= 1024;
    $sfx = 'T';
  }

  return sprintf("%7.0f%s", $num, $sfx);
}
