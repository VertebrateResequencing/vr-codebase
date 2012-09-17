#!/usr/bin/env perl
# print vertres disk usage
use strict;
use warnings;

print "Type       Mounted on               Size     Used    Avail     Use%\n";

my @disks = qw [/lustre/scratch102 /lustre/scratch105 /lustre/scratch106 ];
push (@disks, `ls -d /nfs/vertres*`);
push (@disks, qw [/warehouse/g1k-01 /warehouse/g1k-02 /warehouse/g1k-03_1 /warehouse/g1k-03_2 /warehouse/g1k-04 ]);

my ($Size_t, $Used_t, $Avail_t, $Use_pct_t);

foreach my $d (@disks) {
  chomp $d;
  $_=`cd $d;df -PTh .| grep -v '^Filesystem'`;

  #ons04-g1k:/g1k_data03/g1k-03 nfs   24T   18T  6.7T  73% /warehouse/g1k-03
  my (undef, $Type, $Size, $Used, $Avail, $Use_pct, $Mounted_on) = split;

  printf ("%-10s %-20s %8s %8s %8s %8s\n", $Type, $Mounted_on, $Size, $Used, $Avail, $Use_pct);

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
