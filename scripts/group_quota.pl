#!/usr/bin/env perl
# print vertres group disk quotas
use strict;
use warnings;
use Path::Class;

# my $disks = $ARGV[0] ? $ARGV[0] : file($ENV{CONF}, 'vertres_disks.conf')->stringify;
# open my $dfh, $disks or die "$!: $disks";
# my @disks;
# foreach my $disk (<$dfh>)
# {
#     chomp $disk;
#     next unless ($disk =~ /lustre/);
#     next unless (-d $disk);
# }
# close $dfh;

my @disks = @ARGV;

# my $groups = $ARGV[0] ? $ARGV[0] : file($ENV{CONF}, 'vertres_unix_groups.conf')->stringify;
# open my $gfh, $groups or die "$!: $groups";
# my @groups;
# foreach my $group (<$gfh>)
# {
#     chomp;
#     push @groups, $group;
# }
# close $gfh;

my @groups = split /\s+/, `groups`;
chomp $groups[-1];

foreach my $disk (@disks)
{
    chomp $disk;
    next unless ($disk =~ /lustre/);
    next unless (-d $disk);

    print "disk                 group            disk_used      disk_limit      files      file_limit\n";
    print "------------------------------------------------------------------------------------------\n";

    my ($disk_used_total,$disk_limit_total,$files_total,$file_limit_total);
    foreach my $group (@groups)
    {
        chomp $group;
        my $res = `lfs quota -g $group $disk | tail -1`;
        $res =~ s/^\s+//;
        chomp $res;
        my ($disk_used,undef,$disk_limit,undef,$files,undef,$file_limit,undef) = split /\s+/, $res;
        next unless $disk_used;
        my $disk_pct = $disk_limit ? ($disk_used/$disk_limit) * 100 : 0;
        my $file_pct = $file_limit ? ($files/$file_limit) * 100 : 0;
        printf ("%-20s %-15s %10s %10s %3.0f%% %10d %10d %3.0f%%\n", $disk, $group, displayK($disk_used), displayK($disk_limit), $disk_pct, $files, $file_limit, $file_pct);
        $disk_used_total += $disk_used;
        $disk_limit_total += $disk_limit;
        $files_total += $files;
        $file_limit_total += $file_limit;
    }
    my $disk_pct_total = $disk_limit_total ? ($disk_used_total/$disk_limit_total) * 100 : 0;
    my $file_pct_total = $file_limit_total ? ($files_total/$file_limit_total) * 100 : 0;
    print "------------------------------------------------------------------------------------------\n";
    printf ("%-20s %-15s %10s %10s %3.0f%% %10d %10d %3.0f%%\n\n", "", "TOTAL", displayK($disk_used_total), displayK($disk_limit_total), $disk_pct_total, $files_total, $file_limit_total, $file_pct_total);
}

sub displayK
{
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
