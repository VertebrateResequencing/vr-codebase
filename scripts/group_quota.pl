#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Path::Class;
use Getopt::Long;

my $default_groups_file = file($ENV{CONF}, 'vertres_unix_groups.conf')->stringify;

my $opts = parse_params();
print_quotas($opts);

exit;

#--------------------------------


sub error
{
    my (@msg) = @_;
    if ( scalar @msg )
    {
        croak join('',@msg);
    }
    die
        "About: Print lustre group quotas.\n",
        "Usage: group_quota.pl [OPTIONS] <lustre_disk> [...]\n",
        "Options:\n",
        "   -g, --groups <file>              File listing unix groups for which to check quotas [$default_groups_file].\n",
        "   -h, -?, --help                   This help message.\n",
        "\n";
}

sub parse_params
{
    my $opts = {};
    $$opts{groups} = file($ENV{CONF}, 'vertres_unix_groups.conf')->stringify;
    while (my $arg=shift(@ARGV))
    {
        if ( $arg eq '-g' || $arg eq '--groups' ) { $$opts{groups} = shift(@ARGV); next }
        if ( $arg eq '-?' || $arg eq '-h' || $arg eq '--help' ) { error(); }
        push @{$$opts{disks}}, $arg;
    }
    error() unless @{$$opts{disks}};
    return $opts;
}


sub print_quotas {
    my ($self) = @_;

    open my $gfh, $$self{groups} or die "$!: $$self{groups}";
    my @groups;
    foreach my $group (<$gfh>)
    {
        chomp $group;
        push @groups, $group;
    }
    close $gfh;

    foreach my $disk (@{$$self{disks}})
    {
        chomp $disk;
        next unless ($disk =~ /lustre/);
        next unless (-d $disk);

        print "disk                 group            disk_used      disk_limit      files      file_limit\n";
        print "------------------------------------------------------------------------------------------\n";

        my ($disk_used_total,$disk_limit_total,$files_total,$file_limit_total);
        my (%result, %usage);
        foreach my $group (@groups)
        {
            chomp $group;
            my $res = `sudo /software/lustre_operator/bin/lustre_operator /usr/bin/lfs $disk getquota -g $group | tail -1`;
            chomp $res;
            # filesystem  type    name    size-used   size-softlimit  size-hardlimit  size-remaining  size-grace  inode-used  inode-softlimit inode-hardlimit inode-remaining inode-grace timestamp
            # /lustre/scratch116  group   team145 10827334881280  unlimited   329853488332800 319026153451520 inactive    1098903 unlimited   10000000    8901097 inactive    2015-08-20T11:48:54
            my (undef,undef,undef,$disk_used,undef,$disk_limit,undef,undef,$files,undef,$file_limit,undef) = split /\s+/, $res;
            next unless $disk_used;
            if ($disk_limit eq "unlimited" || $file_limit eq "unlimited") {
                $result{$group} = sprintf ("%-20s %-15s %10s       unlimited %10d       unlimited\n", $disk, $group, displayK($disk_used), $files);
                $usage{$group} = 1000;
            }
            else {
                my $disk_pct = $disk_limit ? ($disk_used/$disk_limit) * 100 : 0;
                my $file_pct = $file_limit ? ($files/$file_limit) * 100 : 0;
                $result{$group} = sprintf ("%-20s %-15s %10s %10s %3.0f%% %10d %10d %3.0f%%\n", $disk, $group, displayK($disk_used), displayK($disk_limit), $disk_pct, $files, $file_limit, $file_pct);
                $usage{$group} = $disk_pct > $file_pct ? $disk_pct : $file_pct;
            }
            $disk_used_total += $disk_used;
            $disk_limit_total += $disk_limit unless ($disk_limit eq "unlimited");
            $files_total += $files;
            $file_limit_total += $file_limit unless ($file_limit eq "unlimited");
        }
        foreach my $group (sort { $usage{$b} <=> $usage{$a} } keys(%usage) ) {
            print $result{$group};
        }
        my $disk_pct_total = $disk_limit_total ? ($disk_used_total/$disk_limit_total) * 100 : 0;
        my $file_pct_total = $file_limit_total ? ($files_total/$file_limit_total) * 100 : 0;
        print "------------------------------------------------------------------------------------------\n";
        printf ("%-20s %-15s %10s %10s %3.0f%% %10d %10d %3.0f%%\n\n", "", "TOTAL", displayK($disk_used_total), displayK($disk_limit_total), $disk_pct_total, $files_total, $file_limit_total, $file_pct_total);
    }
}

sub displayK
{
    my $num = shift;
    my $sfx = '';

    if ($num > 1000) {
      $num /= 1024;
      $sfx = 'K';
    }
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
