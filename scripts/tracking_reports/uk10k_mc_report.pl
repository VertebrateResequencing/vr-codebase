#!/usr/bin/env perl 
# pulls stats for UK10K, totals + prev month.

use strict;
use warnings;
no warnings 'uninitialized';
use VertRes::Utils::VRTrackFactory;
use Getopt::Long;
use Data::Dumper;

my ($year, $month, $help);

GetOptions(
    'y|year=i'  => \$year,
    'm|month=i' => \$month,
    'h|help'    => \$help,
    );

my ($c_sec,$c_min,$c_hour,$c_mday,$c_mon,$c_year,$c_wday,$c_yday,$c_isdst) = localtime(time);

$year ||= $c_year+1900; # current year
$month ||= $c_mon;    # last month, as localtime is zero-based, mysql is 1-based

(!$help) or die <<USAGE;
    Usage: $0   
                [--year  <year for month to report on. Defaults current>]
                [--month <month to report on.  Defaults to last month.>]
                [--help <this message>]

Outputs summary stats with totals + "last month".  Last month is defined by year and month params.

USAGE


my @exome_yield;
my @wgs_yield;

print "Last month is $month/$year\n";

foreach my $db('vrtrack_uk10k_cohort','vrtrack_uk10k_neuro','vrtrack_uk10k_rare','vrtrack_uk10k_obesity'){
    my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $db,
                                                          mode => 'r');
    $vrtrack or die "Can't connect to $vrtrack";
                                                        
    my $dbh = $vrtrack->{_dbh};
    my $samples_with_reqs = get_samples_with_requests($dbh);
    my $samples_with_reqs_lm = get_samples_with_requests_last_month($dbh,$year,$month);
    my $libs = get_libraries($dbh);
    my $libs_lm = get_libraries_last_month($dbh,$year,$month);
    my $bgi  = get_bgi($dbh);
    my $samples_with_seq = get_samples_with_seq($dbh);
    my $samples_with_seq_lm = get_samples_with_seq_last_month($dbh,$year,$month);
    
    my $gbp = $db eq 'vrtrack_uk10k_cohort' ? 12 : 4;
    my $finished = get_finished($dbh,$gbp);
    my $finished_lm = get_finished_last_month($dbh,$gbp,$year,$month);

    if ($db eq 'vrtrack_uk10k_cohort'){
        push @wgs_yield, @{get_lane_yield_last_month($dbh,$year,$month)};
    }
    else {
        push @exome_yield, @{get_lane_yield_last_month($dbh,$year,$month)};
    }
        
    print "$db\n";
    print "Samples with reqs $samples_with_reqs , $samples_with_reqs_lm\n";
    print "Libraries made $libs , $libs_lm \n";
    print "BGI samples $bgi\n";
    print "Samples with seq $samples_with_seq , $samples_with_seq_lm\n";
    print "Finished samples $finished , $finished_lm \n";
    print "\n";
}

my ($sum, $count);
foreach (@wgs_yield){
    $sum += $_->[1];
    $count++;
}

printf ("mean wgs yield $month/$year = %.1f GBp/lane \n", $sum/$count) if $count;

$sum = $count = 0;
foreach (@exome_yield){
    $sum += $_->[1];
    $count++;
}

printf ("mean exome yield $month/$year = %.1f GBp/lane \n", $sum/$count) if $count;


sub get_lane_yield_last_month {
    my ($dbh,$year,$month) = @_;
    my $sql = qq[
        select left(name,6) as lanename, sum(raw_bases)/1e9 as gbp
        from latest_lane l
        where name like '%#%' 
        and year(l.run_date)= $year 
        and month(l.run_date)= $month
        group by lanename;
    ];

    my $sth = $dbh->prepare($sql);
    my $rows = []; 
    if ($sth->execute()) {
        $rows = $sth->fetchall_arrayref();
    }
    else {
	print_and_exit(sprintf('Cannot retrieve yield: %s', $DBI::errstr));
    }
    return $rows;
}

sub get_finished_last_month {
    my ($dbh,$gbp,$year,$month) = @_;

    my $sql = qq[
     select 
        (select count(*) from 
            ( select s.name, round(sum(m.rmdup_bases_mapped)/1e9) as mapped 
                from latest_sample s, latest_library  lib, latest_lane l, latest_mapstats m 
                where s.sample_id = lib.sample_id 
                and lib.library_id = l.library_id 
                and lib.seq_centre_id =1 
                and l.lane_id = m.lane_id 
                and gt_status="confirmed" 
                group by s.name having mapped >= ? 
            ) tot 
        ) 
      - 
        (select count(*) from 
            ( select s.name, round(sum(m.rmdup_bases_mapped)/1e9) as mapped 
                from latest_sample s, latest_library  lib,latest_lane l, latest_mapstats m 
                where s.sample_id = lib.sample_id 
                and lib.library_id = l.library_id 
                and lib.seq_centre_id =1 
                and l.lane_id = m.lane_id 
                and gt_status="confirmed" 
                and l.run_date < '$year-$month-01'  
                group by s.name having mapped >= ? 
            ) last_month) 
     as finished_last_month;
            ];

    my $sth = $dbh->prepare($sql);
    my $count; 
    if ($sth->execute($gbp,$gbp)) {
        ($count) = $sth->fetchrow_array;
    }
    else {
	print_and_exit(sprintf('Cannot retrieve finished: %s', $DBI::errstr));
    }
    return $count;
}


sub get_finished{
    my ($dbh,$gbp) = @_;

    my $sql = qq[
        select count(*) from (
            select s.name, round(sum(m.rmdup_bases_mapped)/1e9) as mapped 
            from latest_sample s, latest_library  lib, latest_lane l, latest_mapstats m,seq_centre sc
            where s.sample_id = lib.sample_id 
            and lib.library_id = l.library_id 
            and lib.seq_centre_id = sc.seq_centre_id 
            and sc.name="SC"
            and l.lane_id = m.lane_id 
            and gt_status="confirmed" 
            group by s.name having mapped >= ?) z
            ];

    my $sth = $dbh->prepare($sql);
    my $count; 
    if ($sth->execute($gbp)) {
        ($count) = $sth->fetchrow_array;
    }
    else {
	print_and_exit(sprintf('Cannot retrieve finished: %s', $DBI::errstr));
    }
    return $count;
}



sub get_samples_with_seq_last_month{
    my ($dbh,$year,$month) = @_;

    my $sql = qq[
        select count(distinct(s.name)) 
        from latest_sample s, latest_library lib, latest_lane l , seq_centre sc 
        where s.sample_id = lib.sample_id 
        and lib.seq_centre_id = sc.seq_centre_id 
        and sc.name="SC" 
        and lib.library_id = l.library_id
        and year(l.run_date)= $year
        and month(l.run_date)=$month
            ];

    my $sth = $dbh->prepare($sql);
    my $count; 
    if ($sth->execute()) {
        ($count) = $sth->fetchrow_array;
    }
    else {
	print_and_exit(sprintf('Cannot retrieve bgi: %s', $DBI::errstr));
    }
    return $count;
}

sub get_samples_with_seq {
    my ($dbh) = @_;

    my $sql = qq[
        select count(distinct(s.name)) 
        from latest_sample s, latest_library lib, latest_lane l , seq_centre sc 
        where s.sample_id = lib.sample_id 
        and lib.seq_centre_id = sc.seq_centre_id 
        and sc.name="SC" 
        and lib.library_id = l.library_id
            ];

    my $sth = $dbh->prepare($sql);
    my $count; 
    if ($sth->execute()) {
        ($count) = $sth->fetchrow_array;
    }
    else {
	print_and_exit(sprintf('Cannot retrieve bgi: %s', $DBI::errstr));
    }
    return $count;
}


sub get_bgi{
    my ($dbh) = @_;

    my $sql = qq[
        select count(distinct(s.name)) 
        from latest_sample s, latest_library lib, latest_lane l , seq_centre sc 
        where s.sample_id = lib.sample_id 
        and lib.seq_centre_id = sc.seq_centre_id 
        and sc.name="BGI" 
        and lib.library_id = l.library_id
            ];

    my $sth = $dbh->prepare($sql);
    my $count; 
    if ($sth->execute()) {
        ($count) = $sth->fetchrow_array;
    }
    else {
	print_and_exit(sprintf('Cannot retrieve bgi: %s', $DBI::errstr));
    }
    return $count;
}



sub get_libraries_last_month{
    my ($dbh,$year,$month) = @_;

    my $sql = qq[
        select count(*) from 
            (select library_id, min(changed)  as orig_date 
            from library lib, seq_centre sc 
            where lib.seq_centre_id = sc.seq_centre_id
            and sc.name="SC"
            group by library_id 
            having year(orig_date) = $year
            and month(orig_date)=$month
            ) z;
        
            ];

    my $sth = $dbh->prepare($sql);
    my $count; 
    if ($sth->execute()) {
        ($count) = $sth->fetchrow_array;
    }
    else {
	print_and_exit(sprintf('Cannot retrieve libraries: %s', $DBI::errstr));
    }
    return $count;
}

sub get_libraries{

    my ($dbh) = @_;

    my $sql = qq[
        select  count(*)
        from latest_library lib, seq_centre sc 
        where lib.seq_centre_id = sc.seq_centre_id
        and sc.name="SC"
        ];

    my $sth = $dbh->prepare($sql);
    my $count; 
    if ($sth->execute()) {
        ($count) = $sth->fetchrow_array;
    }
    else {
	print_and_exit(sprintf('Cannot retrieve libraries: %s', $DBI::errstr));
    }
    return $count;
}


sub get_samples_with_requests{
    my ($dbh) = @_;

    my $sql = qq[
        select  count(distinct(s.name)) 
        from latest_sample s, latest_library_request lr 
        where s.sample_id = lr.sample_id;
            ];

    my $sth = $dbh->prepare($sql);
    my $count; 
    if ($sth->execute()) {
        ($count) = $sth->fetchrow_array;
    }
    else {
	print_and_exit(sprintf('Cannot retrieve samples with requests: %s', $DBI::errstr));
    }
    return $count;
}


sub get_samples_with_requests_last_month{
    my ($dbh,$year,$month) = @_;

    my $sql = qq[
        select count(*) from 
            (select library_request_id, min(changed) as orig_date 
            from library_request 
            group by library_request_id 
            having year(orig_date) = $year
            and month(orig_date)= $month
            ) z;
            ];

    my $sth = $dbh->prepare($sql);
    my $count; 
    if ($sth->execute()) {
        ($count) = $sth->fetchrow_array;
    }
    else {
	print_and_exit(sprintf('Cannot retrieve samples with requests: %s', $DBI::errstr));
    }
    return $count;
}
