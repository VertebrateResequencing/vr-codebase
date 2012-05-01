#!/usr/bin/perl -w
use strict;
use warnings;
use File::Spec;

BEGIN {
    use Test::Most;
    eval {
        require VRTrack::Testconfig;
    };
    if ($@) {
        plan skip_all => "Skipping all tests because VRTrack tests have not been configured";
    }
    else {
        plan tests => 10;
    }

    use_ok('VRTrack::Factory');
}

my %cd = ( database => VRTrack::Testconfig->config('test_db'),
           host     => VRTrack::Testconfig->config('host'),
           port     => VRTrack::Testconfig->config('port'),
           user     => VRTrack::Testconfig->config('user'),
           password => VRTrack::Testconfig->config('password') );

# populate db with full real data to test with
open(my $mysqlfh, "| mysql -h$cd{host} -u$cd{user} -p$cd{password} -P$cd{port}") || die "could not connect to VRTrack database for testing\n";
print $mysqlfh "drop database if exists $cd{database};\n";
print $mysqlfh "create database $cd{database};\n";
print $mysqlfh "use $cd{database};\n";
open(my $sql_fh, File::Spec->catfile(qw(t data vrtrack_pombe_wgs.sql))) || die "Could not open sql file\n";
while (<$sql_fh>) {
    print $mysqlfh $_;
}
close($mysqlfh);
close($sql_fh);

# see if we can use the Factory to connect, instead of new($connection_details)
$ENV{VRTRACK_HOST} = $cd{host};
$ENV{VRTRACK_PORT} = $cd{port};
$ENV{VRTRACK_RW_USER} = $cd{user};
$ENV{VRTRACK_PASSWORD} = $cd{password};
ok my $vrtrack = VRTrack::Factory->instantiate(database => $cd{database}, mode => 'rw'), 'VRTrack::Factory->instantiate worked';
isa_ok $vrtrack, 'VRTrack::VRTrack';

# test VRTrack::Factory->databases
my %dbs = map { $_ => 1 } VRTrack::Factory->databases(1);
is exists $dbs{$cd{database}}, 1, 'VRTrack::Factory->databases returned at least our test database';

# get_lanes
my @lanes = grep { $_->isa('VRTrack::Lane') } $vrtrack->get_lanes();
is scalar(@lanes), 285, 'get_lanes() gave all lanes';
@lanes = grep { $_->isa('VRTrack::Lane') } $vrtrack->get_lanes(project_regex => 'diversity');
is scalar(@lanes), 285, 'get_lanes(project_regex) worked';
@lanes = map { $_->name } grep { $_->isa('VRTrack::Lane') } $vrtrack->get_lanes(sample_regex => '1996STDY5244898');
is_deeply \@lanes, ['7816_3#1', '7211_8#1', '7413_5#1'], 'get_lanes(sample_regex) worked';
@lanes = map { $_->name } grep { $_->isa('VRTrack::Lane') } $vrtrack->get_lanes(library_regex => '4858163|4074489');
is_deeply \@lanes, ['7816_3#1', '7211_8#1', '7413_5#1'], 'get_lanes(library_regex) worked';

# lane_info
my %i = $vrtrack->lane_info('7816_3#1');
is_deeply [$i{centre},
           $i{project},
           $i{study},
           $i{species},
           $i{population},
           $i{individual},
           $i{sample},
           $i{seq_tech},
           $i{library},
           $i{lane},
           $i{withdrawn},
           $i{insert_size},
           $i{individual_coverage}],
          ['SC',
           'S.pombe genetic diversity',
           'ERP000979',
           'Schizosaccharomyces pombe',
           'Population',
           'JB4',
           '1996STDY5244898',
           'SLX',
           '4858163',
           '7816_3#1',
           undef,
           197,
           undef], 'lane_info worked';
%i = $vrtrack->lane_info('7816_3#1', get_coverage => 1, genome_size => 12631379);
is $i{individual_coverage}, 85.95, 'lane_info with get_coverage worked';

exit;
