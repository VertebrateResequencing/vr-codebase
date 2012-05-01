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
        plan tests => 4;
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




exit;
