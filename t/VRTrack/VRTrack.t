#!/usr/bin/perl -w
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 27;

    use_ok('VRTrack::VRTrack');
}

my $connection_details = { database => 'vrtrack_test',
                           host     => 'mcs4a',
                           port     => 3306,
                           user     => 'vreseq_rw',
                           password => 't3aml3ss' };

# drop, create and then populate our testing database with the VRTrack schema
# *** need a better way of doing this...
ok my @schema = VRTrack::VRTrack->schema(), 'schema() returned something';
open(my $mysqlfh, "| mysql -h$connection_details->{host} -u$connection_details->{user} -p$connection_details->{password}") || die "could not connect to database for testing\n";
print $mysqlfh "drop database if exists $connection_details->{database};\n";
print $mysqlfh "create database $connection_details->{database};\n";
print $mysqlfh "use $connection_details->{database};\n";
foreach my $sql (@schema) {
    print $mysqlfh $sql;
}
close($mysqlfh);

# access the db as normal with new()
ok my $vrtrack = VRTrack::VRTrack->new($connection_details), 'new() returned something';
isa_ok $vrtrack, 'VRTrack::VRTrack';

# project
my $name = 'Project1';
my $vrproj = VRTrack::Project->new_by_name($vrtrack, $name);
ok ! $vrproj, 'Project->new_by_name returns undef on empty db';
ok $vrproj = VRTrack::Project->create($vrtrack, $name), 'Project->create returned something';
isa_ok $vrproj, 'VRTrack::Project';
is_deeply [$vrproj->id, $vrproj->row_id, $vrproj->ssid], [1, 1, undef], 'first object in the db gets id 1 and row_id 1, and ssid is not set yet';
ok $vrproj = VRTrack::Project->new_by_name($vrtrack, $name), 'Project->new_by_name now works';
is $vrproj->id, 1, 'and it returns the same project object as we got during creation';
is $vrproj->name, $name, 'and it has the correct name';
$name = 'Project2';
ok $vrproj = VRTrack::Project->create($vrtrack, $name), 'Project->create worked again';
is_deeply [$vrproj->id, $vrproj->row_id], [2, 2], 'ids were incremented for the new object';
my $ssid = 1234;
is $vrproj->ssid($ssid), $ssid, 'ssid could be set';
ok $vrproj->update, 'and written to the db';
$vrproj = VRTrack::Project->new_by_name($vrtrack, $name);
is_deeply [$vrproj->id, $vrproj->row_id, $vrproj->ssid], [2, 3, $ssid], 'getting it back with new_by_name shows the id has not changed, the row_id has incremented and the ssid could be retrieved again';

# history test - can we get row_id 2 again?...

my $study_id = 'test_study';
is $vrproj->study_id($study_id), $study_id, 'study_id could be set';
ok $vrproj->update(), 'and written to the db';
#*** should test all other attributes...



#*** could probably have a generic tester of all core obj types...

# lane
$name = '1111_1';
my $vrlane = VRTrack::Lane->new_by_name($vrtrack, $name);
ok ! $vrlane, 'Lane->new_by_name returns undef before we\'ve added any lanes';
ok $vrlane = VRTrack::Lane->create($vrtrack, $name), 'Lane->create returned something';
isa_ok $vrlane, 'VRTrack::Lane';

# lane children
$name = "${name}_1.fastq";
my $vrfile = $vrlane->get_file_by_name($name);
ok ! $vrfile, 'Lane->get_file_by_name returns undef before we\'ve added any files';
ok $vrfile = $vrlane->add_file($name), 'Lane->add_file did something';

ok my $mapping = $vrlane->add_mapping(), 'Lane->add_mapping did something';
is $mapping->raw_reads(999), 999, 'Mapstats->raw_reads could be set';
is $mapping->raw_bases(111), 111, 'Mapstats->raw_bases could be set';
ok $mapping->update, 'changes could be written to db';

exit;
