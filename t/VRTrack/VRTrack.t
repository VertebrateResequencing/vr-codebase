#!/usr/bin/perl -w
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 88;

    use_ok('VRTrack::VRTrack');
    use_ok('VRTrack::Request');
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

# generic testing of all the core objects
my @core_objects = qw(Project Sample Library Request Lane File Mapstats);
foreach my $class (@core_objects) {
    my $name = $class.'_test';
    my $vrobj;
    if ("VRTrack::$class"->can('new_by_name')) {
        $vrobj = "VRTrack::$class"->new_by_name($vrtrack, $name);
        ok ! $vrobj, $class.'->new_by_name returns undef before we\'ve added any';
    }
    
    $vrobj = "VRTrack::$class"->create($vrtrack, $class ne 'Mapstats' ? $name : undef);
    if ($class eq 'Request') {
        #*** Request objects can't be created? How the fuck do I test them?
        ok ! $vrobj, 'VRTrack::Request->create returns undef';
        next;
    }
    ok $vrobj, $class.'->create returned something';
    isa_ok $vrobj, "VRTrack::$class";
    isa_ok $vrobj, "VRTrack::Core_obj";
    
    # generic methods
    unless ($class eq 'Project') {
        is $vrobj->row_id, 1, 'being the first object of it\'s kind, row_id is 1';
        is $vrobj->id, 1, 'first of an object also gets an id of 1';
    }
    is $vrobj->note_id, undef, 'note_id is undef by default';
    like ($vrobj->changed, qr/^\d{4}\-\d{2}\-\d{2} \d{2}:\d{2}:\d{2}$/, 'changed is set');
    is $vrobj->is_latest, 1, 'is_latest is true';
    
    # shared by some, but not all
    if ($vrobj->can('hierarchy_name')) {
        #isa_ok $vrobj, "VRTrack::Hierarchy_obj";
        
        # Sample hierarchy_name comes from Individual, which we haven't made yet
        unless ($class eq 'Sample') {
            is $vrobj->hierarchy_name, $name, 'hierarchy_name set to name by default';
        }
    }
    if ($vrobj->can('ssid')) {
        is $vrobj->ssid, undef, 'ssid is undef by default';
    }
    if ($vrobj->can('name')) {
        is $vrobj->name, $name, 'name was set';
    }
}

# test bits of the core classes that are unique to themselves
# lane children
$name = 'Lane_test';
my $vrlane = VRTrack::Lane->new_by_name($vrtrack, $name);
$name = "${name}_1.fastq";
my $vrfile = $vrlane->get_file_by_name($name);
ok ! $vrfile, 'Lane->get_file_by_name returns undef before we\'ve added any files';
ok $vrfile = $vrlane->add_file($name), 'Lane->add_file did something';

ok my $mapping = $vrlane->add_mapping(), 'Lane->add_mapping did something';
is $mapping->raw_reads(999), 999, 'Mapstats->raw_reads could be set';
is $mapping->raw_bases(111), 111, 'Mapstats->raw_bases could be set';
ok $mapping->update, 'changes could be written to db';

exit;
