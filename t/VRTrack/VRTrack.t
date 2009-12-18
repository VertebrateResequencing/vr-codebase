#!/usr/bin/perl -w
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 151;

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

# generic testing of all the core objects
my @core_objects = qw(Project Sample Library Request Lane File Mapstats);
my $ssid = 1233;
foreach my $class (@core_objects) {
    my $name = $class.'_test';
    my $new_name = $name.'|new';
    my $new_hname = $name.'_new';
    
    my $vrobj;
    if ("VRTrack::$class"->can('new_by_name')) {
        $vrobj = "VRTrack::$class"->new_by_name($vrtrack, $name);
        ok ! $vrobj, $class.'->new_by_name returns undef before we\'ve added any';
    }
    
    $vrobj = "VRTrack::$class"->create($vrtrack, $class ne 'Mapstats' ? $name : undef);
    if ($class eq 'Request') {
        #*** Request objects can't be created? How the fuck do I test them?
        ok ! $vrobj, 'VRTrack::Request->create returns undef since they can\'t be created anew';
        next;
    }
    ok $vrobj, $class.'->create returned something';
    isa_ok $vrobj, "VRTrack::$class";
    isa_ok $vrobj, "VRTrack::Core_obj";
    
    # generic methods
    is $vrobj->row_id, 1, 'being the first object of it\'s kind, row_id is 1';
    is $vrobj->row_id(99), 99, 'row_id seems as if it can be reset';
    is $vrobj->id, 1, 'first of an object also gets an id of 1';
    is $vrobj->id(999), 999, 'id can be reset';
    is $vrobj->note_id, undef, 'note_id is undef by default';
    is $vrobj->note_id(1), 1, 'note_id could be set';
    my $orig_changed = $vrobj->changed;
    like ($orig_changed, qr/^\d{4}\-\d{2}\-\d{2} \d{2}:\d{2}:\d{2}$/, 'changed is set');
    is $vrobj->is_latest, 1, 'is_latest is true';
    
    # shared by some, but not all
    if ($vrobj->can('hierarchy_name')) {
        #isa_ok $vrobj, "VRTrack::Hierarchy_obj";
        
        # Sample hierarchy_name comes from Individual, which we haven't made yet
        unless ($class eq 'Sample') {
            is $vrobj->hierarchy_name, $name, 'hierarchy_name set to name by default';
            is $vrobj->hierarchy_name($new_name), $new_name, 'hierarchy_name could be reset, but it was not made file-system safe';
        }
    }
    if ($vrobj->can('ssid')) {
        is $vrobj->ssid, undef, 'ssid is undef by default';
        $ssid++;
        is $vrobj->ssid($ssid), $ssid, 'ssid could be set';
    }
    if ($vrobj->can('name')) {
        is $vrobj->name, $name, 'name was set to that requested in the create call';
        is $vrobj->name($new_name), $new_name, 'name could be reset';
    }
    
    # before we can update(), certain things need to be set on certain objects
    if ($class eq 'Library') {
        $vrobj->sample_id(1);
    }
    elsif ($class eq 'Lane') {
        $vrobj->library_id(1);
    }
    elsif ($class eq 'File' || $class eq 'Mapstats') {
        $vrobj->lane_id(1);
    }
    
    # now test if the things we set really got set
    sleep(1); # so that changed() will be changed
    ok $vrobj->update, 'update worked to write changes to db';
    if ("VRTrack::$class"->can('new_by_name')) {
        $vrobj = "VRTrack::$class"->new_by_name($vrtrack, $new_name);
        ok $vrobj, $class.'->new_by_name worked on the new name we set';
        is $vrobj->row_id, 2, 'row_id was actually just autoincremented to 2, not set to 99 like we wanted';
        is $vrobj->id, 999, 'id was actually reset';
        is $vrobj->note_id, 1, 'note_id was actually reset';
        isnt $vrobj->changed, $orig_changed, 'changed was automatically updated';
        is $vrobj->is_latest, 1, 'is_latest is still true';
        is $vrobj->name, $new_name, 'name was actually reset';
        
        if ($vrobj->can('hierarchy_name')) {
            unless ($class eq 'Sample') {
                is $vrobj->hierarchy_name, $new_name, 'hierarchy_name was actually reset';
            }
        }
        if ($vrobj->can('ssid')) {
            is $vrobj->ssid($ssid), $ssid, 'ssid was actually reset';
        }
    }
}

# core objects keep a history
my $name = 'Project_test2';
ok my $vrproj = VRTrack::Project->create($vrtrack, $name), 'Project->create worked again';
is_deeply [$vrproj->id, $vrproj->row_id], [3, 3], 'ids were incremented for the new object';
$vrproj->ssid(123);
ok $vrproj->update(), 'after changing an attribute was able to update';
$vrproj = VRTrack::Project->new_by_name($vrtrack, $name);
is_deeply [$vrproj->id, $vrproj->row_id], [3, 4], 'getting it back with new_by_name shows the id has not changed and the row_id has incremented';
# history test - can we get row_id 2 again?...

# test bits of the core classes that are unique to themselves
{
    # project
    my $study_id = 'test_study';
    is $vrproj->study_id($study_id), $study_id, 'study_id could be set';
    
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
}

# tests for non-core classes
#...

exit;
