#!/usr/bin/perl -w
use strict;
use warnings;

BEGIN {
    use Test::Most;
    eval {
        require VRTrack::Testconfig;
    };
    if ($@) {
        plan skip_all => "Skipping all tests because VRTrack tests have not been configured";
    }
    else {
         plan tests => 227;
    }
    
    use_ok('VRTrack::VRTrack');
    use_ok('VRTrack::Library_request');
    use_ok('VRTrack::Library');
    use_ok('VRTrack::Seq_request');
    use_ok('VRTrack::Lane');
    use_ok('VRTrack::File');
    use_ok('VRTrack::Mapstats');
}

my $connection_details = { database => VRTrack::Testconfig->config('test_db'),
                           host     => VRTrack::Testconfig->config('host'),
                           port     => VRTrack::Testconfig->config('port'),
                           user     => VRTrack::Testconfig->config('user'),
                           password => VRTrack::Testconfig->config('password'), 
                        };

# access the db as normal with new()
ok my $vrtrack = VRTrack::VRTrack->new($connection_details), 'VRTrack new() returned something';
isa_ok $vrtrack, 'VRTrack::VRTrack';

# generic testing of all the core objects
my @core_objects = qw(Project Sample Library_request Library Seq_request Lane File Mapstats);
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
    if ("VRTrack::$class"->can('name')) {
        $vrobj = "VRTrack::$class"->create($vrtrack, $name);
    }
    else {
        $vrobj = "VRTrack::$class"->create($vrtrack);
    }


    if ($class eq 'Request') {
        #*** Request objects can't be created? How the fuck do I test them?
        ok ! $vrobj, 'VRTrack::Request->create returns undef since they can\'t be created anew';
        next;
    }
    ok $vrobj, $class.'->create returned something';
    isa_ok $vrobj, "VRTrack::$class";
    isa_ok $vrobj, "VRTrack::Core_obj";
    isa_ok $vrobj, "VRTrack::Table_obj";
    
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
   
    
    my $values = "VRTrack::$class"->_all_values_by_field($vrtrack,'row_id');
    is($values->[0], 1, " _all_values_by_field should return first row_id value as 1");
     
    # shared by some, but not all
    if ($vrobj->can('hierarchy_name')) {
        isa_ok $vrobj, "VRTrack::Hierarchy_obj";
        
        is $vrobj->hierarchy_name, $name, "hierarchy_name set to name by default $vrobj";
        is $vrobj->hierarchy_name($new_name), $new_name, 'hierarchy_name could be reset, but it was not made file-system safe';
    }
    if ($vrobj->can('ssid')) {
        is $vrobj->ssid, undef, 'ssid is undef by default';
        $ssid++;
        is $vrobj->ssid($ssid), $ssid, 'ssid could be set';
    }
    if ($vrobj->can('name')) {
        is $vrobj->name, $name, "name was set to that requested in the create call $vrobj";
        is $vrobj->name($new_name), $new_name, "name could be reset $vrobj";
    }
    if ($vrobj->can('new_by_name')) {
        # (not all objects that can name() can new_by_name())
        isa_ok $vrobj, "VRTrack::Named_obj";
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
    }
    else {
        $vrobj = "VRTrack::$class"->new($vrtrack, $vrobj->id);
    }
    ok $vrobj, $class.'->new_by_name worked on the new name we set';
    is $vrobj->row_id, 2, 'row_id was actually just autoincremented to 2, not set to 99 like we wanted';
    is $vrobj->id, 999, 'id was actually reset';
    is $vrobj->note_id, 1, 'note_id was actually reset';
    isnt $vrobj->changed, $orig_changed, 'changed was automatically updated';
    is $vrobj->is_latest, 1, 'is_latest is still true';
    if ($vrobj->can('name')) {
        is $vrobj->name, $new_name, 'name was actually reset';
    }
    if ($vrobj->can('hierarchy_name')) {
        is $vrobj->hierarchy_name, $new_name, 'hierarchy_name was actually reset';
    }
    if ($vrobj->can('ssid')) {
        is $vrobj->ssid($ssid), $ssid, 'ssid was actually reset';
    }
}
# cleanup
my $dbh = $vrtrack->{_dbh};
foreach ($dbh->tables()){
    next if /`latest_/;
    next if /schema_version/;
    $dbh->do("TRUNCATE TABLE $_");
}
exit;
