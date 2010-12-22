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
        plan tests => 259;
    }

    use_ok('VRTrack::VRTrack');
    use_ok('VRTrack::History');
    use_ok('VRTrack::Library');
    use_ok('VRTrack::Lane');
}

my $connection_details = { database => VRTrack::Testconfig->config('test_db'),
                           host     => VRTrack::Testconfig->config('host'),
                           port     => VRTrack::Testconfig->config('port'),
                           user     => VRTrack::Testconfig->config('user'),
                           password => VRTrack::Testconfig->config('password'), 
                        };

# access the db as normal with new()
ok my $vrtrack = VRTrack::VRTrack->new($connection_details), 'new() returned something';
isa_ok $vrtrack, 'VRTrack::VRTrack';

# generic testing of all the core objects
my @core_objects = qw(Project Sample Library Request Lane File Mapstats);
my $ssid = 1233;

# core objects keep a history
my $name = 'Project_historytest';
ok my $vrproj = VRTrack::Project->create($vrtrack, $name), 'Project->create worked';
my $vrproj_id = $vrproj->id;
my $orig_row_id = $vrproj->row_id;
my $latest_row_id = $orig_row_id +1;
my $changed = $vrproj->changed;
sleep(1);
$vrproj->ssid(123);
ok $vrproj->update(), 'after changing an attribute was able to update';
$vrproj = VRTrack::Project->new_by_name($vrtrack, $name);
my $changed2 = $vrproj->changed;
is_deeply [$vrproj->id, $vrproj->row_id], [$vrproj_id, $latest_row_id], 'getting it back with new_by_name shows the id has not changed whilst the row_id has incremented';
$vrproj = VRTrack::Project->new($vrtrack, $vrproj_id, $changed2);
is_deeply [$vrproj->id, $vrproj->row_id, $vrproj->changed], [$vrproj_id, $orig_row_id, $changed], 'new(... new date) gives us the old version of the object';

$vrproj = VRTrack::Project->new($vrtrack, $vrproj_id, 'latest');
is_deeply [$vrproj->id, $vrproj->row_id, $vrproj->changed], [$vrproj_id, $latest_row_id, $changed2], 'new(... latest) gives us the latest version of the object';
is $vrproj->global_history_date($changed2), $changed2, 'global_history_date could be set';
$vrproj = VRTrack::Project->new($vrtrack, $vrproj_id);
is_deeply [$vrproj->id, $vrproj->row_id, $vrproj->changed], [$vrproj_id, $orig_row_id, $changed], 'new() gives us the old version of the object when global was set';
$vrproj = VRTrack::Project->new($vrtrack, $vrproj_id, 'latest');
is_deeply [$vrproj->id, $vrproj->row_id, $vrproj->changed], [$vrproj_id, $latest_row_id, $changed2], 'new(... latest) gives us the latset version inspite of global setting';
is $vrproj->global_history_date('latest'), 'latest', 'global_history_date could be set back to latest';
$vrproj = VRTrack::Project->new_by_name($vrtrack, $name, $orig_row_id);
is_deeply [$vrproj->id, $vrproj->row_id, $vrproj->changed], [$vrproj_id, $orig_row_id, $changed], 'new_by_name(... old row_id) gives us the old version of the object';
$vrproj = VRTrack::Project->new_by_name($vrtrack, $name);
is_deeply [$vrproj->id, $vrproj->row_id, $vrproj->changed], [$vrproj_id, $latest_row_id, $changed2], 'new_by_name() still gives us the latest version of the object';
cmp_bag ([$vrproj->row_ids], [$orig_row_id, $latest_row_id], 'row_ids worked');
$vrproj = VRTrack::Project->new_by_name($vrtrack, $name, $orig_row_id);    # get old object again
$vrproj->ssid(321);
is $vrproj->ssid, 321, 'can change attributes on historical objects'; # this is mostly useless, except for is_latest(1), of course
ok ! $vrproj->update(), 'can\'t update on historical objects after changing attribute';

# more in-depth history testing using History object. First create some lanes
# and libraries and see if we can see the state before and after a swap
my $liba = VRTrack::Library->create($vrtrack, 'history_lib_a');
$liba->add_lane('history_lane_a');
$liba->add_lane('history_lane_b');
$liba->update;
my $libb = VRTrack::Library->create($vrtrack, 'history_lib_b');
$libb->add_lane('history_lane_c');
$libb->add_lane('history_lane_d');
$libb->update;
my %initial_lane_changed;
foreach ('a'..'d') {
    my $lane_name = "history_lane_$_";
    my $lane = VRTrack::Lane->new_by_name($vrtrack, $lane_name);
    $initial_lane_changed{$lane_name} = $lane->changed;
    $lane->is_processed('mapped', 1);
    $lane->update;
}
sleep(1);
my $lane = VRTrack::Lane->new_by_name($vrtrack, 'history_lane_b');
$lane->is_processed('swapped', 1);
$lane->library_id($libb->id);
$lane->update;
$lane = VRTrack::Lane->new_by_name($vrtrack, 'history_lane_d');
$lane->is_processed('swapped', 1);
$lane->library_id($liba->id);
$lane->update;
$liba->name('history_lib_a_changed');
$liba->update;
# now start the tests
$lane = VRTrack::Lane->new_by_name($vrtrack, 'history_lane_a');
is_deeply [$lane->library_id, VRTrack::Library->new($vrtrack, $lane->library_id)->name, $lane->is_processed('swapped')], [$liba->id, 'history_lib_a_changed', 0], 'after swap lane_a state was correct';
$lane = VRTrack::Lane->new_by_name($vrtrack, 'history_lane_b');
is_deeply [$lane->library_id, VRTrack::Library->new($vrtrack, $lane->library_id)->name, $lane->is_processed('swapped')], [$libb->id, 'history_lib_b', 1], 'after swap lane_b state was correct';
$lane = VRTrack::Lane->new_by_name($vrtrack, 'history_lane_c');
is_deeply [$lane->library_id, VRTrack::Library->new($vrtrack, $lane->library_id)->name, $lane->is_processed('swapped')], [$libb->id, 'history_lib_b', 0], 'after swap lane_c state was correct';
$lane = VRTrack::Lane->new_by_name($vrtrack, 'history_lane_d');
is_deeply [$lane->library_id, VRTrack::Library->new($vrtrack, $lane->library_id)->name, $lane->is_processed('swapped')], [$liba->id, 'history_lib_a_changed', 1], 'after swap lane_d state was correct';
# do some time travel
my $hist = VRTrack::History->new();
$lane = VRTrack::Lane->new_by_name($vrtrack, 'history_lane_a');
my $datetime = $hist->state_change($lane, 'library_id');
is $datetime, $initial_lane_changed{'history_lane_a'}, 'lane_a never changed its library';
$datetime = $hist->was_processed($lane, 'swapped');
is $datetime, 'latest', 'lane_a was never swapped';
$lane = VRTrack::Lane->new_by_name($vrtrack, 'history_lane_c');
$datetime = $hist->state_change($lane, 'library_id');
is $datetime, $initial_lane_changed{'history_lane_c'}, 'lane_c never changed its library';
$datetime = $hist->was_processed($lane, 'swapped');
is $datetime, 'latest', 'lane_c was never swapped';
$lane = VRTrack::Lane->new_by_name($vrtrack, 'history_lane_b');
$datetime = $hist->state_change($lane, 'library_id');
isnt $datetime, 'latest', 'lane_b changed its library';
$datetime = $hist->was_processed($lane, 'swapped');
isnt $datetime, 'latest', 'lane_b was swapped';
$hist->time_travel($datetime);
$lane = VRTrack::Lane->new_by_name($vrtrack, 'history_lane_b');
is_deeply [$lane->library_id, VRTrack::Library->new($vrtrack, $lane->library_id)->name, $lane->is_processed('swapped')], [$liba->id, 'history_lib_a', 0], 'before swap lane_b state was correct';
$hist->time_travel('latest');
$lane = VRTrack::Lane->new_by_name($vrtrack, 'history_lane_d');
$datetime = $hist->state_change($lane, 'library_id');
isnt $datetime, 'latest', 'lane_d changed its library';
$datetime = $hist->was_processed($lane, 'swapped');
isnt $datetime, 'latest', 'lane_d was swapped';
$hist->time_travel($datetime);
$lane = VRTrack::Lane->new_by_name($vrtrack, 'history_lane_d');
is_deeply [$lane->library_id, VRTrack::Library->new($vrtrack, $lane->library_id)->name, $lane->is_processed('swapped')], [$libb->id, 'history_lib_b', 0], 'before swap lane_d state was correct';
$hist->time_travel('latest');
$lane = VRTrack::Lane->new_by_name($vrtrack, 'history_lane_d');
is_deeply [$lane->library_id, VRTrack::Library->new($vrtrack, $lane->library_id)->name, $lane->is_processed('swapped')], [$liba->id, 'history_lib_a_changed', 1], 'after swap lane_d state still correct';

# Testing unsetting and resetting is_latest
# unset the latest flag to mark an object as 'not current' - e.g. _s_ files before splitting
$name = 'Project_islatest_test';

ok $vrproj = VRTrack::Project->create($vrtrack, $name), 'Project->create worked again';
my $proj_id = $vrproj->id;
sleep(1);
$vrproj->ssid(555);
ok $vrproj->update(), 'after changing an attribute was able to update';
$vrproj = VRTrack::Project->new_by_name($vrtrack, $name); # get latest version again for changed 
$changed = $vrproj->changed;
is $vrproj->is_latest(),1,'latest object has is_latest true';
ok $vrproj = VRTrack::Project->new($vrtrack, $proj_id, $changed), 'can retrieve old object';
ok ! $vrproj->is_latest(),'historical object has is_latest false';
is_deeply [$vrproj->is_latest(0),$vrproj->{'unset_latest'}],[0,undef], 'can\'t unset latest on a historical object';
$vrproj = VRTrack::Project->new_by_name($vrtrack, $name); # latest obj
is_deeply [$vrproj->is_latest(0),$vrproj->{'unset_latest'}],[1,1], 'can unset latest on the latest version of object';
ok $vrproj->update(), 'can update an object after unsetting latest';
$vrproj = VRTrack::Project->new_by_name($vrtrack, $name);
ok ! $vrproj, 'can\'t retrieve object by name after unsetting latest';
$vrproj = VRTrack::Project->new($vrtrack, $proj_id,'latest');
ok ! $vrproj, 'can\'t retrieve latest version of object after unsetting latest';
$vrproj = VRTrack::Project->new($vrtrack, $proj_id,'most_recent');
ok $vrproj, 'can retrieve most_recent version of object after unsetting latest';
ok ! $vrproj->is_latest(),'most_recent object can have is_latest false';
ok $vrproj->is_latest(1), 'can set is_latest on non-latest object';
ok $vrproj->update(), 'can update after resetting latest';
$vrproj = VRTrack::Project->new($vrtrack, $proj_id,'latest');
ok $vrproj, 'can retrieve latest version of object after resetting latest';

# let's do a more in-depth overall test of history, to see if we can pull back
# a whole hierarchy as it was before adding new stuff to the hierarchy:
{
    undef($vrtrack);
    $hist->time_travel('latest');
    new_test_db();
    $vrtrack = VRTrack::VRTrack->new($connection_details);
    
    # setup the state we want to test being able to get back to:
    my %hierarchy;
    build_hierarchy($vrtrack, 2, \%hierarchy);
    
    sub build_hierarchy {
        my ($vrtrack, $objs_per_level, $hash) = @_;
        foreach my $p (1..$objs_per_level) {
            my $p_name = "p$p";
            my $project = $vrtrack->get_project_by_name($p_name) || $vrtrack->add_project($p_name);
            
            foreach my $s (1..$objs_per_level) {
                my $s_name = "$p_name.s$s";
                my $sample = $project->get_sample_by_name($s_name) || $project->add_sample($s_name);
                
                foreach my $l (1..$objs_per_level) {
                    my $l_name = "$s_name.l$l";
                    my $library = $sample->get_library_by_name($l_name) || $sample->add_library($l_name);
                    
                    foreach my $lane (1..$objs_per_level) {
                        my $lane_name = "$l_name.$lane";
                        $library->get_lane_by_name($lane_name) || $library->add_lane($lane_name);
                        $hash->{$p_name}->{$s_name}->{$l_name}->{$lane_name} = 1;
                    }
                }
            }
        }
    }
    
    is check_hierarchy($vrtrack, \%hierarchy), 30, 'initial state for time_travel test ok';
    
    sub check_hierarchy {
        my ($vrtrack, $expected) = @_;
        my $count = 0;
        foreach my $project (@{$vrtrack->projects}) {
            $count++;
            foreach my $sample (@{$project->samples}) {
                $count++;
                foreach my $library (@{$sample->libraries}) {
                    $count++;
                    foreach my $lane (@{$library->lanes}) {
                        $count++;
                        ok exists $expected->{$project->name}->{$sample->name}->{$library->name}->{$lane->name}, '  lane '.$lane->name.' in the hierarchy was expected';
                    }
                }
            }
        }
        return $count;
    }
    
    # change the state by adding new projects, samples, libraries and lanes, and
    # also doing a library swap, where all the lanes of one library move to a
    # new library:
    sleep(2);
    my $datestamp = datestamp();
    sleep(2);
    my %new_hierarchy;
    build_hierarchy($vrtrack, 3, \%new_hierarchy);
    my $swap_lib = VRTrack::Library->new_by_name($vrtrack, 'p2.s2.l2');
    my $old_lib = VRTrack::Library->new_by_name($vrtrack, 'p1.s1.l1');
    my %lane_stamps;
    foreach my $lane (@{$old_lib->lanes}) {
        # while we're here, we'll test that update() results in changed()
        # changing
        my $orig_changed = $lane->changed;
        $lane->library_id($swap_lib->id);
        $lane->is_processed('swapped', 1);
        $lane->update;
        my $new_changed = $lane->changed;
        $lane_stamps{$lane->id} = $new_changed;
        isnt $orig_changed, $new_changed, 'an update() causes changed() to change';
        $new_hierarchy{'p2'}->{'p2.s2'}->{'p2.s2.l2'}->{$lane->name} = 1;
        delete $new_hierarchy{'p1'}->{'p1.s1'}->{'p1.s1.l1'}->{$lane->name};
    }
    is check_hierarchy($vrtrack, \%new_hierarchy), 120, 'post state for time_travel test ok';
    
    # and do an extra update so we can later test was_processed works correctly:
    sleep(2);
    foreach my $lane (@{$old_lib->lanes}) {
        $lane->storage_path('/tmp');
        $lane->update;
    }
    
    # can we get back?
    undef $vrtrack;
    $hist->time_travel($datestamp);
    $vrtrack = VRTrack::VRTrack->new($connection_details);
    is check_hierarchy($vrtrack, \%hierarchy), 30, 'time travelling returned us to the initial state';
    
    # test was_processed
    $hist->time_travel('latest');
    foreach my $lane (@{$old_lib->lanes}) {
        is $hist->was_processed($lane, 'swapped'), $lane_stamps{$lane->id};
    }
    
    # test 2 new methods related to working out if a lane "changed"
    is $hist->datetime_cmp('2010-01-04 10:49:10', '2010-01-04 11:40:34'), -1, 'datetime_cmp test old -> recent';
    is $hist->datetime_cmp('2010-02-04 10:49:10', '2010-01-04 11:40:34'), 1, 'datetime_cmp test recent -> old';
    is $hist->datetime_cmp('2010-01-04 10:49:10', '2010-01-04 10:49:10'), 0, 'datetime_cmp test same';
    
    my ($swapped_lane) = @{$old_lib->lanes};
    my $unchanged_lib = VRTrack::Library->new_by_name($vrtrack, 'p1.s2.l1');
    my ($unchanged_lane) = @{$unchanged_lib->lanes};
    is $hist->lane_changed($unchanged_lane, $datestamp), 0, 'lane_changed on an unchanged lane';
    is $hist->lane_changed($swapped_lane, $datestamp), 1, 'lane_changed on a swapped lane';
    
    # and following one more unrelated change, can we get back to the mid-state?
    sleep(2);
    $datestamp = datestamp();
    sleep(2);
    VRTrack::Sample->create($vrtrack, 'p3.s4', 3);
    undef $vrtrack;
    $hist->time_travel($datestamp);
    $vrtrack = VRTrack::VRTrack->new($connection_details);
    is check_hierarchy($vrtrack, \%new_hierarchy), 120, 'time travelling returned us to the middle state';
    
    sub datestamp {
        my $tz = DateTime::TimeZone->new(name => 'local');
        my $dt = DateTime->now(time_zone => $tz);
        return $dt->year.'-'.sprintf("%02d", $dt->month).'-'.sprintf("%02d", $dt->day).' '.sprintf("%02d", $dt->hour).':'.sprintf("%02d", $dt->minute).':'.sprintf("%02d", $dt->second);
    }
    sub new_test_db {
        my @schema = VRTrack::VRTrack->schema();
        open(my $mysqlfh, "| mysql -h$connection_details->{host} -u$connection_details->{user} -p$connection_details->{password}") || die "could not connect to database for testing\n";
        print $mysqlfh "drop database if exists $connection_details->{database};\n";
        print $mysqlfh "create database $connection_details->{database};\n";
        print $mysqlfh "use $connection_details->{database};\n";
        foreach my $sql (@schema) {
            print $mysqlfh $sql;
        }
        close($mysqlfh);
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
