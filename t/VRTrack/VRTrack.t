#!/usr/bin/env perl
use strict;
use warnings;
use DateTime;

BEGIN {
    use Test::Most tests => 605;

    use_ok('VRTrack::VRTrack');
    use_ok('VRTrack::Request');
    use_ok('VRTrack::Study');
    use_ok('VRTrack::History');
}

my $connection_details = { database => 'vrtrack_test',
                           host     => 'mcs4a',
                           port     => 3306,
                           user     => 'vreseq_rw',
                           password => 't3aml3ss' };

# drop, create and then populate our testing database with the VRTrack schema
# *** need a better way of doing this...
ok my @schema = VRTrack::VRTrack->schema(), 'schema() returned something';
new_test_db();

sub new_test_db {
    open(my $mysqlfh, "| mysql -h$connection_details->{host} -u$connection_details->{user} -p$connection_details->{password}") || die "could not connect to database for testing\n";
    print $mysqlfh "drop database if exists $connection_details->{database};\n";
    print $mysqlfh "create database $connection_details->{database};\n";
    print $mysqlfh "use $connection_details->{database};\n";
    foreach my $sql (@schema) {
        print $mysqlfh $sql;
    }
    close($mysqlfh);
}

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
    
    # shared by some, but not all
    if ($vrobj->can('hierarchy_name')) {
        isa_ok $vrobj, "VRTrack::Hierarchy_obj";
        
        is $vrobj->hierarchy_name, $name, 'hierarchy_name set to name by default';
        is $vrobj->hierarchy_name($new_name), $new_name, 'hierarchy_name could be reset, but it was not made file-system safe';
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

# core objects keep a history
my $name = 'Project_test2';
ok my $vrproj = VRTrack::Project->create($vrtrack, $name), 'Project->create worked again';
is_deeply [$vrproj->id, $vrproj->row_id], [3, 3], 'ids were incremented for the new object';
my $changed = $vrproj->changed;
sleep(1);
$vrproj->ssid(123);
ok $vrproj->update(), 'after changing an attribute was able to update';
$vrproj = VRTrack::Project->new_by_name($vrtrack, $name);
my $changed2 = $vrproj->changed;
is_deeply [$vrproj->id, $vrproj->row_id], [3, 4], 'getting it back with new_by_name shows the id has not changed whilst the row_id has incremented';
$vrproj = VRTrack::Project->new($vrtrack, 3, $changed2);
is_deeply [$vrproj->id, $vrproj->row_id, $vrproj->changed], [3, 3, $changed], 'new(... new date) gives us the old version of the object';
$vrproj = VRTrack::Project->new($vrtrack, 3, 'latest');
is_deeply [$vrproj->id, $vrproj->row_id, $vrproj->changed], [3, 4, $changed2], 'new(... latest) gives us the latest version of the object';
is $vrproj->global_history_date($changed2), $changed2, 'global_history_date could be set';
$vrproj = VRTrack::Project->new($vrtrack, 3);
is_deeply [$vrproj->id, $vrproj->row_id, $vrproj->changed], [3, 3, $changed], 'new() gives us the old version of the object when global was set';
$vrproj = VRTrack::Project->new($vrtrack, 3, 'latest');
is_deeply [$vrproj->id, $vrproj->row_id, $vrproj->changed], [3, 4, $changed2], 'new(... latest) gives us the latset version inspite of global setting';
is $vrproj->global_history_date('latest'), 'latest', 'global_history_date could be set back to latest';
$vrproj = VRTrack::Project->new_by_name($vrtrack, $name, 3);
is_deeply [$vrproj->id, $vrproj->row_id, $vrproj->changed], [3, 3, $changed], 'new_by_name(... old row_id) gives us the old version of the object';
$vrproj = VRTrack::Project->new_by_name($vrtrack, $name);
is_deeply [$vrproj->id, $vrproj->row_id, $vrproj->changed], [3, 4, $changed2], 'new_by_name() still gives us the latest version of the object';
is_deeply [$vrproj->row_ids], [3, 4], 'row_ids worked';
$vrproj = VRTrack::Project->new_by_name($vrtrack, $name, 3);    # get old object again
$vrproj->ssid(321);
is $vrproj->ssid, 321, 'can change attributes on historical objects'; # this is mostly useless, except for is_latest(1), of course
ok ! $vrproj->update(), 'can\'t update on historical objects after changing attribute';


# test bits of the core classes that are unique to themselves
{
    # project
    my $study_id = '5634';
    is $vrproj->study_id($study_id), $study_id, 'study_id could be set';
    my $study_acc = 'Study_test';
    VRTrack::Study->create($vrtrack, $study_acc);
    is $vrproj->study($study_acc)->id, 1, 'study() returned the correct study obj after setting it';
    $study_acc = 'Study_test2';
    is $vrproj->study($study_acc), undef, 'study() returned undef when setting a non-existant study';
    is $vrproj->study_id, 1, 'study_id was updated';
    is $vrproj->add_study($study_acc)->id, 2, 'add_study made a new study';
    is $vrproj->study_id, 2, 'study_id was updated again';
    is $vrproj->add_sample('sample_a')->id, 3, 'add_sample worked';
    is $vrproj->add_sample('sample_b')->id, 5, 'add_sample worked again';
    $vrproj->update;
    is @{$vrproj->samples}, 2, 'samples returned the correct number of samples';
    is_deeply $vrproj->sample_ids, [3, 5], 'sample_ids returned the correct ids';
    my $sample_name = 'sample_a';
    is $vrproj->get_sample_by_name($sample_name)->id, 3, 'get_sample_by_name worked';
    is $vrproj->get_sample_by_id(3)->name, $sample_name, 'get_sample_by_id worked';
    # when a new project is added, it gets a default hname that is file-system
    # safe
    my $added_p = $vrtrack->add_project('added|project');
    is $added_p->name, 'added|project', 'newly added project name is same as set';
    is $added_p->hierarchy_name, 'added_project', 'by default, hname is set to name but made file system safe';
    
    # sample
    ok my $sample = VRTrack::Sample->new_by_name_project($vrtrack, $sample_name, 3), 'sample new_by_name_project worked';
    my $ind_name = 'Individual_test';
    is $sample->individual($ind_name), undef, 'individual() returned undef when setting a non-existant individual';
    is $sample->add_individual($ind_name)->id, 1, 'add_individual worked';
    is $sample->individual_id, 1, 'individual_id was updated';
    is $sample->add_library('lib_a')->id, 3, 'add_library worked';
    is $sample->add_library('lib_b')->id, 5, 'add_library worked again';
    $sample->update;
    is @{$sample->libraries}, 2, 'libraries returned the correct number of libraries';
    is_deeply $sample->library_ids, [3, 5], 'library_ids returned the correct ids';
    my $library_name = 'lib_a';
    is $sample->get_library_by_name($library_name)->id, 3, 'get_library_by_name worked';
    is $sample->get_library_by_id(3)->name, $library_name, 'get_library_by_id worked';
    
    # library
    ok my $library = VRTrack::Library->new_by_name($vrtrack, $library_name), 'library new_by_name worked';
    my $ltype_name = 'library_type_test';
    is $library->library_type($ltype_name), undef, 'library_type() returned undef when setting a non-existant library type';
    is $library->library_type_id, undef, 'library_type_id starts undefined';
    is $library->add_library_type($ltype_name)->id, 1, 'add_library_type made a new library type';
    is $library->library_type_id, 1, 'library_type_id was updated';
    my $sc_name = 'seq_centre_test';
    is $library->seq_centre($sc_name), undef, 'seq_centre() returned undef when setting a non-existant seq_centre';
    is $library->seq_centre_id, undef, 'seq_centre_id starts undefined';
    is $library->add_seq_centre($sc_name)->id, 1, 'add_seq_centre made a new seq_centre';
    is $library->seq_centre_id, 1, 'seq_centre_id was updated';
    my $st_name = 'seq_tech_test';
    is $library->seq_tech($st_name), undef, 'seq_tech() returned undef when setting a non-existant seq_tech';
    is $library->seq_tech_id, undef, 'seq_tech_id starts undefined';
    is $library->add_seq_tech($st_name)->id, 1, 'add_seq_tech made a new seq_tech';
    is $library->seq_tech_id, 1, 'seq_tech_id was updated';
    is $library->add_lane('lane_a')->id, 3, 'add_lane worked';
    is $library->add_lane('lane_b')->id, 5, 'add_lane worked again';
    $library->update;
    is @{$library->lanes}, 2, 'lanes returned the correct number of lanes';
    is_deeply $library->lane_ids, [3, 5], 'lane_ids returned the correct ids';
    my $lane_name = 'lane_a';
    is $library->get_lane_by_name($lane_name)->id, 3, 'get_lane_by_name worked';
    is $library->get_lane_by_id(3)->name, $lane_name, 'get_lane_by_id worked';
    
    # request
    # *** (can't test since can't create...)
    
    # lane
    my $vrlane = VRTrack::Lane->new_by_name($vrtrack, $lane_name);
    is $vrlane->qc_status(), undef, 'lane qc_status defaults undef';
    is $vrlane->qc_status('passed'), 'passed', 'lane qc_status could be set to passed';
    is $vrlane->is_processed('import'), 0, 'is_processed import starts off';
    is $vrlane->is_processed('import', 1), 1, 'is_processed import could be turned on';
    is $vrlane->is_processed('mapped', 1), 1, 'is_processed mapped could be turned on';
    is $vrlane->processed(), 5, 'actual processed number is correct';
    my $sub_name = 'submission_test';
    is $vrlane->submission($sub_name), undef, 'submission() returned undef when setting a non-existant submission';
    is $vrlane->submission_id, undef, 'submission_id starts undefined';
    is $vrlane->add_submission($sub_name)->id, 1, 'add_submission made a new submission';
    is $vrlane->submission_id, 1, 'submission_id was updated';
    is $vrlane->add_mapping()->id, 3, 'add_mapping worked';
    is $vrlane->add_mapping()->id, 5, 'add_mapping worked again';
    $vrlane->update;
    is @{$vrlane->mappings}, 2, 'mappings returned the correct number of mappings';
    is_deeply $vrlane->mapping_ids, [3, 5], 'mapping_ids returned the correct ids';
    is $vrlane->latest_mapping->id, 5, 'latest_mapping got the last mapstats';
    my $file_name = 'File/test_1.fastq';
    ok ! $vrlane->get_file_by_name($file_name), 'Lane->get_file_by_name returns undef before we\'ve added any files';
    my $file = $vrlane->add_file($file_name);
    is $file->id, 3, 'add_file worked';
    is_deeply [$file->name, $file->hierarchy_name], [$file_name, 'File_test_1.fastq'], 'added files have the correct name and hierarchy_name';
    is $vrlane->get_file_by_name($file_name)->id, 3, 'get_file_by_name now works';
    is $vrlane->add_file('File_test_2.fastq')->id, 5, 'add_file worked again';
    $vrlane->update;
    is @{$vrlane->files}, 2, 'files returned the correct number of files';
    is_deeply $vrlane->file_ids, [3, 5], 'file_ids returned the correct ids';
    
    # hierarchy_path_of_lane
    $ENV{DATA_HIERARCHY} = '';
    is $vrtrack->hierarchy_path_of_lane($vrlane), 'Project_test2/sample_a/seq_tech_test/lib_a/lane_a', 'hierarchy_path_of_lane works with DATA_HIERARCHY unset';
    my $lane_b = VRTrack::Lane->new_by_name($vrtrack, 'lane_b');
    my $lib_b = $sample->get_library_by_id(5);
    $lib_b->add_seq_tech('seq_tech_test_b');
    $lib_b->update;
    $lane_b->library_id(5);
    $lane_b->update;
    is $vrtrack->hierarchy_path_of_lane($lane_b), 'Project_test2/sample_a/seq_tech_test_b/lib_b/lane_b', 'hierarchy_path_of_lane works with DATA_HIERARCHY unset and on a different lane';
    $ENV{DATA_HIERARCHY} = 'species:foo:library';
    is $vrtrack->hierarchy_path_of_lane($vrlane), 'species/foo/lib_a', 'hierarchy_path_of_lane works with DATA_HIERARCHY set to species:library';
    
    # file
    ok my $vrfile = VRTrack::File->new($vrtrack, 3), 'was able to retrieve file object by id';
    is $vrfile->type(1), 1, 'File type could be set to 1';
    is $vrfile->type(0), 0, 'File type could be set to 0';
    # *** write some more tests... (they're just simple get/setters)
    ok $vrfile->update, 'File update was successful';
    $vrfile = VRTrack::File->new($vrtrack, 3);
    is $vrfile->type, 0, 'File type really was set to 0';
    
    # mapstats
    ok my $mapping = VRTrack::Mapstats->new($vrtrack, 3), 'was able to retrieve mapstats object by id';
    is $mapping->raw_reads(999), 999, 'Mapstats->raw_reads could be set';
    is $mapping->raw_bases(111), 111, 'Mapstats->raw_bases could be set';
    my $img_filename = File::Spec->catfile(qw(t data io_test.txt));
    is $mapping->add_image_by_filename($img_filename)->id, 1, 'was able to add an image by filename to the mapstats object';
    is $mapping->add_image_by_filename($img_filename)->id, 2, 'add_image_by_filename works again, even given the same file';
    ok $mapping->update, 'changes could be written to db';
    is @{$mapping->images}, 2, 'images returned the correct number of images';
    is_deeply $mapping->image_ids, [1, 2], 'image_ids returned the correct ids';
    my $mapper_name = 'mapper_test';
    my $mapper_version = 0.5;
    is $mapping->mapper($mapper_name, $mapper_version), undef, 'mapper() returned undef when setting a non-existant mapper';
    is $mapping->mapper_id, undef, 'mapper_id starts undefined';
    is $mapping->add_mapper($mapper_name, $mapper_version)->id, 1, 'add_mapper made a new mapper';
    is $mapping->mapper_id, 1, 'mapper_id was updated';
    my $assembly_name = 'assembly_test';
    is $mapping->assembly($assembly_name), undef, 'assembly() returned undef when setting a non-existant assembly';
    is $mapping->assembly_id, undef, 'assembly_id starts undefined';
    is $mapping->add_assembly($assembly_name)->id, 1, 'add_assembly made a new assembly';
    is $mapping->assembly_id, 1, 'assembly_id was updated';
    my $assembly_name2 = $assembly_name.'_2';
    is $mapping->add_assembly($assembly_name2)->id, 2, 'add_assembly used again replaces the old assembly';
    is $mapping->assembly_id, 2, 'assembly_id was updated';
    is $mapping->assembly->id, 2, 'assembly() returns the new second assembly';
    is $mapping->assembly_id(1), 1, 'assembly_id could be changed back to the first assembly';
    is $mapping->assembly->id, 1, 'assembly() returns the first assembly';
    
    # test descendants using lane
    $vrlane->update;
    my @children = @{$vrlane->descendants};
    is_deeply [sort map { ref($_) } @children], ['VRTrack::File', 'VRTrack::File', 'VRTrack::Image', 'VRTrack::Image', 'VRTrack::Mapstats', 'VRTrack::Mapstats'], 'got correct descendants for the lane';
}

# tests for non-core classes
{
    # allocation
    # *** not a Table_obj (yet?), don't feel like testing it...
    
    # image
    my $image = VRTrack::Image->new($vrtrack, 1);
    is $image->id, 1, 'Image new returned the first image we made before';
    my $image_name = 'image_test_changed';
    is $image->name($image_name), $image_name, 'image name could be get/set';
    is $image->caption('foo'), 'foo', 'image caption could be get/set';
    ok $image->update, 'update worked on image';
    $image = VRTrack::Image->new($vrtrack, 1);
    is_deeply [$image->id, $image->name, $image->caption], [1, $image_name, 'foo'], 'after update, the things really were set';
    $image_name .= '2';
    ok $image = VRTrack::Image->create($vrtrack, $image_name, 'bar'), 'image create worked';
    $image = VRTrack::Image->new($vrtrack, 3);
    is_deeply [$image->id, $image->name], [3, $image_name], 'Image new worked again';
    
    # individual
    # setting species_id and population id for the first individual
    my $individual = VRTrack::Individual->new($vrtrack, 1);
    is $individual->id, 1, 'Individual new returned the first assembly we made before';
    my $individual_name = 'individual_test|changed';
    is $individual->name($individual_name), $individual_name, 'Individual name could be get/set';
    is $individual->sex('M'), 'M', 'Individual sex could be get/set to M';
    is $individual->sex('F'), 'F', 'Individual sex could be get/set to F';
    is $individual->sex('foobar'), 'unknown', 'Individual sex is set to unknown for non M/F';
    ok $individual->update, 'update worked on Individual';
    $individual = VRTrack::Individual->new_by_name($vrtrack, $individual_name);
    is_deeply [$individual->id, $individual->name, $individual->sex], [1, $individual_name, 'unknown'], 'Individual new_by_name worked';
    my $species_name = 'genus species subspecies';
    is $individual->species($species_name), undef, 'species() returned undef when setting a non-existant species';
    is $individual->species_id, undef, 'species_id starts undefined';
    is $individual->add_species($species_name)->id, 1, 'add_species made a new species';
    is $individual->species_id, 1, 'species_id was updated';
    my $pop_name = 'population_test';
    is $individual->population($pop_name), undef, 'population() returned undef when setting a non-existant population';
    is $individual->population_id, undef, 'population_id starts undefined';
    is $individual->add_population($pop_name)->id, 1, 'add_population made a new population';
    is $individual->population_id, 1, 'population_id was updated';

    $individual->update();

    $individual_name .= '2';
    ok $individual = VRTrack::Individual->create($vrtrack, $individual_name), 'Individual create worked';
    $individual = VRTrack::Individual->new_by_name($vrtrack, $individual_name);
    is_deeply [$individual->id, $individual->name, $individual->hierarchy_name], [2, $individual_name, 'individual_test_changed2'], 'Individual new_by_name worked again';
    
    # more tests for hierarchy_path_of_lane 
    my $vrlane = VRTrack::Lane->new_by_name($vrtrack, 'lane_a'); #Get our favourite lane back. At this point all the bits of data exist in the database
    $ENV{DATA_HIERARCHY} = 'species:foo:library';
    is $vrtrack->hierarchy_path_of_lane($vrlane), 'species/foo/lib_a', 'hierarchy_path_of_lane works with DATA_HIERARCHY set to species:foo:library';
    $ENV{DATA_HIERARCHY} = 'genus:species-subspecies:project:strain:sample:technology:library:lane';
    is $vrtrack->hierarchy_path_of_lane($vrlane) ,'genus/species_subspecies/Project_test2/Individual_test/sample_a/seq_tech_test/lib_a/lane_a', 'hierarchy_path_of_lane works with DATA_HIERARCHY set to all items';
    
    # mapper
    my $mapper = VRTrack::Mapper->new($vrtrack, 1);
    is $mapper->id, 1, 'Mapper new returned the first Mapper we made before';
    my $mapper_name = 'mapper_test_changed';
    is $mapper->name($mapper_name), $mapper_name, 'Mapper name could be get/set';
    ok $mapper->update, 'update worked on Mapper';
    $mapper = VRTrack::Mapper->new_by_name_version($vrtrack, $mapper_name, 0.5);
    is_deeply [$mapper->id, $mapper->name], [1, $mapper_name], 'Mapper new_by_name_version worked';
    $mapper_name .= '2';
    ok $mapper = VRTrack::Mapper->create($vrtrack, $mapper_name, 3), 'Mapper create worked';
    $mapper = VRTrack::Mapper->new_by_name_version($vrtrack, $mapper_name, 3);
    is_deeply [$mapper->id, $mapper->name, $mapper->version], [2, $mapper_name, 3], 'Mapper new_by_name_version worked again';
    
    # study
    my $study = VRTrack::Study->new($vrtrack, 1);
    is $study->id, 1, 'Study new returned the first study we made before';
    my $study_acc = 'study_test_changed';
    is $study->acc($study_acc), $study_acc, 'study acc could be get/set';
    ok $study->update, 'update worked on study';
    $study = VRTrack::Study->new_by_acc($vrtrack, $study_acc);
    is_deeply [$study->id, $study->acc], [1, $study_acc], 'Study new_by_acc worked';
    $study_acc .= '2';
    ok $study = VRTrack::Study->create($vrtrack, $study_acc), 'study create worked';
    $study = VRTrack::Study->new_by_acc($vrtrack, $study_acc);
    is_deeply [$study->id, $study->acc], [3, $study_acc], 'Study new_by_acc worked again';

    # species
    my $species = VRTrack::Species->create($vrtrack,"organism");
    is $species->taxon_id(),0, 'Created a species entry with taxon id 0 as expected';
    my $other_species = VRTrack::Species->create($vrtrack,"organism_2",456);
    is $other_species->taxon_id(),456, 'Created a species entry with taxon id 456 as expected';
    
    
    # common tests for the classes that have just an id and name
    my %class_to_next_id = (Assembly => 3,
                            Library_type => 2,
                            Population => 2,
                            Seq_centre => 2,
                            Seq_tech => 3,
                            # *** actually the next two have some simple get/setters as well
                            Species => 4, 
                            Submission => 2);
    
    while (my ($class, $next_id) = each %class_to_next_id) {
        my $obj = "VRTrack::$class"->new($vrtrack, 1);
        is $obj->id, 1, "$class new returned the first $class we made before";
        my $obj_name = "${class}_test_changed";
        is $obj->name($obj_name), $obj_name, "$class name could be get/set";
        ok $obj->update, "update worked on $class";
        $obj = "VRTrack::$class"->new_by_name($vrtrack, $obj_name);
        is_deeply [$obj->id, $obj->name], [1, $obj_name], "$class new_by_name worked";
        $obj_name .= '2';
        ok $obj = "VRTrack::$class"->create($vrtrack, $obj_name), "$class create worked";
        $obj = "VRTrack::$class"->new_by_name($vrtrack, $obj_name);
        is_deeply [$obj->id, $obj->name], [$next_id, $obj_name], "$class new_by_name worked again";
    }
}

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
foreach my $lane_name (qw(history_lane_a history_lane_b history_lane_c history_lane_d)) {
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
$vrproj->ssid(123);
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
}

exit;
