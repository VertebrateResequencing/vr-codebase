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
        plan tests => 134;
    }

    use_ok('VRTrack::VRTrack');
    use_ok('VRTrack::Study');
    use_ok('VRTrack::Sample');
    use_ok('VRTrack::Library');
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

# project
my $project_name = 'Project htest';
my $project_hname = 'Project_htest';
my $project_ssid = 1111;
ok my $vrproj = $vrtrack->add_project($project_name), 'Add project worked';
isa_ok $vrproj, 'VRTrack::Project';
ok my $project_id = $vrproj->id, 'project->id call worked';
is $vrproj->hierarchy_name, $project_hname, 'by default, hname is set to name but made file system safe';
is $vrproj->ssid(), undef, 'project->ssid() returned undef on new project';
is $vrproj->ssid($project_ssid), $project_ssid, 'project->ssid can be set';
ok $vrproj->dirty(),'project dirty after setting ssid';
ok $vrproj->update(),'project updated OK';

{ # check we can get the project back
my $tmp_proj;
ok $tmp_proj = $vrtrack->get_project_by_name($project_name), 'retrieve project by name worked';
is $tmp_proj->id, $project_id, 'project_by_name is the right project';

ok $tmp_proj = $vrtrack->get_project_by_ssid($project_ssid), 'retrieve project by ssid worked';
is $tmp_proj->id, $project_id, 'project_by_ssid is the right project';

ok $tmp_proj = $vrtrack->get_project_by_id($project_id), 'retrieve project by id worked';
is $tmp_proj->id, $project_id, 'project_by_id is the right project';
}

# Study
my $study_acc = 'Study_test';
is $vrproj->study(), undef, 'No study on new project';
is $vrproj->study($study_acc), undef, 'study() returned undef when setting a non-existant study';

is $vrproj->study_id('1234'), '1234', 'study_id could be set';

ok my $study = VRTrack::Study->create($vrtrack, $study_acc), 'Study->create worked';
my $study_id = $study->id;
is $vrproj->study($study_acc)->id, $study_id, 'study() returned the correct study obj after setting it';
is $vrproj->study_id, $study_id, 'study_id was updated';

$study_acc = 'Study_test2';
ok my $study2 = $vrproj->add_study($study_acc), 'add_study worked';
isnt $study2->id, $study_id, 'add_study made a new study';
is $vrproj->study_id, $study2->id, 'study_id was updated again';

#SAMPLE
my $sample_1 = $vrproj->add_sample('sample_a');
isa_ok $sample_1, 'VRTrack::Sample', 'add_sample returned a sample';
my $sample_2 = $vrproj->add_sample('sample_b');
isa_ok $sample_2, 'VRTrack::Sample', 'add_sample returned a sample again';
isnt $sample_1->id, $sample_2->id, 'samples are different';
$vrproj->update;

is @{$vrproj->samples}, 2, 'samples returned the correct number of samples';
cmp_bag ($vrproj->sample_ids, [$sample_1->id, $sample_2->id], 'sample_ids returned the correct ids');
my $sample_name = 'sample_a';
is $vrproj->get_sample_by_name($sample_name)->id, $sample_1->id, 'get_sample_by_name worked';
is $vrproj->get_sample_by_id($sample_1->id)->name, $sample_name, 'get_sample_by_id worked';
ok my $sample = VRTrack::Sample->new_by_name_project($vrtrack, $sample_name, $project_id), 'sample new_by_name_project worked';
is $sample->id(),$sample_1->id,'sample id is correct';

my $ind_name = 'Individual_test';
is $sample->individual($ind_name), undef, 'individual() returned undef when setting a non-existant individual';
ok my $ind = $sample->add_individual($ind_name), 'add_individual worked';
is $sample->individual_id, $ind->id, 'individual_id was updated';
my $lib_1 = $sample->add_library('lib_a');
isa_ok $lib_1, 'VRTrack::Library', 'add_library returned a library';
my $lib_2 = $sample->add_library('lib_b');
isa_ok $lib_2, 'VRTrack::Library', 'add_library returned a library again';
isnt $lib_1->id, $lib_2->id, 'libs are different';

ok $sample->dirty(),'sample dirty after adding libs';
ok $sample->update(), 'sample updated OK';

is @{$sample->libraries}, 2, 'libraries returned the correct number of libraries';
cmp_bag ($sample->library_ids, [$lib_1->id, $lib_2->id], 'library_ids returned the correct ids');
my $library_name = 'lib_a';
is $sample->get_library_by_name($library_name)->id, $lib_1->id, 'get_library_by_name worked';
is $sample->get_library_by_id($lib_1->id)->name, $library_name, 'get_library_by_id worked';
ok my $library_request=$sample->add_library_request('1234'),'add library_request to sample worked';

my $sample_id = $sample->id;

# library
ok my $library = VRTrack::Library->new_by_name($vrtrack, $library_name), 'library new_by_name worked';
is $library->id, $lib_1->id, 'library is the correct library';
is $library->sample_id(), $sample_id, "library linked to sample $sample_id correctly";

# lib_type
my $ltype_name = 'library_type_test';
is $library->library_type($ltype_name), undef, 'library_type() returned undef when setting a non-existant library type';
is $library->library_type_id, undef, 'library_type_id starts undefined';
my $lib_type = $library->add_library_type($ltype_name);
isa_ok $lib_type, 'VRTrack::Library_type', 'add_library_type returned a library_type';
is $library->library_type_id, $lib_type->id, 'library_type_id was updated';

# lib_tag, group, sequence (multiplex index tag)
is $library->library_tag(2),2,'setting library_tag worked';
is $library->library_tag_group(3),3,'setting library_tag_group worked';
is $library->library_tag_sequence('ACGTCTCT'),'ACGTCTCT','setting library_tag_sequence worked';

# seq_centre
my $sc_name = 'seq_centre_test';
is $library->seq_centre($sc_name), undef, 'seq_centre() returned undef when setting a non-existant seq_centre';
is $library->seq_centre_id, undef, 'seq_centre_id starts undefined';
my $seq_centre = $library->add_seq_centre($sc_name);
isa_ok $seq_centre, 'VRTrack::Seq_centre', 'add_seq_centre made a new seq_centre';
is $library->seq_centre_id, $seq_centre->id, 'seq_centre_id was updated';

# seq tech
my $st_name = 'seq_tech_test';
is $library->seq_tech($st_name), undef, 'seq_tech() returned undef when setting a non-existant seq_tech';
is $library->seq_tech_id, undef, 'seq_tech_id starts undefined';
my $seq_tech = $library->add_seq_tech($st_name);
isa_ok $seq_tech, 'VRTrack::Seq_tech', 'add_seq_tech made a new seq_tech';
is $library->seq_tech_id, $seq_tech->id, 'seq_tech_id was updated';

# LANE

my $lane_1 = $library->add_lane('lane_a');
isa_ok $lane_1, 'VRTrack::Lane', 'add_lane returned a lane';
my $lane_2 = $library->add_lane('lane_b');
isa_ok $lane_2, 'VRTrack::Lane', 'add_lane returned a lane again';
isnt $lane_1->id, $lane_2->id, 'lanes are different';

ok $library->update;

is @{$library->lanes}, 2, 'lanes returned the correct number of lanes';
cmp_bag ($library->lane_ids, [$lane_1->id, $lane_2->id], 'lane_ids returned the correct ids');
my $lane_name = 'lane_a';
is $library->get_lane_by_name($lane_name)->id, $lane_1->id, 'get_lane_by_name worked';
is $library->get_lane_by_id($lane_1->id)->name, $lane_name, 'get_lane_by_id worked';

ok my $vrlane = VRTrack::Lane->new_by_name($vrtrack, $lane_name), 'lane new_by_name worked';
is $vrlane->id, $lane_1->id, 'lane is the right lane';
is $vrlane->qc_status(), undef, 'lane qc_status defaults undef';
is $vrlane->qc_status('passed'), 'passed', 'lane qc_status could be set to passed';
is $vrlane->is_processed('import'), 0, 'is_processed import starts off';
is $vrlane->is_processed('import', 1), 1, 'is_processed import could be turned on';
is $vrlane->is_processed('mapped', 1), 1, 'is_processed mapped could be turned on';
is $vrlane->processed(), 5, 'actual processed number is correct';

my $sub_name = 'submission_test';
is $vrlane->submission($sub_name), undef, 'submission() returned undef when setting a non-existant submission';
is $vrlane->submission_id, undef, 'submission_id starts undefined';

# submission
my $sub = $vrlane->add_submission($sub_name);
isa_ok $sub, 'VRTrack::Submission', 'add_submission returns a submission';
is $vrlane->submission_id, $sub->id, 'submission_id was updated';

# mapping
my $map_1 = $vrlane->add_mapping();
isa_ok $map_1, 'VRTrack::Mapstats', 'add_mapping returns a mapstats';
my $map_2 = $vrlane->add_mapping();
isa_ok $map_2, 'VRTrack::Mapstats', 'add_mapping returns a mapstats again';
ok $vrlane->update;
is @{$vrlane->mappings}, 2, 'mappings returned the correct number of mappings';
cmp_bag ($vrlane->mapping_ids, [$map_1->id, $map_2->id], 'mapping_ids returned the correct ids');
is $vrlane->latest_mapping->id, $map_2->id, 'latest_mapping got the last mapstats';

ok my $mapping = VRTrack::Mapstats->new($vrtrack, $map_1->id), 'was able to retrieve mapstats object by id';
is $mapping->raw_reads(999), 999, 'Mapstats->raw_reads could be set';
is $mapping->raw_bases(111), 111, 'Mapstats->raw_bases could be set';
my $img_filename = File::Spec->catfile(qw(t data io_test.txt));
my $img_1 = $mapping->add_image_by_filename($img_filename);
isa_ok $img_1, 'VRTrack::Image','was able to add an image by filename to the mapstats object';
my $img_2 = $mapping->add_image_by_filename($img_filename);
isa_ok $img_2, 'VRTrack::Image','add_image_by_filename works again, even given the same file';
isnt $img_1->id, $img_2->id, 'images are different';

ok $mapping->update, 'changes could be written to db';
is @{$mapping->images}, 2, 'images returned the correct number of images';
cmp_bag ($mapping->image_ids, [$img_1->id, $img_2->id], 'image_ids returned the correct ids');
my $mapper_name = 'mapper_test';
my $mapper_version = 0.5;
is $mapping->mapper($mapper_name, $mapper_version), undef, 'mapper() returned undef when setting a non-existant mapper';
is $mapping->mapper_id, undef, 'mapper_id starts undefined';
my $mapper = $mapping->add_mapper($mapper_name, $mapper_version);
isa_ok $mapper, 'VRTrack::Mapper', 'add_mapper made a new mapper';
is $mapping->mapper_id, $mapper->id, 'mapper_id was updated';
my $assembly_name = 'assembly_test';
is $mapping->assembly($assembly_name), undef, 'assembly() returned undef when setting a non-existant assembly';
is $mapping->assembly_id, undef, 'assembly_id starts undefined';
my $assy = $mapping->add_assembly($assembly_name);
isa_ok $assy, 'VRTrack::Assembly', 'add_assembly made a new assembly';
is $mapping->assembly_id, $assy->id, 'assembly_id was updated';
my $assembly_name2 = $assembly_name.'_2';
my $assy2 = $mapping->add_assembly($assembly_name2);
isnt $assy2->id, $assy->id, 'add_assembly used again replaces the old assembly';
is $mapping->assembly_id, $assy2->id, 'assembly_id was updated';
is $mapping->assembly->id, $assy2->id, 'assembly() returns the new second assembly';
is $mapping->assembly_id($assy->id), $assy->id, 'assembly_id could be changed back to the first assembly';
is $mapping->assembly->id, $assy->id, 'assembly() returns the first assembly';

# FILE
my $file_name = 'File/test_1.fastq';
ok ! $vrlane->get_file_by_name($file_name), 'Lane->get_file_by_name returns undef before we\'ve added any files';
my $file_1 = $vrlane->add_file($file_name);
isa_ok $file_1, 'VRTrack::File', 'add_file returns file';
is_deeply [$file_1->name, $file_1->hierarchy_name], [$file_name, 'File_test_1.fastq'], 'added files have the correct name and hierarchy_name';

is $vrlane->get_file_by_name($file_name)->id, $file_1->id, 'get_file_by_name now works';

my $file_2 = $vrlane->add_file('File_test_2.fastq');
isa_ok $file_2, 'VRTrack::File', 'add_file worked again';
isnt $file_1->id, $file_2->id, 'files are different';
$vrlane->update;

is @{$vrlane->files}, 2, 'files returned the correct number of files';
cmp_bag ($vrlane->file_ids, [$file_1->id, $file_2->id], 'file_ids returned the correct ids');

ok my $vrfile = VRTrack::File->new($vrtrack, $file_1->id), 'was able to retrieve file object by id';

is $vrfile->type(0), 0, 'File type could be set to 0';
is $vrfile->type(1), 1, 'File type could be set to 1';
ok $vrfile->update, 'File update was successful';
undef($vrfile);
$vrfile = VRTrack::File->new($vrtrack, $file_1->id);
is $vrfile->type, 1, 'File type really was set to 1';

# hierarchy_path_of_lane
$ENV{DATA_HIERARCHY} = '';
is $vrtrack->hierarchy_path_of_lane($vrlane), 'Project_htest/sample_a/seq_tech_test/lib_a/lane_a', 'hierarchy_path_of_lane works with DATA_HIERARCHY unset';
my $lane_b = VRTrack::Lane->new_by_name($vrtrack, 'lane_b');
my $lib_b = $sample->get_library_by_id($lib_2->id);
$lib_b->add_seq_tech('seq_tech_test_b');
$lib_b->update;
$lane_b->library_id($lib_b->id);
$lane_b->update;
is $vrtrack->hierarchy_path_of_lane($lane_b), 'Project_htest/sample_a/seq_tech_test_b/lib_b/lane_b', 'hierarchy_path_of_lane works with DATA_HIERARCHY unset and on a different lane';
$ENV{DATA_HIERARCHY} = 'species:foo:library';
is $vrtrack->hierarchy_path_of_lane($vrlane), 'species/foo/lib_a', 'hierarchy_path_of_lane works with DATA_HIERARCHY set to species:library';

# more tests for hierarchy_path_of_lane 
# my $vrlane = VRTrack::Lane->new_by_name($vrtrack, 'lane_a'); #Get our favourite lane back. At this point all the bits of data exist in the database
# $ENV{DATA_HIERARCHY} = 'species:foo:library';
# is $vrtrack->hierarchy_path_of_lane($vrlane), 'species/foo/lib_a', 'hierarchy_path_of_lane works with DATA_HIERARCHY set to species:foo:library';
# $ENV{DATA_HIERARCHY} = 'genus:species-subspecies:project:strain:sample:technology:library:lane';
# is $vrtrack->hierarchy_path_of_lane($vrlane) ,'genus/species_subspecies/Project_test2/Individual_test/sample_a/seq_tech_test/lib_a/lane_a', 'hierarchy_path_of_lane works with DATA_HIERARCHY set to all items';

# test descendants using lane
$vrlane->update;
my @children = @{$vrlane->descendants};
is_deeply [sort map { ref($_) } @children], ['VRTrack::File', 'VRTrack::File', 'VRTrack::Image', 'VRTrack::Image', 'VRTrack::Mapstats', 'VRTrack::Mapstats'], 'got correct descendants for the lane';

# cleanup
my $dbh = $vrtrack->{_dbh};
foreach ($dbh->tables()){
    next if /`latest_/;
    next if /schema_version/;
    $dbh->do("TRUNCATE TABLE $_");
}
exit;
