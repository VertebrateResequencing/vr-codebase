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
        plan tests => 63;
    }

    use_ok('VRTrack::VRTrack');
    use_ok('VRTrack::Image');
    use_ok('VRTrack::AutoQC');
    use_ok('VRTrack::Individual');
    use_ok('VRTrack::Lane');
    use_ok('VRTrack::Mapper');
    use_ok('VRTrack::Study');
    use_ok('VRTrack::Species');
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
# allocation
# TODO

# image
my $image_name = 'image_test_orig';
my $imagedata = 'test_binaryimagedata';
ok my $image = VRTrack::Image->create($vrtrack, $image_name, $imagedata), 'image create worked';
my $image_id = $image->id;
$image = VRTrack::Image->new($vrtrack, $image_id);
is $image->id, $image_id, 'Image new returned the image we made before';
is $image->image(),$imagedata, 'image data stored correctly';
$image_name = 'image_test_changed';
is $image->name($image_name), $image_name, 'image name could be get/set';
is $image->caption('foo'), 'foo', 'image caption could be get/set';
ok $image->update, 'update worked on image';
$image = VRTrack::Image->new($vrtrack, $image_id);
is_deeply [$image->id, $image->name, $image->image], [$image_id, $image_name, $imagedata], 'after update, the things really were set';

# AutoQC 
my $reason='The lane failed the NPG QC check, so auto-fail';
ok my $autoqc = VRTrack::AutoQC->create($vrtrack, 'NPG QC status', 0, $reason), 'autoqc create worked';
my $autoqc_id = $autoqc->id;
$autoqc = VRTrack::AutoQC->new($vrtrack, $autoqc_id);
is $autoqc->id, $autoqc_id, 'AutoQC new returned the autoqc we made before';
is $autoqc->reason(),$reason, 'autoqc reason stored correctly';

# test autoqc getters/setters
my $qc_test = 'autoqc_qc_test_name_changed';
is $autoqc->test($qc_test), $qc_test, 'autoqc test name could be get/set';
my $qc_reason= 'modified test result reason';
is $autoqc->reason($qc_reason), $qc_reason, 'autoqc reason name could be get/set';
my $qc_result= 0;
is $autoqc->result($qc_result), $qc_result, 'autoqc result could be updated';
ok $autoqc->update, 'update worked on autoqc';

# Add autoqcs to mapping stats
my $mapstats = "VRTrack::Mapstats"->create($vrtrack);
ok $autoqc = $mapstats->add_autoqc('NPG QC status',1,'The lane failed the NPG QC check'),'added autoqc to mapping stats';
my @autoqc_statuses = @{ $mapstats->autoqcs() };
is scalar @autoqc_statuses,1, 'got autoqc array ref back from db';

# update autoqc test result
ok $autoqc = $mapstats->add_autoqc('NPG QC status',0,'The lane passed the NPG QC check'),'re-added (updated) an autoqc test';
is $autoqc->result(), 0, 'autoqc result was updated';

# individual
# setting species_id and population id for the first individual
my $ind_name = "test ind";
my $ind_hname = "test_ind";
ok my $individual = VRTrack::Individual->create($vrtrack, $ind_name), 'created individual';
my $ind_id = $individual->id;
$individual = VRTrack::Individual->new($vrtrack, $ind_id);
is $individual->id, $ind_id, 'Individual new returned the individual we made before';
is $individual->name, $ind_name, 'Individual name was stored correctly';
is $individual->hierarchy_name, $ind_hname, 'Individual hierarchy name was generated correctly';
is $individual->sex(), 'unknown', 'Individual sex defaults to unknown if not set';
$ind_name = 'individual_test|changed';
is $individual->name($ind_name), $ind_name, 'Individual name could be get/set';
is $individual->sex('M'), 'M', 'Individual sex could be get/set to M';
is $individual->sex('foobar'), 'unknown', 'Individual sex is set to unknown for non M/F';
is $individual->sex('F'), 'F', 'Individual sex could be get/set to F';
ok $individual->update, 'update worked on Individual';
$individual = VRTrack::Individual->new_by_name($vrtrack, $ind_name);
is_deeply [$individual->id, $individual->name, $individual->sex], [1, $ind_name, 'F'], 'Individual new_by_name worked';
my $species_name = 'genus species subspecies';
is $individual->species($species_name), undef, 'species() returned undef when setting a non-existant species';
is $individual->species_id, undef, 'species_id starts undefined';
ok my $spp = $individual->add_species($species_name), 'add_species made a new species';
ok my $spp_id = $spp->id, 'species id has an id';
is $individual->species_id, $spp_id, 'species_id was updated';
my $pop_name = 'population_test';
is $individual->population($pop_name), undef, 'population() returned undef when setting a non-existant population';
is $individual->population_id, undef, 'population_id starts undefined';
ok my $pop = $individual->add_population($pop_name), 'add_population made a new population';
ok my $pop_id = $pop->id, 'population has an id';

is $individual->population_id, $pop_id, 'population_id was updated';

$individual->update();

# mapper
my $mapper_name = 'test_mapper';
ok my $mapper = VRTrack::Mapper->create($vrtrack, $mapper_name, 1.1), 'Mapper create worked';
my $mapper_id = $mapper->id;
$mapper = VRTrack::Mapper->new($vrtrack, $mapper_id);
is $mapper->id, $mapper_id, 'Mapper new returned the first Mapper we made before';
is $mapper->name, $mapper_name, 'Mapper name set correctly';
is $mapper->version, 1.1, 'Mapper version set correctly';
$mapper_name = 'mapper_test_changed';
is $mapper->name($mapper_name), $mapper_name, 'Mapper name could be get/set';
ok $mapper->update, 'update worked on Mapper';
$mapper = VRTrack::Mapper->new_by_name_version($vrtrack, $mapper_name, 1.1);
is_deeply [$mapper->id, $mapper->name, $mapper->version], [$mapper_id, $mapper_name, 1.1], 'Mapper new_by_name_version worked';

# study
my $study_acc = 'test_study';
ok my $study = VRTrack::Study->create($vrtrack, $study_acc), 'study create worked';
my $study_id = $study->id;
$study = VRTrack::Study->new_by_acc($vrtrack, $study_acc);
is_deeply [$study->id, $study->acc], [$study_id, $study_acc], 'Study new_by_acc worked';
$study_acc .= '2';
is $study->acc($study_acc), $study_acc, 'study acc could be get/set';
ok $study->update, 'update worked on study';
$study = VRTrack::Study->new_by_acc($vrtrack, $study_acc);
is_deeply [$study->id, $study->acc], [$study_id, $study_acc], 'Study new_by_acc worked again';

# species
my $species = VRTrack::Species->create($vrtrack,"organism");
is $species->taxon_id(),0, 'Created a species entry with taxon id 0 as expected';
my $other_species = VRTrack::Species->create($vrtrack,"organism_2",456);
is $other_species->taxon_id(),456, 'Created a species entry with taxon id 456 as expected';

# cleanup
my $dbh = $vrtrack->{_dbh};
foreach ($dbh->tables()){
    next if /`latest_/;
    next if /schema_version/;
    $dbh->do("TRUNCATE TABLE $_");
}
exit;
