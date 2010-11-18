#!/usr/bin/perl -w
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 56;

    use_ok('VRTrack::VRTrack');
    use_ok('VRTrack::Library');
    use_ok('VRTrack::Multiplex_pool');
    use_ok('VRTrack::Library_Multiplex_pool');
}
use VRTrack::Testconfig;

my $connection_details = { database => VRTrack::Testconfig->config('test_db'),
                           host     => VRTrack::Testconfig->config('host'),
                           port     => VRTrack::Testconfig->config('port'),
                           user     => VRTrack::Testconfig->config('user'),
                           password => VRTrack::Testconfig->config('password'), 
                        };

# access the db as normal with new()
ok my $vrtrack = VRTrack::VRTrack->new($connection_details), 'new() returned something';
isa_ok $vrtrack, 'VRTrack::VRTrack';

# Need to set up project & sample so we can test library request.
# project
my $project_name = 'Project reqtest';
my $vrproj = $vrtrack->add_project('Project reqtest');
my $sample = $vrproj->add_sample('sample_reqsample');
my $sample_id = $sample->id;
$vrproj->update;

# LIBRARY REQUEST
my $libreq_ssid = 1234;
ok my $library_request=$sample->add_library_request($libreq_ssid),'add library_request to sample worked';
is $library_request->sample_id(),$sample_id,"library request linked to sample $sample_id correctly";
is $library_request->ssid(),$libreq_ssid,'library request ssid correct';
$libreq_ssid++;
is $library_request->ssid($libreq_ssid),$libreq_ssid,'library request set ssid worked';
is $library_request->prep_status('passed'),'passed','library request get/set prep_status worked';
$library_request->update;

{
    my $tmp_lr;
    warning_like {$tmp_lr = $sample->add_library_request($libreq_ssid)} {carped => qr/already present in the database/i}, 'add library_request with existing ssid throws warning';
    is $tmp_lr,undef,'add library_request with existing ssid correctly failed';
}

# check retrieval
{
    ok my $library_request1=$sample->get_library_request_by_ssid($libreq_ssid),'retrieve library_request by ssid worked';
    is $library_request1->id, $library_request->id, 'library_request is the correct library request';
    ok my $library_request2=$sample->get_library_request_by_id($library_request->id),'retrieve library_request by id worked';
}

# LIBRARY

my $libname = 'lib_from_req_a';
my $library = $library_request->add_library($libname);
isa_ok $library, 'VRTrack::Library', 'add library to request returned a library'; 
is $library->sample_id, $sample->id, 'library linked to sample correctly';

# check retrieval
{
    ok my $lib = VRTrack::Library->new_by_name($vrtrack, $libname), 'library new_by_name worked';
    is $lib->id, $library->id, 'library is the correct library';
}

# link to sample?
{
    my $tmp_sample = $vrproj->get_sample_by_name($sample->name);
    is @{$tmp_sample->libraries}, 1, 'sample returned the correct number of libraries';
    cmp_bag ($tmp_sample->library_ids, [$library->id], 'sample is linked to library correctly');
}

# SEQ_REQUEST

my $seqreq_ssid = 1111;
my $seqreq2_ssid = 2222;
my $seq_request=$library->add_seq_request($seqreq_ssid);
isa_ok $seq_request, 'VRTrack::Seq_request','add seq_request to library returned a Seq_request'; 
is $seq_request->ssid, $seqreq_ssid, 'SSID set correctly on seq_request';

my $seq_request2=$library->add_seq_request($seqreq2_ssid);
isa_ok $seq_request2, 'VRTrack::Seq_request','add seq_request to library returned a Seq_request again'; 
is $seq_request2->ssid, $seqreq2_ssid, 'SSID set correctly on seq_request';

isnt $seq_request->id, $seq_request2->id, 'seq requests are different';

$seqreq_ssid++;
is $seq_request->ssid($seqreq_ssid),$seqreq_ssid,'seq request set ssid worked';
is $seq_request->seq_status('passed'),'passed','seq request get/set seq_status worked';
is $seq_request->seq_type('Single ended sequencing'),'Single ended sequencing','seq request get/set seq_type worked';
$seq_request->update;
is @{$library->seq_requests}, 2, 'seq_requests returned correct number of requests';
cmp_bag ($library->seq_request_ids, [$seq_request->id, $seq_request2->id], 'seq_request_ids returned the correct ids');

# check retrieval
{
   ok my $tmp_req = $library->get_seq_request_by_ssid($seqreq_ssid),'retrieve seq_request by ssid worked';
   is $tmp_req->id, $seq_request->id, 'seq_request is the correct seq_request';
   ok my $tmp_req2 = $library->get_seq_request_by_id($seq_request->id),'retrieve seq_request by id worked';
   is $tmp_req2->ssid, $seq_request->ssid, 'seq_request is the correct seq_request again';
}

# LANE

my $lane_1 = $seq_request->add_lane('lane_from_req_a');
isa_ok $lane_1, 'VRTrack::Lane', 'add_lane returned a lane';

my $lane_2 = $seq_request->add_lane('lane_from_req_b');
isa_ok $lane_2, 'VRTrack::Lane', 'add_lane returned a lane again';

$seq_request->update;
isnt $lane_1->id, $lane_2->id, 'lanes are different';

is $lane_1->name, 'lane_from_req_a', 'lane name set correctly';

is @{$library->lanes}, 2, 'library->lanes returned the correct number of lanes';
cmp_bag( $library->lane_ids, [$lane_1->id, $lane_2->id], 'library returns lane_ids correctly');

is $lane_1->library_id, $library->id, 'lane correctly linked to library';

# check retrieval by lib
{
    ok my $tmp_lane = $library->get_lane_by_name($lane_1->name),'get_lane_by_name worked';
    is $tmp_lane->id, $lane_1->id, 'lane is the correct lane';
    ok my $tmp_lane2 = $library->get_lane_by_id($lane_1->id), 'get_lane_by_id worked';
    is $tmp_lane2->id, $lane_1->id, 'lane is the correct lane again';
}

# check retrieval by seq_req
{
    my $tmp_req = $library->get_seq_request_by_id($seq_request->id);
    cmp_bag ($tmp_req->lane_ids, [$lane_1->id, $lane_2->id], 'lanes are linked to seq_req correctly');
}

# Multiplex Seq request
# setup multiplexing
my $mplex_pool=VRTrack::Multiplex_pool->new_by_ssid($vrtrack,1);
unless ($mplex_pool){
    $mplex_pool=VRTrack::Multiplex_pool->create($vrtrack,1);
}

# library_multiplex_pool
my $lib_mplex_pool=VRTrack::Library_Multiplex_pool->new_by_library_id_multiplex_pool_id($vrtrack,$library->id,$mplex_pool->id);
unless ($lib_mplex_pool){
    $lib_mplex_pool=VRTrack::Library_Multiplex_pool->create($vrtrack,$library->id,$mplex_pool->id);
}
                    
my $mplex_seqreq_id = 3333;
ok my $mplex_seq_request=$mplex_pool->add_seq_request($mplex_seqreq_id),'add seq_request to multiplex_pool worked';
isa_ok $mplex_seq_request, 'VRTrack::Seq_request', 'add_multiplex_seq_request returns a Seq_request';
is $mplex_seq_request->ssid, $mplex_seqreq_id, 'ssid set correctly on seq req';
is $mplex_seq_request->multiplex_pool_id, $mplex_pool->id, 'seq req linked correctly to mplex_pool';

$mplex_seqreq_id++;
is $mplex_seq_request->ssid($mplex_seqreq_id),$mplex_seqreq_id,'multiplex seq request set ssid worked';
is $mplex_seq_request->seq_status('passed'),'passed','multiplex seq request get/set seq_status worked';
is $mplex_seq_request->seq_type('Single ended sequencing'),'Single ended sequencing','seq request get/set seq_type worked';
$mplex_seq_request->update;

# Need some more tests on adding lanes to mplex_seqreq, but they won't be
# added to the library correctly, as this has to be done manually after
# deplexing.  Not sure how to do/test this.

# delete
ok $vrproj->delete, 'deleting project works';

# cleanup
my $dbh = $vrtrack->{_dbh};
foreach ($dbh->tables()){
    next if /`latest_/;
    next if /schema_version/;
    $dbh->do("TRUNCATE TABLE $_");
}
exit;
