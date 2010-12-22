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
        plan tests => 42;
    }
    
    use_ok('VRTrack::VRTrack');
    use_ok('VRTrack::Multiplex_pool');
    use_ok('VRTrack::Library');
    use_ok('VRTrack::Library_Multiplex_pool');
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

# setup libs to go into multiplex
my $lib_a = VRTrack::Library->create($vrtrack,'lib_a');
my $lib_b = VRTrack::Library->create($vrtrack,'lib_b');

# multiplex pool
my $mplex_ssid = 9999;
my $multiplex_pool=VRTrack::Multiplex_pool->new_by_ssid($vrtrack,$mplex_ssid);
is $multiplex_pool,undef,'multiplex_pool new_by_ssid on unknown ssid worked';
ok $multiplex_pool=VRTrack::Multiplex_pool->create($vrtrack,$mplex_ssid),'create a new multiplex_pool by ssid worked';
ok my $mplex_id = $multiplex_pool->id, 'pool has an id';
is $multiplex_pool->ssid,$mplex_ssid, 'mplex_pool ssid set correctly';

{
    my $tmp_mplex;
    warning_like {$tmp_mplex = VRTrack::Multiplex_pool->create($vrtrack,$mplex_ssid)} {carped => qr/already present in the database/i}, 'Create mplex_pool with existing ssid throws warning';
    is $tmp_mplex,undef,'Create mplex_pool with existing ssid correctly failed';
}

my $mplex_name= "My multiplex";
is $multiplex_pool->name($mplex_name), $mplex_name, 'mplex_pool name set works';
ok $multiplex_pool->dirty, 'mplex_pool gets dirty when name changed';
ok $multiplex_pool->update, 'mplex_pool update works';

#check retrieval
{
    ok my $tmp_mplex = VRTrack::Multiplex_pool->new_by_ssid($vrtrack,$mplex_ssid),'mplex_pool new_by_ssid works correctly';
    is $tmp_mplex->id,$multiplex_pool->id,'mplex_pool is the correct mplex_pool';
    
    ok my $tmp_mplex2 = VRTrack::Multiplex_pool->new_by_name($vrtrack,$mplex_name),'mplex_pool new_by_name works correctly';
    is $tmp_mplex2->id,$multiplex_pool->id,'mplex_pool is the correct mplex_pool';

}

$mplex_ssid++;
is $multiplex_pool->ssid($mplex_ssid),$mplex_ssid, 'mplex_pool ssid get/set works';

ok $multiplex_pool->dirty, 'mplex_pool gets dirty when ssid changed';
ok $multiplex_pool->update, 'mplex_pool update works again';

$multiplex_pool=VRTrack::Multiplex_pool->new_by_ssid($vrtrack,$mplex_ssid);
is $multiplex_pool->id,$mplex_id,'multiplex_pool new_by_ssid on known ssid worked';

# library_multiplex_pool
my $lib_multiplex_pool=VRTrack::Library_Multiplex_pool->new_by_library_id_multiplex_pool_id($vrtrack,$lib_a->id,$multiplex_pool->id);
is $lib_multiplex_pool,undef,'library_multiplex_pool new_by_library_id_multiplex_pool_id on unknown library_id,multiplex_pool_id';

ok $lib_multiplex_pool=VRTrack::Library_Multiplex_pool->create($vrtrack,$lib_a->id,$multiplex_pool->id),'create a library_multiplex_pool worked';
ok my $lib_mpp_id = $lib_multiplex_pool->id, 'lib mplex has an id';
is $lib_multiplex_pool->multiplex_pool_id,$mplex_id,'library_multiplex_pool get multiplex_pool_id works';
is $lib_multiplex_pool->library_id,$lib_a->id,'library_multiplex_pool get library_id works';

$lib_multiplex_pool=VRTrack::Library_Multiplex_pool->new_by_library_id_multiplex_pool_id($vrtrack,$lib_a->id,$multiplex_pool->id);
is $lib_multiplex_pool->id,$lib_mpp_id,'library_multiplex_pool retrieval works';

ok my $lib_multiplex_pool2=VRTrack::Library_Multiplex_pool->create($vrtrack,$lib_b->id,$multiplex_pool->id),'can add two libraries to a multiplex_pool';

# check relationships
ok my $lib_mps = $multiplex_pool->library_multiplex_pools(), 'can retrieve library links off multiplex_pool';
my @lib_ids = map {$_->library_id} @$lib_mps;
cmp_bag(\@lib_ids, [$lib_a->id, $lib_b->id],'multiplex linked to libraries correctly');

ok $lib_mps = $lib_a->library_multiplex_pools(), 'can retrieve mplex links off library';
my @mplex_ids = map {$_->multiplex_pool_id} @$lib_mps;
cmp_bag(\@mplex_ids, [$multiplex_pool->id],'library linked to multiplex correctly');

# multiplex seq request
ok my $seq_request=$multiplex_pool->add_seq_request('12345'),'add seq_request to multiplex_pool worked';
is $seq_request->ssid(),12345,'seq request has correct ssid';
is $seq_request->multiplex_pool_id(),$multiplex_pool->id,'seq request has correct mplex_pool id';
is $seq_request->library_id(),undef,'seq request has correctly not got a library id';
is $seq_request->ssid(4321),4321,'multiplex seq request set ssid worked';
is $seq_request->seq_status('passed'),'passed','multiplex seq request get/set seq_status worked';
is $seq_request->seq_type('Single ended sequencing'),'Single ended sequencing','seq request get/set seq_type worked';
ok $seq_request->update, 'mplex_seq_req updated';    
    
# cleanup
my $dbh = $vrtrack->{_dbh};
foreach ($dbh->tables()){
    next if /`latest_/;
    next if /schema_version/;
    $dbh->do("TRUNCATE TABLE $_");
}

exit;
