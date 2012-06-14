#!/usr/bin/perl -w
use strict;
use warnings;

BEGIN {
    use Test::Most;
    use VRTrack::VRTrack;
    eval {
        require VRTrack::Testconfig;
    };
    if ($@) {
        plan skip_all => "Skipping all tests because VRTrack tests have not been configured";
    }
    use_ok('VRTrack::Lane');
}

my $connection_details = { database => VRTrack::Testconfig->config('test_db'),
                           host     => VRTrack::Testconfig->config('host'),
                           port     => VRTrack::Testconfig->config('port'),
                           user     => VRTrack::Testconfig->config('user'),
                           password => VRTrack::Testconfig->config('password') };
my $vrtrack = VRTrack::VRTrack->new($connection_details);
clean_test_db($vrtrack);

my $vr_lane = VRTrack::Lane->create( $vrtrack, "1234_5");
my $mapping_1 = $vr_lane->add_mapping();

is $mapping_1->is_qc, 0, 'not qc by default';
is $mapping_1->prefix, '_', 'defaults to normal prefix';
# populate with some test data
my $mapping_qc = $vr_lane->add_mapping();
$mapping_qc->is_qc(1);
$mapping_qc->update();
push(@expected_qc_mappings, $mapping_qc);
# add some other mappings
my @expected_mappings;
my @expected_qc_mappings;
push(@expected_mappings,$mapping_1 );
push(@expected_mappings,$vr_lane->add_mapping());
push(@expected_mappings,$vr_lane->add_mapping());
# a final QC mapping
my $mapping_qc_2 = $vr_lane->add_mapping();
$mapping_qc_2->is_qc(1);
$mapping_qc_2->update();
push(@expected_qc_mappings, $mapping_qc_2);

is_deeply $vr_lane->mappings_excluding_qc, \@expected_mappings, 'only non qc mappings returned';
is_deeply $vr_lane->qc_mappings, \@expected_qc_mappings, 'only qc mappings returned';

clean_test_db($vrtrack);
done_testing();

sub clean_test_db
{
    my ($vrtrack) = @_;

    my $dbh = $vrtrack->{_dbh};
    foreach ($dbh->tables()){
    next if /\`latest_/;
    next if /schema_version/;
    $dbh->do("TRUNCATE TABLE $_");
    }
}
