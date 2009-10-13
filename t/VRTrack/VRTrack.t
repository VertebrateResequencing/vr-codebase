#!/usr/bin/perl -w
use strict;
use warnings;

BEGIN {
    use Test::Most; # tests => 56;

    use_ok('VRTrack::VRTrack');
}

my $db =
{
    database => 'vrtrack_test',
    host     => 'mcs4a',
    port     => 3306,
    user     => 'vreseq_rw',
    password => 't3aml3ss',
};

my $vrtrack = VRTrack::VRTrack->new($db);
isa_ok $vrtrack, 'VRTrack::VRTrack';

my $name = '1111_1';
my $vrlane  = VRTrack::Lane->new_by_name($vrtrack,$name);
if ( !$vrlane )
{
    $vrlane = VRTrack::Lane->create($vrtrack,$name);
    $vrlane->update();
}

$name = "${name}_1.fastq";
my $vrfile = 0; #$vrlane->get_file_by_name($name);
if ( !$vrfile )
{
    $vrfile = $vrlane->add_file($name); # This will add two rows
}

my $mapping = $vrlane->add_mapping();   # This will add two rows
$mapping->raw_reads(999);
$mapping->raw_bases(111);
$mapping->update;                       # This will add another row


