#!/usr/bin/perl -w
use strict;
use warnings;
use File::Spec;

BEGIN {
    use Test::Most tests => 26;
    
    use_ok('VertRes::Parser::bas');
}

my $basp = VertRes::Parser::bas->new();
isa_ok $basp, 'VertRes::Parser::ParserI';
isa_ok $basp, 'VertRes::IO';
isa_ok $basp, 'VertRes::Base';

ok my $rh = $basp->result_holder(), 'result_holder returned something';
is ref($rh), 'ARRAY', 'result_holder returns an array';
is @{$rh}, 0, 'the result_holder starts off empty';

ok ! $basp->next_result, 'next_result returns false when we have no file set';

my $bas_file = File::Spec->catfile('t', 'data', 'example2.bas');
ok -e $bas_file, 'file we will test with exists';
ok $basp->file($bas_file), 'file set into parser';

ok $basp->next_result, 'next_result now works';

is_deeply $rh, [qw(NA00001.ILLUMINA.bwa.unknown_population.unknown_analysisgroup.20100208 f203566d66a1751151823a36ff0cfc1d SRP000001 NA00001 ILLUMINA alib SRR00001 115000 62288 2000 1084 1084 1070 2.05 23.32 286 74.10 275 48 2)], 'first result was correct';

while ($basp->next_result) {
    next;
}

is_deeply $rh, [qw(NA00003.ILLUMINA.bwa.unknown_population.unknown_analysisgroup.20100208 f203566d66a1751151823a36ff0cfc1d SRP000003 NA00003 ILLUMINA alib3 SRR00003 115000 62288 2000 1084 1084 1070 2.05 23.32 286 74.10 275 46 2)], 'last result was correct';

# test the total-related methods, and that they work mid-file
$basp = VertRes::Parser::bas->new(file => $bas_file);
$basp->next_result;
$basp->next_result;
$rh = $basp->result_holder();
is_deeply $rh, [qw(NA00002.ILLUMINA.bwa.unknown_population.unknown_analysisgroup.20100208 f203566d66a1751151823a36ff0cfc1d SRP000002 NA00002 ILLUMINA alib2 SRR00002 115000 62288 2000 1084 1084 1070	2.05 23.32 286 74.10 275 47 2)], 'second result was correct';
is $basp->total_reads, 6000, 'total_reads was correct';
is $basp->mapped_reads, 3252, 'mapped_reads was correct';
is sprintf("%0.1f", $basp->percent_mapped), 54.2, 'percent_mapped was correct';
is_deeply $rh, [qw(NA00002.ILLUMINA.bwa.unknown_population.unknown_analysisgroup.20100208 f203566d66a1751151823a36ff0cfc1d SRP000002 NA00002 ILLUMINA alib2 SRR00002 115000 62288 2000 1084 1084 1070	2.05 23.32 286 74.10 275 47 2)], 'second result still correct';
$basp->next_result;
is_deeply $rh, [qw(NA00003.ILLUMINA.bwa.unknown_population.unknown_analysisgroup.20100208 f203566d66a1751151823a36ff0cfc1d SRP000003 NA00003 ILLUMINA alib3 SRR00003 115000 62288 2000 1084 1084 1070 2.05 23.32 286 74.10 275 46 2)], 'last result was correct';

# test parsing an empty (header-only) bas file: should return unknowns and 0s,
# not the header column names
$bas_file = File::Spec->catfile('t', 'data', 'empty.bas');
ok -e $bas_file, 'empty file we will test with exists';
ok $basp->file($bas_file), 'file set into parser';
ok ! $basp->next_result, 'next_result never worked';
$rh = $basp->result_holder();
is_deeply $rh, [qw(unknown unknown unknown unknown unknown unknown unknown 0 0 0 0 0 0 0 0 0 0 0 0 0)], 'result holder contains unknows and 0s';

# test new 20-column bas files
$bas_file = File::Spec->catfile('t', 'data', 'duplicates.bas');
ok -e $bas_file, 'duplicates file we will test with exists';
ok $basp->file($bas_file), 'file set into parser';
$rh = $basp->result_holder();
$basp->next_result;
is $rh->[19], 5, 'number of duplicate reads could be retrieved';

exit;
