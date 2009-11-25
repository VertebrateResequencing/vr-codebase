#!/usr/bin/perl -w
use strict;
use warnings;
use File::Spec;

BEGIN {
    use Test::Most tests => 27;
    
    use_ok('VertRes::Parser::flagstat');
}

my $fsp = VertRes::Parser::flagstat->new();
isa_ok $fsp, 'VertRes::Parser::ParserI';
isa_ok $fsp, 'VertRes::IO';
isa_ok $fsp, 'VertRes::Base';

ok my $rh = $fsp->result_holder(), 'result_holder returned something';
is ref($rh), 'ARRAY', 'result_holder returns an array';
is @{$rh}, 0, 'the result_holder starts off empty';

ok ! $fsp->next_result, 'next_result returns false when we have no file set';
is $fsp->total_reads(), undef, 'total_reads returns undef before we set file';

my $fs_file = File::Spec->catfile('t', 'data', 'bam.flagstat');
ok -e $fs_file, 'file we will test with exists';
ok $fsp->file($fs_file), 'file set into parser';

ok ! $fsp->next_result, 'next_result returns false even after setting the file';

is $fsp->total_reads(), 2000, 'total_reads test';
is $fsp->qc_failures(), 0, 'qc_failures test';
is $fsp->duplicates(), 0, 'duplicates test';
is $fsp->mapped_reads(), 50, 'mapped_reads test';
is $fsp->paired_reads(), 2000, 'paired_reads test';
is $fsp->read1_reads(), 1000, 'read1_reads test';
is $fsp->read2_reads(), 1000, 'read2_reads test';
is $fsp->mapped_proper_paired_reads(), 26, 'mapped_proper_paired_reads test';
is $fsp->mapped_paired_reads(), 28, 'mapped_paired_reads test';
is $fsp->singletons(), 1950, 'singletons test';
is $fsp->mate_mapped_to_other_chr(), 2, 'mate_mapped_to_other_chr test';
is $fsp->mate_mapped_to_other_chr(1), 1, 'mate_mapped_to_other_chr(1) test';

# check we can change files
$fs_file = File::Spec->catfile('t', 'data', 'bam2.flagstat');
ok -e $fs_file, 'second file we will test with exists';
ok $fsp->file($fs_file), 'second file set into parser';

is $fsp->mapped_proper_paired_reads(), 1696, 'correct result after changing file';

exit;
