#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;


BEGIN { unshift(@INC, './modules') }
BEGIN {
    use Test::Most;
    use_ok('VertRes::Pipelines::RNASeqExpression');
}

my $obj = VertRes::Pipelines::RNASeqExpression->new();

is($obj->does_reference_in_bam_match_annotation("t/data/rnaseq_expression.bam", "/lustre/scratch108/pathogen/pathpipe/refs/Mus/musculus/Mus_musculus_mm10.gff"), 1, 'annotation base name matches reference');
is($obj->does_reference_in_bam_match_annotation("t/data/rnaseq_expression.bam", "Mus_musculus_mm10.gff"), 1, 'annotation base name matches reference where theres no directory');
is($obj->does_reference_in_bam_match_annotation("t/data/rnaseq_expression.bam", "/lustre/scratch108/pathogen/pathpipe/refs/Mus/musculus/Mus_musculus_mm10.fa.gff"), 1, 'annotation base name matches reference where the extension is still on');
is($obj->does_reference_in_bam_match_annotation("t/data/rnaseq_expression.bam", "/path/to/something_which_doesnt_match.gff"), 0, 'annotation base name doesnt matche reference');

done_testing();
