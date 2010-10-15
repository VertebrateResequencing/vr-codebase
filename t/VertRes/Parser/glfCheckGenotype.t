#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;

BEGIN {
    use Test::Most tests => 46;
    
    use_ok('VertRes::Parser::glfCheckGenotype');
}

my $gcg = VertRes::Parser::glfCheckGenotype->new();
isa_ok $gcg, 'VertRes::Parser::ParserI';
isa_ok $gcg, 'VertRes::IO';
isa_ok $gcg, 'VertRes::Base';

ok my $rh = $gcg->result_holder(), 'result_holder returned something';
is ref($rh), 'ARRAY', 'result_holder returns an array';
is @{$rh}, 0, 'the result_holder starts off empty';

ok ! $gcg->next_result, 'next_result returns false when we have no file set';

my $check_genotype_file = File::Spec->catfile('t', 'data', 'glf_checkgenotype.out');
ok -e $check_genotype_file, 'file we will test with exists';
ok $gcg->file($check_genotype_file), 'file set into parser';

my @expected_data = ([qw(snps-hapmap3/NA19468.snp 341694 76226 0.101695)],
                     [qw(snps-hapmap3/NA19046.snp 341762 74898 0.099907)],
                     [qw(snps-hapmap3/NA20589.snp 443997 76985 0.102706)],
                     [qw(snps-hapmap3/NA20810.snp 781899 79255 0.105717)]);

while ($gcg->next_result) {
    my $expected = shift @expected_data;
    foreach my $i (0..3) {
        is $rh->[$i], $expected->[$i], 'parsed data test';
    }
}
is @expected_data, 0, 'parsed out all the results';

is $gcg->entropy, 32.8, 'entropy test';

# can we change file and still parse OK?
$check_genotype_file = File::Spec->catfile('t', 'data', 'glf_checkgenotype2.out');
$gcg->file($check_genotype_file);

is $gcg->entropy, 30.2, 'entropy test';

@expected_data = ([qw(snps-hapmap3/NA19468.snp 5 76226 0.101695)],
                  [qw(snps-hapmap3/NA19046.snp 10 74898 0.099907)],
                  [qw(snps-hapmap3/NA20589.snp 15 76985 0.102706)],
                  [qw(snps-hapmap3/NA20810.snp 20 79255 0.105717)]);

while ($gcg->next_result) {
    my $expected = shift @expected_data;
    foreach my $i (0..3) {
        is $rh->[$i], $expected->[$i], 'parsed data test';
    }
}
is @expected_data, 0, 'parsed out all the results';

exit;
