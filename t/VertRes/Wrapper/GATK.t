#!/usr/bin/perl -w
use strict;
use warnings;
use File::Copy;
use File::Spec;

BEGIN {
    use Test::Most tests => 12;
    
    use_ok('VertRes::Wrapper::GATK');
    use_ok('VertRes::Utils::FileSystem');
    use_ok('VertRes::Parser::sam');
    use_ok('VertRes::Utils::FastQ');
    use_ok('Math::NumberCruncher');
}

my $gatk = VertRes::Wrapper::GATK->new(quiet => 1);
isa_ok $gatk, 'VertRes::Wrapper::WrapperI';

# setup in/out files
my $fsu = VertRes::Utils::FileSystem->new();
my $ref = File::Spec->catfile('t', 'data', 'S_suis_P17.fa');
my $in_bam = File::Spec->catfile('t', 'data', 'simple.bam');
ok -e $in_bam, 'bam file we will test with exists';
my $temp_dir = $fsu->tempdir();
my $temp_in_bam = File::Spec->catfile($temp_dir, 'in.bam');
copy($in_bam, $temp_in_bam);
ok -s $temp_in_bam, 'copied in bam ready to test with';
$in_bam = $temp_in_bam;
my $out_bam = File::Spec->catfile($temp_dir, 'out.bam');
my $out_csv_prefix = File::Spec->catfile($temp_dir, 'in.bam');
my $out_csv = $out_csv_prefix.'.recal_data.csv';

# individual method tests
$gatk->count_covariates($in_bam, $out_csv_prefix, R => $ref);
ok -s $out_csv, 'count_covariates generated a csv';

$gatk->table_recalibration($in_bam, $out_csv, $out_bam, R => $ref);
ok -s $out_bam, 'table_recalibration generated a bam';
unlink($out_bam);

# the multi-step method, and more careful testing of the output
$gatk = VertRes::Wrapper::GATK->new(quiet => 1, reference => $ref);
$gatk->recalibrate($in_bam, $out_bam);
ok -s $out_bam, 'recalibrate generated a bam';
my @recal_qs = get_qualities($out_bam);
my $recal_mean_q = Math::NumberCruncher::Mean(\@recal_qs);
# (orig_mean_q == 20.8908434782609)
ok $recal_mean_q > 10 && $recal_mean_q < 20, 'recalibrated qualities changed, and are in the right ball-park';

exit;


sub get_qualities {
    my $sp = VertRes::Parser::sam->new(file => shift);
    my $fqu = VertRes::Utils::FastQ->new();
    
    my @qs;
    while (my ($qual_str) = $sp->get_fields('QUAL')) {
        push(@qs, $fqu->qual_to_ints($qual_str));
    }
    
    return @qs;
}
