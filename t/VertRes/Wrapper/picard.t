#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;

BEGIN {
    use Test::Most tests => 25;
    
    use_ok('VertRes::Wrapper::picard');
    use_ok('VertRes::Wrapper::samtools');
    use_ok('VertRes::Parser::sam');
    use_ok('VertRes::Utils::FileSystem');
}

my $pt = VertRes::Wrapper::picard->new(quiet => 1, java_memory => 500);
isa_ok $pt, 'VertRes::Wrapper::WrapperI';
is $pt->quiet, 1, 'quiet set via new';

# setup files
my $fsu = VertRes::Utils::FileSystem->new();
my $bam_input_file = File::Spec->catfile('t', 'data', 'simple.bam');
ok -s $bam_input_file, 'input bam file ready to test on';
my $temp_dir = $fsu->tempdir();
my $bam_out_file = File::Spec->catfile($temp_dir, 'merged.bam');
my $rmdup_file = File::Spec->catfile($temp_dir, 'rmdup.bam');
my $markdup_file = File::Spec->catfile($temp_dir, 'markdup.bam');

# MergeSamFiles
$pt->MergeSamFiles($bam_out_file, ($bam_input_file, $bam_input_file),
                   VALIDATION_STRINGENCY => 'SILENT',
                   TMP_DIR => $temp_dir);
cmp_ok $pt->run_status, '>=', 1, 'MergeSamFiles ran ok';
is get_bam_lines($bam_out_file), 4000, 'merged bam had the correct number of lines';

# we can also set stringency and tmp via new
unlink($bam_out_file);
$pt = VertRes::Wrapper::picard->new(quiet => 2,
                                    java_memory => 500,
                                    validation_stringency => 'silent',
                                    tmp_dir => $temp_dir);
$pt->MergeSamFiles($bam_out_file, ($bam_input_file, $bam_input_file));
cmp_ok $pt->run_status, '>=', 1, 'MergeSamFiles ran ok with settings via new';
is get_bam_lines($bam_out_file), 4000, 'merged bam had the correct number of lines with settings via new';

# and do this same thing but with checks:
unlink($bam_out_file);
$pt->merge_and_check($bam_out_file, [$bam_input_file, $bam_input_file]);
cmp_ok $pt->run_status, '>=', 1, 'merge_and_check ran ok';
is get_bam_lines($bam_out_file), 4000, 'merged bam had the correct number of lines with merge_and_check';

# test of rmdup should bring that back down to 2000... actually, 3972 (!).
# samtools rmdup gets rid of more, but still only gets it to 2927. Both samtools
# and picard agree that simple.bam contains 6 duplicates. Whatever...
$pt->rmdup($bam_out_file, $rmdup_file);
is $pt->run_status, 2, 'rmdup ran ok';
is get_bam_lines($rmdup_file), 3972, 'rmdup bam had the correct number of lines';

$pt->markdup($bam_out_file, $markdup_file);
is $pt->run_status, 2, 'markdup ran ok';
is get_bam_lines($markdup_file), 4000, 'markdup bam had the correct number of lines';
my $sp = VertRes::Parser::sam->new(file => $markdup_file);
my ($normals, $dups) = (0, 0);
while (my ($flag) = $sp->get_fields('FLAG')) {
    if ($sp->is_duplicate($flag)) {
        $dups++;
    }
    else {
        $normals++;
    }
}
is_deeply [$normals, $dups], [3972, 28], 'markdup actually marked the correct number of lines as dups';

# merge of two headed bams with the same PG header line shouldn't rename the
# id, but it does!
my $headed1_bam = File::Spec->catfile('t', 'data', 'headed1.bam');
ok -s $headed1_bam, 'headed1 bam file ready to test on';
my $headed2_bam = File::Spec->catfile('t', 'data', 'headed2.bam');
ok -s $headed2_bam, 'headed2 bam file ready to test on';
unlink($bam_out_file);
$pt->merge_and_check($bam_out_file, [$headed1_bam, $headed2_bam]);
cmp_ok $pt->run_status, '>=', 1, 'merge_and_check ran ok on 2 headed bams';
is get_bam_lines($bam_out_file), 4000, 'merged bam had the correct number of lines with merge_and_check on 2 headed bams';
my $st = VertRes::Wrapper::samtools->new(quiet => 1, run_method => 'open');
my $fh = $st->view($headed1_bam, undef, h => 1);
$sp = VertRes::Parser::sam->new(fh => $fh);
is $sp->program, 'bwa', 'input bam1 has expected PG ID';
$fh = $st->view($headed2_bam, undef, h => 1);
$sp->fh($fh);
is $sp->program, 'bwa', 'input bam2 has expected PG ID';
$fh = $st->view($bam_out_file, undef, h => 1);
$sp->fh($fh);
is $sp->program, 'bwa', 'merge didn\'t change the PG ID';

exit;

# verify the merged bam is ok, by converting it to sam with samtools
sub get_bam_lines {
    my $bam = shift;
    
    my $st = VertRes::Wrapper::samtools->new(quiet => 1, run_method => 'open');
    my $fh = $st->view($bam, undef);
    
    my $lines = 0;
    while (<$fh>) {
        $lines++;
    }
    close($fh);
    
    return $lines;
}
