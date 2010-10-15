#!/usr/bin/env perl
use strict;
use warnings;
use Cwd 'abs_path';
use File::Spec;

BEGIN {
    use Test::Most tests => 16;
    
    use_ok('VertRes::Utils::FastQ');
    use_ok('VertRes::Utils::FileSystem');
    use_ok('VertRes::IO');
}

my $fastq_util = VertRes::Utils::FastQ->new(verbose => 1);
isa_ok $fastq_util, 'VertRes::Base';

# setup our input fastqs
my $io = VertRes::IO->new();
my $fq1_file = File::Spec->catfile('t', 'data', '2822_6_1_1000.fastq');
ok -s $fq1_file, 'first fastq file ready to test with';
# (it has 61bp reads and 1000 sequences)
my $fq2_file = File::Spec->catfile('t', 'data', '2822_6_2_1000.fastq');
ok -s $fq2_file, 'first fastq file ready to test with';
# (it has 54bp reads and 1000 sequences)

# test split
my $fsu = VertRes::Utils::FileSystem->new();
my $split_dir = $fsu->tempdir();
is $fastq_util->split([$fq1_file], chunk_size => 6100, split_dir => $split_dir), 10, 'split with one fastq worked';
my @expected = (400, 400, 400, 400, 400, 400, 400, 400, 400, 400);
my @got = split_check($split_dir);
is_deeply \@got, \@expected, 'each split had the correct number of lines';

$split_dir = $fsu->tempdir();
is $fastq_util->split([$fq1_file, $fq2_file], chunk_size => 11500, split_dir => $split_dir), 10, 'split with two fastqs worked';
push(@expected, @expected);
@got = split_check($split_dir);
is_deeply \@got, \@expected, 'each split had the correct number of lines, and there were twice as many files';

# when just 1 split we only make symlinks
$split_dir = $fsu->tempdir();
is $fastq_util->split([$fq1_file], chunk_size => 61000, split_dir => $split_dir), 1, 'split with one fastq -> 1 split worked';
my $split_file = File::Spec->catfile($split_dir, '2822_6_1_1000.1.fastq');
ok -l $split_file, 'split file was made and it was a symlink';
is abs_path($split_file), abs_path($fq1_file), 'symlink points to the original fastq';

# test qual_to_ints
is_deeply [$fastq_util->qual_to_ints('!"#$%&\'()*+,5?DIS]~')], [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 20, 30, 35, 40, 50, 60, 93], 'qual_to_ints test';

# filter_reads and clip_point
TODO: {
    local $TODO = "filter_reads and clip_point need tests...";
    ok 0;
    ok 0;
}

exit;

sub split_check {
    my $split_dir = shift;
    
    opendir(my $dirfh, $split_dir);
    
    my @got = ();
    foreach my $file (readdir($dirfh)) {
        if ($file =~ /fastq.gz$/) {
            my $split_file = File::Spec->catfile($split_dir, $file);
            $io->file($split_file);
            push(@got, $io->num_lines);
        }
    }
    
    return @got;
}
