#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd 'cwd';

BEGIN {
    use Test::Most tests => 40;
    
    use_ok('VertRes::Utils::FileSystem');
}

my $fsu = VertRes::Utils::FileSystem->new();
isa_ok $fsu, 'VertRes::Base';
isa_ok $fsu, 'VertRes::Utils::FileSystem';

ok -d $fsu->catfile('t', 'data'), 'catfile on two dirs works';

# get a temp file
my ($fh, $file) = $fsu->tempfile;
close($fh);
# write to it

# get a temp dir, test rmtree and get_filepaths as well
my $tmp_dir = $fsu->tempdir;
ok -d $tmp_dir, 'tmpdir created ok';
my $test_dir = $fsu->catfile($tmp_dir, 'test_dir');
mkdir($test_dir);
ok -d $test_dir, 'subdir created in tempdir';
my $foo_file = $fsu->catfile($test_dir, 'foo.txt');
system("touch $foo_file");
ok -e $foo_file, 'file created in subdir of tempdir';
my $bar_file = $fsu->catfile($test_dir, 'bar.gif');
system("touch $bar_file");
my $gz_file = $fsu->catfile($test_dir, 'llama.ps.gz');
system("touch $gz_file");
my $dot_file = $fsu->catfile($test_dir, '.dot');
system("touch $dot_file");
is_deeply [$fsu->get_filepaths($tmp_dir)], [$foo_file, $bar_file, $gz_file, $dot_file], 'get_filepaths no extra args test';
is_deeply [$fsu->get_filepaths($tmp_dir, suffix => 'gif')], [$bar_file], 'get_filepaths suffix test';
is_deeply [$fsu->get_filepaths($tmp_dir, suffix => 'ps.gz')], [$gz_file], 'get_filepaths suffix with .gz test';
is_deeply [$fsu->get_filepaths($tmp_dir, prefix => 'foo')], [$foo_file], 'get_filepaths prefix test';
is_deeply [$fsu->get_filepaths($tmp_dir, filename => 'f.+xt')], [$foo_file], 'get_filepaths filename test';
is_deeply [$fsu->get_filepaths($tmp_dir, subdir => 'test')], [$foo_file, $bar_file, $gz_file, $dot_file], 'get_filepaths subdir test';
is_deeply [$fsu->get_filepaths($tmp_dir, dir => 'test_dir')], [$test_dir], 'get_filepaths dir test';
is_deeply [$fsu->get_filepaths($tmp_dir, subdir => 'test', dir => 'test_dir')], [$test_dir], 'get_filepaths dir + subdir test';
is_deeply [$fsu->get_filepaths($tmp_dir, subdir => 'moo', dir => 'test_dir')], [], 'get_filepaths dir + bad subdir test';
is_deeply [$fsu->get_filepaths($tmp_dir, dir => 'moo')], [], 'get_filepaths bad dir';

# copy... impossible to test the diff part of copy failure?
$file = $fsu->catfile('t', 'data', 'io_test.txt');
my $copy = $fsu->catfile($tmp_dir, 'copy_test');
ok $fsu->copy($file, $copy), 'simple copy test';
ok -s $copy, 'copy test really did work';
my $devnull = $fsu->catfile('dev', 'null', 'copy_test');
my $ok = $fsu->copy($file, $devnull);
is $ok, 0, 'copy fails to /dev/null';

my $sub_dir = $fsu->catfile($test_dir, 'sub_dir');
mkdir($sub_dir);
my $sub_file = $fsu->catfile($sub_dir, 'sub_file');
system("touch $sub_file");
my $copy_dir = $fsu->catfile($tmp_dir, 'copy_dir');
ok $fsu->copy($test_dir, $copy_dir), 'copy on a directory';
my $all_present = 1;
foreach my $file ($foo_file, $bar_file, $gz_file, $dot_file, $sub_file) {
    my $copy_file = $file;
    $copy_file =~ s/test_dir/copy_dir/;
    unless (-e $copy_file) {
        $all_present = 0;
        last;
    }
}
ok $all_present, 'copy on a directory really did work';

ok $fsu->directory_structure_same($test_dir, $copy_dir), 'directory_structure_same positive test';
ok $fsu->directory_structure_same($test_dir, $copy_dir, consider_files => 1), 'directory_structure_same positive test (considering files)';
sleep(2);
system("touch $sub_file");
ok ! $fsu->directory_structure_same($test_dir, $copy_dir, consider_files => 1), 'directory_structure_same negative test (considering files)';
my %leaf_mtimes;
$fsu->directory_structure_same($test_dir, $copy_dir, leaf_mtimes => \%leaf_mtimes);
my @sds = stat($sub_dir);
is_deeply \%leaf_mtimes, {$test_dir => { sub_dir => $sds[9] }}, 'directory_structure_same leaf_mtimes gave the right hash structure';
ok $fsu->rmtree($sub_dir), 'rmtree removed the test sub dir';
ok ! $fsu->directory_structure_same($test_dir, $copy_dir), 'directory_structure_same negative test (missing directory)';

$fsu->rmtree($test_dir);
ok ! -d $test_dir, 'rmtree removed a directory that contained files';

ok $fsu->move($copy_dir, $test_dir), 'move on a directory';
ok ! -d $copy_dir, 'source dir has gone';
$all_present = 1;
foreach my $file ($foo_file, $bar_file, $gz_file, $dot_file, $sub_file) {
    unless (-e $file) {
        $all_present = 0;
        last;
    }
}
ok $all_present, 'move on a directory really did work';
my $move = $copy.'.move';
ok $fsu->move($copy, $move), 'move on a file';
ok ! -e $copy, 'move on a file deleted the source';
ok -e $move, 'move on a file created the destination';

undef $fsu;
ok ! -d $tmp_dir, 'tmpdir destroyed ok';

$fsu = VertRes::Utils::FileSystem->new();

# md5 tests
my $test_file = $fsu->catfile(qw(t data bad_1.fastq.gz));
ok -s $test_file, "test file present for md5 tests";
is $fsu->calculate_md5($test_file), '6fd5c2f7f105b15d44ed374916d0cd34', 'calculate_md5 worked';
ok $fsu->verify_md5($test_file, '6fd5c2f7f105b15d44ed374916d0cd34'), 'verify_md5 worked in a positive test';
ok ! $fsu->verify_md5($test_file, '6fd5c2f7f105b15d44ed374916d0cd35'), 'verify_md5 worked in a negative test';




