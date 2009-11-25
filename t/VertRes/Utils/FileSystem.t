#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd 'cwd';

BEGIN {
    use Test::Most tests => 22;
    
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
$fsu->rmtree($test_dir);
ok ! -d $test_dir, 'rmtree removed a directory that contained a file';
# copy... impossible to test the diff part of copy failure?
$file = $fsu->catfile('t', 'data', 'io_test.txt');
my $copy = $fsu->catfile($tmp_dir, 'copy_test');
ok $fsu->copy($file, $copy), 'simple copy test';
ok -s $copy, 'copy test really did work';
my $devnull = $fsu->catfile('dev', 'null', 'copy_test');
my $ok = $fsu->copy($file, $devnull);
is $ok, 0, 'copy fails to /dev/null';

# verify_md5 calculate_md5 ...

undef $fsu;
ok ! -d $tmp_dir, 'tmpdir destroyed ok';
