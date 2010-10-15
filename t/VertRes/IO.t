#!/usr/bin/env perl
use strict;
use warnings;
use Cwd 'cwd';

BEGIN {
    use Test::Most tests => 41;
    
    use_ok('VertRes::IO');
    use_ok('VertRes::Utils::FileSystem');
}

my $io = VertRes::IO->new();
isa_ok $io, 'VertRes::Base';
isa_ok $io, 'VertRes::IO';
my $fsu = VertRes::Utils::FileSystem->new();

# open and read from a normal file
my $file = $fsu->catfile('t', 'data', 'io_test.txt');
ok -f $io->file($file), 'file returns a functioning path to the file';
ok my $fh = $io->fh(), 'fh returned something';
is ref($fh), 'GLOB', 'fh returned a glob';
my @expected_lines = qw(foo bar);
while (<$fh>) {
    chomp;
    next unless $_;
    my $exp = shift @expected_lines;
    is $_, $exp, 'filehandle must have been ok';
}
$io->close();

# and from a gzip compressed file
$file = $fsu->catfile('t', 'data', 'io_test.txt.gz');
ok -f $io->file($file), 'file returns a functioning path to the file';
ok $fh = $io->fh(), 'fh returned something for a gzip .gz';
is ref($fh), 'IO::Uncompress::Gunzip', 'fh returned a IO::Uncompress::Gunzip';
@expected_lines = qw(foo bar);
while (<$fh>) {
    chomp;
    next unless $_;
    my $exp = shift @expected_lines;
    is $_, $exp, 'filehandle must have been ok for a gzip .gz';
}
$io->close();

# and from a bgzip compressed file
$file = $fsu->catfile('t', 'data', 'io_test.bgzip.gz');
ok -f $io->file($file), 'file returns a functioning path to the file';
ok $fh = $io->fh(), 'fh returned something for a bgzip .gz';
is ref($fh), 'GLOB', 'fh returned a GLOB';
my $lines = 0;
while (<$fh>) {
    chomp;
    next unless $_;
    $lines++;
}
is $lines, 2000, 'all 2000 lines of the bgzip file were readable';
$io->close();

# get a temp file
($fh, $file) = $fsu->tempfile;
close($fh);
# write to it
is $io->file(">$file"), $file, 'filename not broken when opening for write';
$fh = $io->fh();
ok print($fh "foo bar\nbar foo\nboo far"), 'could print to an output file';
$io->close();
$io->file($file);
$fh = $io->fh();
is <$fh>, "foo bar\n", 'could read back what we wrote';
is $io->num_lines, 3, 'number of lines correct, even after manually using the filehandle';
is <$fh>, "bar foo\n", 'could read the next line even after getting the number of all lines';
$io->close();
ok -e $file, 'file exists after a close';

$io = VertRes::IO->new();
my $tmp_dir = $fsu->tempdir;
# write compressed
my $gz_file = $fsu->catfile($tmp_dir, 'test.gz');
ok $io->file(">$gz_file"), 'could set up write to .gz file';
$fh = $io->fh();
ok print($fh "compressed test\n"), 'could print to a gz file';
$io->close();
ok $io->file($gz_file), 'the .gz file was created';
$fh = $io->fh();
is <$fh>, "compressed test\n", 'could read back what we wrote';

# parse a fod file
$io = VertRes::IO->new();
$file = $fsu->catfile('t', 'data', 'fod.txt');
my $cwd = cwd();
my @expected;
foreach my $dir ('data', 'VertRes') {
    push(@expected, $fsu->catfile($cwd, 't', $dir));
}
is_deeply [$io->parse_fod($file)], [sort @expected], 'parse_fod test';

# parse a fofn file  ***ideally should have a symlink inside the fofn as well...
$file = $fsu->catfile('t', 'data', 'fofn.txt');
@expected = ();
foreach my $file ('2822_6_1_1000.fastq', 'S_suis_P17.fa.fai', 'fastq.gz.fastqcheck') {
    push(@expected, $fsu->catfile($cwd, 't', 'data', $file));
}
is_deeply [$io->parse_fofn($file)], [sort @expected], 'parse_fofn test';

# download a remote file
$io = VertRes::IO->new();
my $remote_md5 = '4b62815f42eeadf93eb7835a99d9fbb5';
my @remote_content = ("A test file for VertRes::IO->get_remote_file()\n", "it\n", "has\n", "4 lines.\n");
# via http
my $http_url = 'http://wwwdev.sanger.ac.uk/modelorgs/mousegenomes/VertResIO.get_remote_file.test_file.txt';
ok my $http_file = $io->get_remote_file($http_url), 'get_remote_file on http returned a downloaded file path';
test_file($http_file);
unlink($http_file);
# via ftp, with md5 check and save to custom location
my $ftp_url = 'ftp://ftp.sanger.ac.uk/pub/1000genomes/sb10/VertResIO.get_remote_file.test_file.txt';
$tmp_dir = $fsu->tempdir;
my $custom_save_location = $fsu->catfile($tmp_dir, 'ftp_download.txt');
ok my $ftp_file = $io->get_remote_file($ftp_url, md5 => $remote_md5, save => $custom_save_location), 'get_remote_file on ftp returned a downloaded file path';
is $ftp_file, $custom_save_location, 'get_remote_file with save option save file to correct place';
test_file($ftp_file);
unlink($ftp_file);
# with the file() shorthand
ok $ftp_file = $io->file($ftp_url), 'file() on a url returned something';
test_file($ftp_file);
$fh = $io->fh;
is <$fh>, $remote_content[0], 'fh() readline on a url file worked transparently';

exit;

sub test_file {
    my $file = shift;
    ok open(my $fh, $file), 'downloaded file could be opened';
    my @download_content = <$fh>;
    is_deeply \@download_content, \@remote_content, 'downloaded file content was correct';
    close($fh);
}
