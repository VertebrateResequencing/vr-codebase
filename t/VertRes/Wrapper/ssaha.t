#!/usr/bin/perl -w
use strict;
use warnings;
use File::Copy;

BEGIN {
    use Test::Most tests => 16;
    
    use_ok('VertRes::Wrapper::ssaha');
    use_ok('VertRes::IO');
}

my $debug = 0;
my $quiet = $debug ? 0 : 1;
my $verbose = $debug ? 1 : 0;
my $ssaha = VertRes::Wrapper::ssaha->new(quiet => $quiet, verbose => $verbose);
isa_ok $ssaha, 'VertRes::Wrapper::WrapperI';
is $ssaha->quiet, $quiet, 'quiet set via new';

is $ssaha->exe, 'ssaha2', 'exe ok';
like $ssaha->version, qr/\d\.\d.\d/, 'version ok';

# prepare our test files; copy them to a temp dir where we will do the tests
my $io = VertRes::IO->new();
my $temp_dir = $io->tempdir;
my $read1_basename = '2822_6_1_1000.fastq';
my $read2_basename = '2822_6_2_1000.fastq';
my $ref_basename = 'S_suis_P17.dna';
my $read1 = $io->catfile($temp_dir, $read1_basename);
my $read2 = $io->catfile($temp_dir, $read2_basename);
my $ref = $io->catfile($temp_dir, $ref_basename);
copy($io->catfile('t', 'data', $read1_basename), $read1);
copy($io->catfile('t', 'data', $read2_basename), $read2);
copy($io->catfile('t', 'data', $ref_basename), $ref);
ok -s $read1, 'test file 1 ready to use';
ok -s $read2, 'test file 2 ready to use';
ok -s $ref, 'test file 3 ready to use';

# the files we expect to be created
my $mapping = $io->catfile($temp_dir, 'mapping.sam');
my @ref_files;
foreach my $suffix ('head', 'body', 'name', 'base', 'size') {
    push(@ref_files, $io->catfile($temp_dir, "$ref_basename.$suffix"));
}

# build the reference index hash files
$ssaha->ssaha2Build($ref, $ref_basename, skip => 3);
is $ssaha->run_status, 1, 'status after Build is ok';
my $ok_ref_files = 0;
my @mtimes;
foreach my $ref_file (@ref_files) {
    -s $ref_file || next;
    $ok_ref_files++;
    my @stat = stat($ref_file);
    push(@mtimes, $stat[9]);
}
is $ok_ref_files, 5, 'made the correct number of ref index hash files';

# run the whole mapping
$ssaha->do_mapping(ref => $ref,
                   read1 => $read1,
                   read2 => $read2,
                   output => $mapping,
                   insert_size => 2000,
                   local_cache => $temp_dir);
is $ssaha->run_status, 1, 'status after mapping is ok';
ok -s $mapping, 'output file exists';
my ($lines, $mapped) = check_sam($mapping);
is $lines, 2000, 'output sam not truncated';
cmp_ok $mapped, '>=', 1596, 'mapped enough reads';

# doing the mapping didn't re-create the index hash files
my @new_mtimes;
foreach my $ref_file (@ref_files) {
    my @stat = stat($ref_file);
    push(@new_mtimes, $stat[9]);
}
is_deeply \@new_mtimes, \@mtimes, 'do_mapping with local_cache equal to ref dir didn\'t repeat the Build';


exit;

sub check_sam {
    my $sam = shift;
    open(my $sfh, $sam) || return (0, 0);
    my $lines = 0;
    my $mapped = 0;
    while (<$sfh>) {
        next unless $_ =~ /\S/;
        $lines++;
        my @a = split;
        $mapped++ if $a[3];
    }
    return ($lines, $mapped);
}
