#!/usr/bin/perl -w
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 11;
    
    use_ok('VertRes::Utils::Hierarchy');
    use_ok('VertRes::IO');
    use_ok('File::Copy');
}

my $h_util = VertRes::Utils::Hierarchy->new();
isa_ok $h_util, 'VertRes::Base';

# setup our files and paths
my $io = VertRes::IO->new();
my $si_file = $io->catfile('t', 'data', 'sequence.index');
ok -s $si_file, 'sequence.index file ready to test with';
my $temp_dir = $io->tempdir();

my $lane_path = '/path/to/LowCov-CEU/NA06986/SLX/Solexa-5459/SRR003670';
is_deeply {$h_util->parse_lane($lane_path)}, {study => 'LowCov-CEU',
                                              sample => 'NA06986',
                                              platform => 'SLX',
                                              library => 'Solexa-5459',
                                              lane => 'SRR003670'}, 'parse_lane test';

ok $h_util->check_lanes_vs_sequence_index([$lane_path], $si_file), 'check_lanes_vs_sequence_index test';

# setup for and test create_release_hierarchy()
my $release_dir = $io->catfile($temp_dir, 'REL');
my $mapping_tar_gz = $io->catfile($temp_dir, 'mapping.tar.gz');
copy($io->catfile('t', 'data', 'mapping_hierarchy_with_bams.tar.gz'), $mapping_tar_gz);
ok -s $mapping_tar_gz, 'mapping_hierarchy_with_bams.tar.gz copied and ready to test with';
system("tar -xz -C $temp_dir -f $mapping_tar_gz");
ok -d $io->catfile($temp_dir, 'mapping'), 'mapping hierarchy extracted and ready to test with';
my @mapping_paths;
my @expected_rel_paths;
foreach my $lane_path ('mapping/proj1/ind1/SLX/lib1/lane1',
                       'mapping/proj1/ind1/SLX/lib1/lane2',
                       'mapping/proj2/ind2/454/lib3/lane3/',
                       'mapping/proj2/ind2/454/lib3/lane4/') {
    push(@mapping_paths, $io->catfile($temp_dir, $lane_path));
    my $rel_path = $lane_path;
    $rel_path =~ s/^mapping\///;
    push(@expected_rel_paths, $io->catfile($temp_dir, 'REL', $rel_path));
}
my @expected_rel_bams = ('pe_raw.sorted.bam', 'se_raw.sorted.bam', 'se_raw.sorted.bam', 'pe_raw.sorted.bam');
my @expected_rel_bam_paths;
foreach my $rel_path (@expected_rel_paths) {
    my $expected_bam = shift @expected_rel_bams;
    push(@expected_rel_bam_paths, $io->catfile($rel_path, $expected_bam));
}

is_deeply [$h_util->create_release_hierarchy(\@mapping_paths, $release_dir)], \@expected_rel_bam_paths, 'create_release_hierarchy test';

my $created_links = 0;
foreach my $link (@expected_rel_bam_paths) {
    $created_links++ if -e $link;
}
is $created_links, 4, 'create_release_hierarchy created the correct bam symlinks';

exit;
