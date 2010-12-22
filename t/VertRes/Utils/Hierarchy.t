#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;
use File::Copy;

BEGIN {
    use Test::Most tests => 25;
    
    use_ok('VertRes::Utils::Hierarchy');
    use_ok('VRTrack::Lane');
}

my $h_util = VertRes::Utils::Hierarchy->new();
isa_ok $h_util, 'VertRes::Base';

# setup our files and paths
my $si_file = File::Spec->catfile('t', 'data', 'sequence.index.2010');
ok -s $si_file, 'sequence.index file ready to test with';
my $bam_file = File::Spec->catfile('t', 'data', 'headed1.bam');
ok -s $bam_file, 'bam file ready to test with';
my $fsu = VertRes::Utils::FileSystem->new();
my $temp_dir = $fsu->tempdir();

# dcc_filename test
is $h_util->dcc_filename($bam_file, '20100208', $si_file), 'NA11918.ILLUMINA.bwa.CEU.low_coverage.20100208', 'dcc_filename test';
my $chrom_bam = File::Spec->catfile('t', 'data', 'NA06985.chrom7.ILLUMINA.bwa.SRP000031.20091216.bam');
is $h_util->dcc_filename($chrom_bam, '20100208', $si_file), 'NA06985.chrom7.ILLUMINA.bwa.CEU.low_coverage.20100208', 'dcc_filename chrom test';
is $h_util->dcc_filename($chrom_bam, '20100208'), 'NA06985.chrom7.ILLUMINA.bwa.unknown_population.unknown_analysisgroup.20100208', 'dcc_filename chrom test without sequence.index';

# lane details
my $lane_path = '/path/to/META/CEU_low_coverage/NA06986/SLX/Solexa-5459/SRR003670';
is_deeply {$h_util->parse_lane($lane_path)}, {study => 'CEU_low_coverage',
                                              sample => 'NA06986',
                                              platform => 'SLX',
                                              library => 'Solexa-5459',
                                              lane => 'SRR003670'}, 'parse_lane test';

ok $h_util->check_lanes_vs_sequence_index([$lane_path], $si_file), 'check_lanes_vs_sequence_index test';

SKIP: {
    eval { require VRTrack::Testconfig; };
    skip "Skipping some tests because VRTrack db has not been configured", 15 if $@;
    
    my $connection_details = { database => VRTrack::Testconfig->config('test_db'),
                               host     => VRTrack::Testconfig->config('host'),
                               port     => VRTrack::Testconfig->config('port'),
                               user     => VRTrack::Testconfig->config('user'),
                               password => VRTrack::Testconfig->config('password') };
    
    # set up the test VRTrack db to be like the one from 1000 genomes, for a
    # single lane
    ok my $vrtrack = VRTrack::VRTrack->new($connection_details), 'setup vrtrack accessing the test db';
    clear_vrtrack_db($vrtrack);
    my $update_vrmeta_script = File::Spec->catfile(qw(scripts update_vrmeta.pl));
    my $sample_file = File::Spec->catfile(qw(t data vrtrack_vertres_1kg.samples));
    my $seq_ind_file = File::Spec->catfile(qw(t data vrtrack_vertres_1kg.si));
    system("$update_vrmeta_script --database $connection_details->{database} --samples $sample_file --index $seq_ind_file > /dev/null 2> /dev/null");
    $vrtrack = VRTrack::VRTrack->new($connection_details);
    
    $lane_path = '/path/to/META/SRP000031/NA06986/SLX/Solexa_5459/SRR003670';
    ok $h_util->check_lanes_vs_database([$lane_path], $vrtrack), 'check_lanes_vs_database test';
    #ok ! $h_util->check_lanes_vs_database([$lane_path], $vrtrack, 1), 'check_lanes_vs_database check all test'; # too slow...
    
    # setup for and test create_release_hierarchy(), with new db contents
    clear_vrtrack_db($vrtrack);
    $sample_file = File::Spec->catfile(qw(t data vrtrack_vertres_test.samples));
    $seq_ind_file = File::Spec->catfile(qw(t data vrtrack_vertres_test.si));
    system("$update_vrmeta_script --database $connection_details->{database} --samples $sample_file --index $seq_ind_file > /dev/null 2> /dev/null");
    $vrtrack = VRTrack::VRTrack->new($connection_details);
    my $all_lanes = $vrtrack->processed_lane_hnames();
    my $count = 0;
    foreach my $lane_hname (sort @{$all_lanes}) {
        my $lane = VRTrack::Lane->new_by_hierarchy_name($vrtrack, $lane_hname);
        $lane->is_processed("import", 1);
        $lane->is_processed("mapped", 1);
        foreach my $file (@{$lane->files()}) {
            $file->is_processed("import", 1); $file->update;
        }
        my $mapping = $lane->add_mapping();
        my $assembly = $mapping->assembly("NCBI37");
        if ( ! $assembly) {
            $assembly = $mapping->add_assembly("NCBI37");
        }
        my $exe;
        $count++;
        if ($count <= 2) {
            $exe = "bwa";
        }
        else {
            $exe = "ssaha2";
        }
        my $mapper = $mapping->mapper($exe, "1.5");
        if ( ! $mapper) {
            $mapper = $mapping->add_mapper($exe, "1.5");
        }
        $mapping->update;
        $lane->update;
    }
    
    my $release_dir = File::Spec->catfile($temp_dir, 'REL');
    my $mapping_tar_gz = File::Spec->catfile($temp_dir, 'mapping.tar.gz');
    copy(File::Spec->catfile('t', 'data', 'mapping_hierarchy_with_bams.tar.gz'), $mapping_tar_gz);
    ok -s $mapping_tar_gz, 'mapping_hierarchy_with_bams.tar.gz copied and ready to test with';
    system("tar -xz -C $temp_dir -f $mapping_tar_gz");
    ok -d File::Spec->catfile($temp_dir, 'mapping'), 'mapping hierarchy extracted and ready to test with';
    my @mapping_paths;
    my @expected_rel_paths;
    foreach my $lane_path ('mapping/proj1/ind1/SLX/lib1/lane1', # 1.pe.raw.sorted.bam
                           'mapping/proj1/ind1/SLX/lib1/lane2', # 4.se.raw.sorted.bam
                           'mapping/proj2/ind2/454/lib3/lane3/', # 7.se.raw.sorted.bam
                           'mapping/proj2/ind2/454/lib3/lane4/') { # 10.pe.raw.sorted.bam
        push(@mapping_paths, File::Spec->catfile($temp_dir, $lane_path));
        my $rel_path = $lane_path;
        $rel_path =~ s/^mapping\///;
        push(@expected_rel_paths, File::Spec->catfile($temp_dir, 'REL', $rel_path));
    }
    my @expected_rel_bams = ('pe.bam', 'se.bam', 'se.bam', 'pe.bam');
    my @expected_rel_bam_paths;
    foreach my $rel_path (@expected_rel_paths) {
        my $expected_bam = shift @expected_rel_bams;
        push(@expected_rel_bam_paths, File::Spec->catfile($rel_path, $expected_bam));
    }
    
    is_deeply [$h_util->create_release_hierarchy(\@mapping_paths, $release_dir,
                                                 vrtrack => $vrtrack,
                                                 slx_mapper => 'bwa',
                                                 '454_mapper' => 'ssaha',
                                                 assembly_name => 'NCBI37')], \@expected_rel_bam_paths, 'create_release_hierarchy test';
    
    my $created_links = 0;
    foreach my $link (@expected_rel_bam_paths) {
        $created_links++ if -e $link;
    }
    is $created_links, 4, 'create_release_hierarchy created the correct bam symlinks';
    clear_vrtrack_db($vrtrack);
    
    # the next few tests use the complete tracking db from the mouse genomes
    # project, which is too big to distribute in t/data, so these tests run
    # internally only
    SKIP: {
        my $node = $ENV{HOST};
        skip "Skipping some tests because they only work internal to the Sanger Institute", 9 unless $node =~ /^farm2-head\d+/;
        
        $connection_details = { host => $ENV{VRTRACK_HOST}, port => $ENV{VRTRACK_PORT}, user => $ENV{VRTRACK_RO_USER}, database => 'mouse_reseq_track' };
        
        # lane_info test
        ok my %info = $h_util->lane_info('3034_8', db => $connection_details, qc_passed => 1, mapped => 1), 'lane_info ran ok';
        is_deeply [$info{technology}, $info{seq_tech}, $info{project}, $info{study}, $info{sample}], ['ILLUMINA', 'SLX', '129P2 Mouse Genome', 'ERP000034', '129P2_1'], 'lane_info gave correct technology and seq_tech, project, study and sample';
        cmp_ok $info{individual_coverage}, '>=', 22.39, 'lane_info had the correct individual_coverage';
        # needs better tests for all the different ways of calling lane_info...
        
        # hierarchy_coverage test
        my $cov = $h_util->hierarchy_coverage(individual => ['129P2/OlaHsd'], db => $connection_details, qc_passed => 1, mapped => 1);
        cmp_ok $cov, '>=', 22.39, 'hierarchy_coverage direct test';
        my $cov2 = $h_util->hierarchy_coverage(lane => '3034_8', level => 'individual', db => $connection_details, qc_passed => 1, mapped => 1);
        cmp_ok $cov, '>=', 22.39, 'hierarchy_coverage lane level test';
        is $cov, $cov2, 'coverage was the same in both tests';
        
        # netapp-related methods
        is_deeply [$h_util->nfs_disks], [qw(/nfs/vertreseq01 /nfs/vertreseq02 /nfs/vertreseq03 /nfs/vertreseq04 /nfs/vertreseq05 /nfs/vertreseq06 /nfs/vertreseq07 /nfs/vertreseq08 /nfs/vertreseq09 /nfs/vertreseq10 /nfs/vertreseq11 /nfs/vertreseq12 /nfs/vertreseq13 /nfs/vertreseq14 /nfs/vertreseq15 /nfs/vertreseq16)], 'nfs_disks returned the expected disks';
        like $h_util->nfs_disk, qr{/nfs/vertreseq\d\d}, 'nfs_disk returns one of the disks'; # *** hard to test if it's the correct one though...
        my $lane = VRTrack::Lane->new_by_hierarchy_name($vrtrack, 'lane1');
        $ENV{DATA_HIERARCHY} = 'project:sample:technology:library:lane';
        like $h_util->lane_storage_path($lane), qr{/nfs/vertreseq\d\d/hashed_lanes/vrtrack_vertres_test/6/b/8/2/lane1}, 'lane_storage_path gave the correct path';
    }
}

# *** store_lane can't really be tested? implemented with other well-tested
#     things though

# *** fix_simple_swaps test??


exit;

sub clear_vrtrack_db {
    my $vrtrack = shift;
    
    my $dbh = $vrtrack->{_dbh};
    foreach ($dbh->tables()){
        next if /`latest_/;
        next if /schema_version/;
        $dbh->do("TRUNCATE TABLE $_");
    }
}
