#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use File::Spec;
use File::HomeDir;
use File::Slurp qw(read_dir);
use Data::Dumper;

BEGIN {
    use Test::Most;
    
    use_ok('VertRes::Wrapper::tophat');
    use_ok('VertRes::Utils::Mappers::tophat');
    use VertRes::Wrapper::fastqcheck;
    use VertRes::Utils::FileSystem;
}


my $tophat = VertRes::Wrapper::tophat->new(quiet => 1);
isa_ok $tophat, 'VertRes::Wrapper::WrapperI';
is $tophat->quiet, 1, 'quiet set via new';

is $tophat->exe, 'tophat', 'exe ok';
like $tophat->version, qr/\d+\.\d+/, 'version ok';


# prepare our test files; copy them to a temp dir where we will do the tests
my $fsu = VertRes::Utils::FileSystem->new();

#If running this test, you're home folder will be used for the temp folder where everything will happen
my $temp_dir = $fsu->tempdir( dir => File::HomeDir->my_home );

my $read1_basename = '2822_6_1_1000.fastq';
my $read2_basename = '2822_6_2_1000.fastq';
my $ref_basename = 'S_suis_P17.fa';
my $read1 = File::Spec->catfile($temp_dir, $read1_basename);
my $read2 = File::Spec->catfile($temp_dir, $read2_basename);
my $ref = File::Spec->catfile($temp_dir, $ref_basename);

copy(File::Spec->catfile('t', 'data', $read1_basename), $read1);
copy(File::Spec->catfile('t', 'data', $read2_basename), $read2);
copy(File::Spec->catfile('t', 'data', $ref_basename), $ref);

ok -s $read1, 'test file 1 ready to use';
ok -s $read2, 'test file 2 ready to use';
ok -s $ref, 'test file 3 ready to use';

# the files we expect to be created
my $sai1 = File::Spec->catfile($temp_dir, '2822_6_1_1000.sai');
my $sai2 = File::Spec->catfile($temp_dir, '2822_6_2_1000.sai');
my $mapping = File::Spec->catfile($temp_dir, 'mapping.sam');

my @ref_index_files;
foreach my $suffix (qw(1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2)) {
    push(@ref_index_files, File::Spec->catfile($temp_dir, 'S_suis_P17.fa.'.$suffix));
}

my $fw = VertRes::Wrapper::fastqcheck->new(quiet => 1);
my $out_file1 = File::Spec->catfile($temp_dir, '2822_6_1_1000.fastq.fastqcheck');
my $out_file2 = File::Spec->catfile($temp_dir, '2822_6_2_1000.fastq.fastqcheck');
$fw->run($read1, $out_file1);
$fw->run($read2, $out_file2);


# run the whole mapping
$tophat->do_mapping(ref => $ref,
                 read1 => $read1,
                 read2 => $read2,
                 output => $mapping,
                 insert_size => 500);
is $tophat->run_status, 2, 'status after mapping is ok';
print "MAPPING: $mapping\n";

my $mapping_dir = get_mapping_dir( $temp_dir );

my %merging_op_file_list = (
			    accepted_sam => {
					     filepath => File::Spec->catfile( $mapping_dir, 'accepted_hits.sam' )
					    },
			    accepted_bam => {
					     filepath => File::Spec->catfile( $mapping_dir, 'accepted_hits.bam' )
					    },
			    unmapped_bam => {
					     filepath => File::Spec->catfile( $mapping_dir, 'unmapped.bam' )
					    },
			    qfixed_unmapped_sam => {
						    filepath => File::Spec->catfile( $mapping_dir, 'qfixed_unmapped.sam' )
						   },
			    qfixed_unmapped_bam => {
						    filepath => File::Spec->catfile( $mapping_dir, 'qfixed_unmapped.bam' )
						   },
			    merged_bam => {
					   filepath => File::Spec->catfile( $mapping_dir, 'merged.bam' )
					  },
			    merge_sorted_bam => {
						 filepath => File::Spec->catfile( $mapping_dir, 'merged.sorted.bam' )
						},
			    merge_sorted_fixed_bam => {
						       filepath => File::Spec->catfile( $mapping_dir, 'merged.sorted.fixed.bam' )
						      }
			   );

for my $filetype ( keys %merging_op_file_list ) {

  my $filepath = $merging_op_file_list{$filetype}{filepath};

  ok -e $filepath, "$filepath exists";
  ok -s $filepath, "$filepath has data";


  $merging_op_file_list{$filetype}{read_count} = `samtools view -c $filepath` if $filetype =~ m/_bam$/;
  chomp($merging_op_file_list{$filetype}{read_count}) if $filetype =~ m/_bam$/;

}

is $merging_op_file_list{accepted_bam}{read_count}, 492, "matches number of expected reads in $merging_op_file_list{accepted_bam}{filepath}";

is $merging_op_file_list{unmapped_bam}{read_count}, 1508, "matches number of expected reads in $merging_op_file_list{unmapped_bam}{filepath}";

is $merging_op_file_list{qfixed_unmapped_bam}{read_count}, 1508, "matches number of expected reads in $merging_op_file_list{qfixed_unmapped_bam}{filepath}";

is $merging_op_file_list{merged_bam}{read_count}, 2000, "matches number of expected reads in $merging_op_file_list{merged_bam}{filepath}";

is $merging_op_file_list{merge_sorted_bam}{read_count}, 2000, "matches number of expected reads in $merging_op_file_list{merge_sorted_bam}{filepath}";

is $merging_op_file_list{merge_sorted_fixed_bam}{read_count}, 2000, "matches number of expected reads in $merging_op_file_list{merge_sorted_fixed_bam}{filepath}";

is $merging_op_file_list{accepted_bam}{read_count} + $merging_op_file_list{unmapped_bam}{read_count}, $merging_op_file_list{merge_sorted_fixed_bam}{read_count}, 
  "Number of reads in $merging_op_file_list{accepted_bam}{filepath} PLUS reads in $merging_op_file_list{unmapped_bam}{filepath}\n\t match the number of reads in $merging_op_file_list{merge_sorted_fixed_bam}{filepath}";


my ($lines, $mapped) = check_sam($mapping);
is $lines, 2003, 'output sam not truncated';
cmp_ok $mapped, '>=', 964, 'mapped enough reads';
foreach my $file ($mapping, @ref_index_files) {
    ok -s $file, "output file exists $file";
    unlink($file);
}

# individual commands
$tophat->setup_reference($ref);
is $tophat->run_status, 2, 'status after index is ok';
foreach my $file (@ref_index_files) {
    ok -s $file, "output file exists $file";
}

done_testing;
exit;

sub get_mapping_dir {

  my ($temp_dir) = @_;

  my $mapping_dir;
  for my $dir ( grep { -d "$temp_dir/$_" } read_dir( $temp_dir ) ) {
    $mapping_dir = "$temp_dir/$dir";
  }

  return $mapping_dir;

}

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

