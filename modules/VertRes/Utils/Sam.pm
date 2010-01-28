=head1 NAME

VertRes::Utils::Sam - sam utility functions

=head1 SYNOPSIS

use VertRes::Utils::Sam;

my $sam_util = VertRes::Utils::Sam->new();

# use any of the utility functions described here, eg.
my $are_similar = $sam_util->bams_are_similar(@bam_files);

=head1 DESCRIPTION

General utility functions for working on or with sam/bam files.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Utils::Sam;

use strict;
use warnings;
use File::Copy;
use File::Basename;
use File::Spec;
use VertRes::IO;
use VertRes::Utils::FileSystem;
use VertRes::Wrapper::samtools;
use VertRes::Wrapper::picard;
use HierarchyUtilities;
use VertRes::Parser::dict;
use VertRes::Parser::sequence_index;
use VertRes::Parser::sam;
use VertRes::Utils::FastQ;
use VertRes::Utils::Math;
use VertRes::Utils::Hierarchy;
use Digest::MD5;

use base qw(VertRes::Base);

use VertRes::Parser::sam;
our %flags = (paired_tech    => '0x0001',
              paired_map     => '0x0002',
              self_unmapped  => '0x0004',
              mate_unmapped  => '0x0008',
              self_reverse   => '0x0010',
              mate_reverse   => '0x0020',
              '1st_in_pair'  => '0x0040',
              '2nd_in_pair'  => '0x0080',
              not_primary    => '0x0100',
              failed_qc      => '0x0200',
              duplicate      => '0x0400');

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Utils::Sam->new();
 Function: Create a new VertRes::Utils::Sam object.
 Returns : VertRes::Utils::Sam object
 Args    : java_memory => int (for methods that call picard, this will be passed
                               on to the picard wrapper, which has a default)

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    return $self;
}

=head2 bams_are_similar

 Title   : bams_are_similar
 Usage   : my $are_similar = $obj->bams_are_similar(@bam_files);
 Function: Check if multiple name-sorted bam files are the same, ignoring
           mapping quality 0 reads.
 Returns : boolean (false if they differ)
 Args    : list of bam files

=cut

sub bams_are_similar {
    my ($self, @files) = @_;
    
    @files >= 2 or $self->throw("At least 2 bam files are needed");
    
    my $sw = VertRes::Wrapper::samtools->new(verbose => $self->verbose,
                                             run_method => 'open',
                                             quiet => 1);
    
    # open
    my @in_fhs;
    foreach my $file (@files) {
        my $bamfh = $sw->view($file, undef, h => 1);
        push(@in_fhs, $bamfh);
    }
    
    my $first_fh = $in_fhs[0];
    while (! eof($first_fh)) {
        # if a line in one of the input bams has mapping quality 0, we skip it
        # along with the corresponding line in all other input bams, even if
        # they have better mapping qualities
        my @lines;
        while (1) {
            @lines = ();
            my $good_mapping_qs = 0;
            foreach my $fh (@in_fhs) {
                my $line = <$fh>;
                push(@lines, $line);
                
                if ($line =~ /^\@/) {
                    $good_mapping_qs++;
                    next;
                }
                
                my (undef, undef, undef, undef, $mapping_q) = split("\t", $line);
                $good_mapping_qs++ if $mapping_q > 0;
            }
            
            last if $good_mapping_qs == @in_fhs;
            last if eof($first_fh);
        }
        
        my %lines;
        foreach my $line (@lines) {
            if ($line =~ /^\@/) {
                $lines{$line}++;
                next;
            }
            
            my ($col1, $col2, $col3, $col4, $mapping_q) = split("\t", $line);
            $lines{"$col1\t[...]\t$col3\t$col4\t[...]"}++;
        }
        
        unless (keys %lines == 1) {
            $self->debug("Files started differing here:\n".join("\n", keys %lines));
            return 0;
        }
        
        my ($line, $count) = each %lines;
        return 0 unless $count == @in_fhs;
    }
    
    foreach my $fh (@in_fhs) {
        return 0 unless eof($fh);
    }
    
    return 1;
}

=head2 split_bam_by_sequence

 Title   : split_bam_by_sequence
 Usage   : my @bams = $obj->split_bam_by_sequence($bam_file);
 Function: Splits a bam file into multiple bam files, one for each sequence the
           reads were aligned to.
 Returns : list of the bam files it made
 Args    : bam filename, optionally these hash options:
           ignore => 'regex' to ignore sequences with ids that match the regex,
                     eg. ignore => '^N[TC]_\d+' to avoid making splits for
                     contigs
           output_dir => 'path' to specify where the split bams are created;
                         default is the same dir as the input bam

=cut

sub split_bam_by_sequence {
    my ($self, $bam, %opts) = @_;
    
    my $sw = VertRes::Wrapper::samtools->new(verbose => $self->verbose,
                                             run_method => 'open',
                                             quiet => 1);
    
    my $basename = basename($bam);
    my $output_dir;
    if (defined $opts{output_dir}) {
        $output_dir = $opts{output_dir};
    }
    else {
        $output_dir = $bam;
        $output_dir =~ s/$basename$//;
    }
    
    # find out what sequences there are
    my $bamfh = $sw->view($bam, undef, H => 1);
    my $sp = VertRes::Parser::sam->new(fh => $bamfh);
    my %all_sequences = $sp->sequence_info();
    
    $sw->run_method('system');
    
    # we must have a bam index file
    my $bai = $bam.'.bai';
    my $created_bai = 0;;
    unless (-s $bai) {
        $sw->index($bam, $bai);
        $sw->run_status >= 1 || $self->throw("Failed to create $bai");
        $created_bai = 1;
    }
    
    # make a split for each sequence
    my @out_bams;
    foreach my $seq (keys %all_sequences) {
        next if ($opts{ignore} && $seq =~ /$opts{ignore}/);
        
        my $out_bam = File::Spec->catfile($output_dir, $seq.'.'.$basename);
        
        $sw->view($bam, $out_bam, regions => [$seq], h => 1, b => 1);
        $sw->run_status >= 1 || $self->throw("Failed to create split $out_bam");
        
        push(@out_bams, $out_bam);
    }
    
    if ($created_bai) {
        unlink($bai);
    }
    
    return @out_bams;
}

=head2 add_sam_header

 Title   : add_sam_header
 Usage   : $obj->add_sam_header('headerless.sam',
                                sequence_index => 'sequence.index',
                                lane_path => '/path/to/lane');
 Function: Adds more meta infomation to the header of a lane-level sam file.
           Also, if the records don't have RG tags, these will be added (or
           corrected).
 Returns : boolean (true on success - the sam filename doesn't change)
 Args    : path to sam and:
           either just sequence_index => 'sequence.index' (to calculate the
           below from it)
           -or-
           all of these key => value pairs:
           sample_name => string, eg. NA00000;
           library => string
           platform => string
           centre => string
           insert_size => int
           project => string, the study id, eg. SRP000001
 
           -and-

           either just lane_path => '/path/to/lane' (to calculate the below
           using HierarchyUtilities::lane_info)
           -or-
           all of these key => value pairs:
           lane => run_id, eg. SRR00000
           ref_fa => '/full/path/to/reference.fasta'
           ref_dict => '/full/path/to/reference.dict'
           ref_name => standard name for the reference, eg. NCBI36

           -optionally-
           program => name of the program used to do the mapping
           program_version => version of the program
           run_name => string, the platform unit - the original accession before
                       DCC gave it the [ES]RR id

=cut

sub add_sam_header {
    my ($self, $raw_sam_file, %args) = @_;
    -s $raw_sam_file || $self->throw("sam file '$raw_sam_file' doesn't exist, can't add header to it!");
    
    # gather information
    my $lane_path = $args{lane_path};
    my %lane_info = %{HierarchyUtilities::lane_info($lane_path)} if $lane_path;
    my $lane = $args{lane} || $lane_info{lane} || $self->throw("lane must be supplied");
    my $ref_fa = $args{ref_fa} || $lane_info{fa_ref} || $self->throw("the reference fasta must be supplied");
    my $ref_dict = $args{ref_dict} || $lane_info{dict_ref} || $self->throw("the reference dict must be supplied");
    my $ref_name = $args{ref_name} || $lane_info{ref_name} || $self->throw("the reference assembly name must be supplied");
    
    my $seq_index = $args{sequence_index} || '';
    my $parser = VertRes::Parser::sequence_index->new(file => $seq_index) if $seq_index;
    my $individual = $args{sample_name} || $parser->lane_info($lane, 'sample_name') || $lane_info{sample} || $self->throw("Unable to get sample_name for $lane from sequence index '$seq_index' or other supplied args");
    my $library = $args{library} || $parser->lane_info($lane, 'LIBRARY_NAME') || $lane_info{library};
    my $insert_size = $args{insert_size};
    unless (defined $insert_size) {
        $insert_size = $parser->lane_info($lane, 'INSERT_SIZE') || $lane_info{insert_size};
    }
    my $run_name = $args{run_name};
    if (! defined $run_name && $parser) {
        $run_name = $parser->lane_info($lane, 'run_name');
    }
    my $platform = $args{platform} || $parser->lane_info($lane, 'INSTRUMENT_PLATFORM');
    my $centre = $args{centre} || $parser->lane_info($lane, 'CENTER_NAME');
    my $project = $args{project} || $parser->lane_info($lane, 'STUDY_ID');
    if ($parser && (! $project || $project eq '-')) {
        my $sn = $parser->lane_info($lane, 'STUDY_NAME');
        if ($sn && $sn ne '-') {
            $project = $sn;
        }
    }
    
    # write the sam header
    my $headed_sam = $raw_sam_file.'.withheader';
    my $header_lines = 0;
    open(my $shfh, '>', $headed_sam) or $self->throw("Cannot create $headed_sam: $!");
    print $shfh "\@HD\tVN:1.0\tSO:coordinate\n";
    $header_lines++;
    
    # md5 in the sam spec is supposed to be of the sequence in the uppercase
    # with no spaces, which can be found with my sequence_dicter.pl script,
    # which was used to generate the .dict files for our references
    my $dict_parser = VertRes::Parser::dict->new(file => $ref_dict);
    my $rh = $dict_parser->result_holder;
    while ($dict_parser->next_result) {
        print $shfh "\@SQ\tSN:$rh->[0]\tLN:$rh->[1]\tAS:$ref_name\tM5:$rh->[3]\tUR:file:$ref_fa\n";
        $header_lines++;
    }
    
    print $shfh "\@RG\tID:$lane\tLB:$library\tSM:$individual";
    print $shfh "\tPU:$run_name" if (defined $run_name);
    print $shfh "\tPI:$insert_size" if (defined $insert_size);
    print $shfh "\tCN:$centre" if (defined $centre);
    print $shfh "\tPL:$platform" if (defined $platform);
    print $shfh "\tDS:$project" if (defined $project);
    
    print $shfh "\n";
    $header_lines++;
    
    if ($args{program}) {
        print $shfh "\@PG\tID:$args{program}";
        print $shfh "\tVN:$args{program_version}" if (defined $args{program_version});
        print $shfh "\n";
        $header_lines++;
    }
    
    # combine the header with the raw sam file, adding/correcting RG tag if
    # necessary, ignoring existing header
    open(my $rsfh, '<', $raw_sam_file) or $self->throw("Couldn't open raw sam file '$raw_sam_file'");
    my $expected_lines = 0;
    while (<$rsfh>) {
        next if /^@/;
        $expected_lines++;
        unless (/\tRG:Z:/) {
            s/\n$/\tRG:Z:$lane\n/;
        }
        else {
            s/\tRG:Z:\S+/\tRG:Z:$lane/;
        }
        print $shfh $_;
    }
    close($rsfh);
    close($shfh);
    
    # check and mv
    my $io = VertRes::IO->new(file => $raw_sam_file);
    my $raw_lines = $io->num_lines();
    $expected_lines += $header_lines;
    $io->file($headed_sam);
    my $actual_lines = $io->num_lines();
    
    if ($expected_lines == $actual_lines) {
        move($headed_sam, $raw_sam_file) || $self->throw("Failed to move $headed_sam to $raw_sam_file: $!");
        return 1;
    }
    else {
        unlink($headed_sam);
        $self->throw("Failed to prepend sam header to sam file '$raw_sam_file' (expected $expected_lines lines but got $actual_lines");
    }
}

=head2 sam_to_fixed_sorted_bam

 Title   : sam_to_fixed_sorted_bam
 Usage   : $obj->sam_to_fixed_sorted_bam('headed.sam', 'sorted.bam', 'ref.fa');
 Function: Converts a sam file to a mate-fixed and sorted bam file. Also runs
           the bam via fillmd to get accurate MD and NM tags. (Both maq and bwa
           get these wrong; ssaha doesn't give it at all.)
           NB: the input sam must already have headers, eg. by using
           add_sam_header() on it. 
 Returns : boolean (true on success)
 Args    : starting sam file, output name for bam file, reference fasta used to
           make the sam. Optionally, args to pass to VertRes::Wrapper::samtools
           (eg. quiet => 1).

=cut

sub sam_to_fixed_sorted_bam {
    my ($self, $in_sam, $out_bam, $ref_fa, @args) = @_;
    
    my $io = VertRes::IO->new();
    my $fsu = VertRes::Utils::FileSystem->new();
    my $temp_dir = $fsu->tempdir();
    my $tmp_bam = $fsu->catfile($temp_dir, 'fixed_sorted.bam');
    
    # sam -> fixed, sorted bam
    my $wrapper = VertRes::Wrapper::samtools->new(verbose => $self->verbose, @args);
    $wrapper->sam_to_fixed_sorted_bam($in_sam, $tmp_bam);
    unless ($wrapper->run_status() >= 1) {
        $self->warn("Failed the initial sam->bam step");
        return 0;
    }
    
    # bam -> fillmd'd bam
    my $tmp2_bam = $fsu->catfile($temp_dir, 'fillmd.bam');
    $wrapper->quiet(1);
    $wrapper->fillmd($tmp_bam, $ref_fa, $tmp2_bam, b => 1);
    
    # check it
    $wrapper->run_method('open');
    my $fh = $wrapper->view($tmp2_bam, undef, h => 1);
    my $bam_count = 0;
    while (<$fh>) {
        $bam_count++;
    }
    close($fh);
    $io->file($in_sam);
    my $sam_count = $io->num_lines();
    if ($bam_count >= $sam_count) {
        move($tmp2_bam, $out_bam) || $self->throw("Failed to move $tmp2_bam to $out_bam: $!");;
        return 1;
    }
    else {
        $self->warn("$tmp2_bam is bad (only $bam_count lines vs $sam_count)");
        return 0;
    }
}

=head2 rmdup

 Title   : rmdup
 Usage   : $obj->rmdup('in.bam', 'rmduped.bam');
 Function: Removes duplicates from a bam file, automatically taking into account
           if these are single ended or paired reads.
 Returns : boolean (true on success)
 Args    : starting bam file, output name for bam file, and either
           lane_path => '/path/to/lane' (to decide the below using
           HierarchyUtilities::getFastqInfo)
           -or-
           single_ended => boolean (false by default, ie. paired)

           Optionally, args to pass to VertRes::Wrapper::samtools (eg.
           quiet => 1)

=cut

sub rmdup {
    my ($self, $in_bam, $out_bam, %args) = @_;
    
    my $command = 'rmdup';
    
    my @args = ();
    my $lane_path = delete $args{lane_path};
    if ($lane_path) {
        my @fastq_info = @{HierarchyUtilities::getFastqInfo($lane_path)};
        if ($fastq_info[0]->[0] && ! $fastq_info[1]->[0] && ! $fastq_info[2]->[0]) {
            @args = (s => 1);
        }
    }
    if (delete $args{single_ended}) {
        @args = (s => 1);
    }
    
    my $wrapper = VertRes::Wrapper::samtools->new(verbose => $self->verbose, %args);
    $wrapper->rmdup($in_bam, $out_bam, @args);
    
    return $wrapper->run_status() >= 1;
}

=head2 merge

 Title   : merge
 Usage   : $obj->merge('out.bam', @bams_for_merging);
 Function: Merges bam files using picard-tools. If only one bam supplied for
           merging, just symlinks it to out.bam.
 Returns : boolean (true on success)
 Args    : output bam filename, list of input bam filenames

=cut

sub merge {
    my ($self, $out_bam, @in_bams) = @_;
    return unless @in_bams;
    
    unlink($out_bam);
    
    if (@in_bams == 1) {
        # we want a relative symlink, so get the path of out_bam and make the
        # symlink relative to that
        my ($basename, $base) = fileparse($out_bam);
        my $rel_path = File::Spec->abs2rel($in_bams[0], $base);
        return symlink($rel_path, $out_bam);
    }
    
    my $fsu = VertRes::Utils::FileSystem->new();
    my $verbose = $self->verbose();
    my $wrapper = VertRes::Wrapper::picard->new(verbose => $verbose,
                                                quiet => $verbose ? 0 : 1,
                                                $self->{java_memory} ? (java_memory => $self->{java_memory}) : (),
                                                validation_stringency => 'silent',
                                                tmp_dir => $fsu->tempdir());
    
    $wrapper->merge_and_check($out_bam, \@in_bams);
    
    return $wrapper->run_status() >= 1;
}

=head2 stats

 Title   : stats
 Usage   : $obj->stats('in.bam', 'in2.bam', ...);
 Function: Generate both flagstat and bas files for the given bam files.
 Returns : boolean (true on success)
 Args    : list of bam files (for each, a .bas and .flagstat file will be
           created)

=cut

sub stats {
    my ($self, @in_bams) = @_;
    
    my $wrapper = VertRes::Wrapper::samtools->new(verbose => $self->verbose);
    my $all_ok = 1;
    foreach my $in_bam (@in_bams) {
        my $flagstat = $in_bam.'.flagstat';
        $wrapper->flagstat($in_bam, $flagstat.'.tmp');
        unless ($wrapper->run_status() >= 1) {
            $all_ok = 0;
            unlink($flagstat);
        }
        else {
            move($flagstat.'.tmp', $flagstat) || ($all_ok = 0);
        }
        
        my $bas = $in_bam.'.bas';
        my $ok = $self->bas($in_bam, $bas.'.tmp');
        unless ($ok) {
            $all_ok = 0;
            unlink($bas);
        }
        else {
            move($bas.'.tmp', $bas) || ($all_ok = 0);
        }
    }
    
    return $all_ok;
}

=head2 bas

 Title   : bas
 Usage   : $obj->bas('in.bam', 'out.bas');
 Function: Generate a 'bas' file. These are bam statistic files that provide
           various stats about each readgroup in a bam.
 Returns : boolean (true on success)
 Args    : input bam file (must have been run via samtools fillmd so that there
           are accurate NM tags for each record), output filename, optionally a
           sequence.index file from the DCC if working on bams with poor headers

=cut

sub bas {
    my ($self, $in_bam, $out_bas, $seq_index) = @_;
    
    open(my $bas_fh, '>', $out_bas) || $self->throw("Couldn't write to '$out_bas': $!");
    
    my $expected_lines = 0;
    
    # print header
    print $bas_fh join("\t", 'bam_filename', 'md5', 'study', 'sample', 'platform',
                             'library', 'readgroup', '#_total_bases',
                             '#_mapped_bases', '#_total_reads',
                             '#_mapped_reads',
                             '#_mapped_reads_paired_in_sequencing',
                             '#_mapped_reads_properly_paired',
                             '%_of_mismatched_bases',
                             'average_quality_of_mapped_bases',
                             'mean_insert_size', 'insert_size_sd',
                             'median_insert_size',
                             'insert_size_median_absolute_deviation'), "\n";
    $expected_lines++;
    
    # get the stats for each read group
    my %readgroup_data = $self->bam_statistics($in_bam);
    
    # get the meta data
    my $hu = VertRes::Utils::Hierarchy->new(verbose => $self->verbose);
    my $dcc_filename = $hu->dcc_filename($in_bam);
    
    my $md5;
    my $md5_file = $in_bam.'.md5';
    if (-s $md5_file) {
        open(my $md5fh, $md5_file) || $self->throw("Couldn't open md5 file '$md5_file': $!");
        my $line = <$md5fh>;
        ($md5) = split(' ', $line);
    }
    else {
        my $dmd5 = Digest::MD5->new();
        open(FILE, $in_bam) or $self->throw("Couldn't open bam '$in_bam': $!");
        binmode(FILE);
        $dmd5->addfile(*FILE);
        $md5 = $dmd5->hexdigest;
    }
    
    my $sip = VertRes::Parser::sequence_index->new(file => $seq_index) if $seq_index;
    
    my $stw = VertRes::Wrapper::samtools->new(quiet => 1);
    $stw->run_method('open');
    my $view_fh = $stw->view($in_bam, undef, H => 1);
    $view_fh || $self->throw("Failed to samtools view '$in_bam'");
    my $ps = VertRes::Parser::sam->new(fh => $view_fh);
    
    # add in the meta data
    while (my ($rg, $data) = each %readgroup_data) {
        $readgroup_data{$rg}->{dcc_filename} = $dcc_filename;
        $readgroup_data{$rg}->{study} = $ps->readgroup_info($rg, 'DS') || 'unknown_study';
        $readgroup_data{$rg}->{sample} = $ps->readgroup_info($rg, 'SM') || 'unknown_sample';
        $readgroup_data{$rg}->{platform} = $ps->readgroup_info($rg, 'PL') || 'unknown_platform';
        $readgroup_data{$rg}->{library} = $ps->readgroup_info($rg, 'LB') || 'unknown_library';
        
        # fall back on the sequence.index if we have to
        if ($sip) {
            my %convert = (study => 'STUDY_ID',
                           sample => 'SAMPLE_NAME',
                           platform => 'INSTRUMENT_PLATFORM',
                           library => 'LIBRARY_NAME');
            foreach my $field (keys %convert) {
                if (! $readgroup_data{$rg}->{$field} || $readgroup_data{$rg}->{$field} eq "unknown_$field" || $readgroup_data{$rg}->{$field} eq '-') {
                    $readgroup_data{$rg}->{$field} = $sip->lane_info($rg, $convert{$field}) || "unknown_$field";
                }
            }
        }
        
        $readgroup_data{$rg}->{md5} = $md5;
    }
    
    # print stats per-readgroup
    foreach my $rg (sort keys %readgroup_data) {
        my %data = %{$readgroup_data{$rg}};
        print $bas_fh join("\t", $data{dcc_filename},
                                 $data{md5},
                                 $data{study},
                                 $data{sample},
                                 $data{platform},
                                 $data{library},
                                 $rg,
                                 $data{total_bases} || 0,
                                 $data{mapped_bases} || 0,
                                 $data{total_reads} || 0,
                                 $data{mapped_reads} || 0,
                                 $data{mapped_reads_paired_in_seq} || 0,
                                 $data{mapped_reads_properly_paired} || 0,
                                 $data{percent_mismatch},
                                 $data{avg_qual},
                                 $data{avg_isize},
                                 $data{sd_isize},
                                 $data{median_isize},
                                 $data{mad}), "\n";
        $expected_lines++;
    }
    close($bas_fh);
    
    my $io = VertRes::IO->new(file => $out_bas);
    my $actual_lines = $io->num_lines;
    if ($actual_lines == $expected_lines) {
        return 1;
    }
    else {
        unlink($out_bas);
        $self->warn("Wrote $expected_lines to $out_bas, but only read back $actual_lines! Deleted the output.");
        return 0;
    }
}

=head2 bam_statistics

 Title   : bam_statistics
 Usage   : my %readgroup_stats = $obj->bam_statistics('in.bam');
 Function: Calculate all the stats per-readgroup for a bam needed for the bas
           format.
 Returns : hash with keys as readgroup ids and values as hash refs. Those refs
           have the keys:
           total_bases, mapped_bases, total_reads, mapped_reads,
           mapped_reads_paired_in_seq, mapped_reads_properly_paired,
           percent_mismatch, avg_qual, avg_isize, sd_isize, median_isize, mad
 Args    : input bam file (must have been run via samtools fillmd so that there
           are accurate NM tags for each record)

=cut

sub bam_statistics {
    my ($self, $bam) = @_;
    
    # go through the bam and accumulate all the raw stats in little memory
    my $ps = VertRes::Parser::sam->new(file => $bam, verbose => -1);
    my $fqu = VertRes::Utils::FastQ->new();
    
    my %readgroup_data;
    my $previous_rg = 'unknown_readgroup';
    while (my ($qname, $read_length, $mapped_length, $flag, $qual_str, $mapq, $isize, $rg, $nm) =
           $ps->get_fields('QNAME', 'SEQ_LENGTH', 'MAPPED_SEQ_LENGTH', 'FLAG', 'QUAL', 'MAPQ', 'ISIZE', 'RG', 'NM')) {
        unless ($rg) {
            $self->warn("$qname had no RG tag, using previous RG tag '$previous_rg'");
            $rg = $previous_rg;
        }
        $previous_rg = $rg;
        
        my @this_rg_data = @{$readgroup_data{$rg} || []};
        $this_rg_data[0]++;
        $this_rg_data[1] += $read_length;
        
        if ($ps->is_mapped($flag)) {
            $this_rg_data[2]++;
            $this_rg_data[3] += $mapped_length;
            $this_rg_data[4]++ if $ps->is_sequencing_paired($flag);
            
            # avg quality of mapped bases
            foreach my $qual ($fqu->qual_to_ints($qual_str)) {
                $this_rg_data[5]++;
                
                if ($this_rg_data[5] == 1) {
                    $this_rg_data[6] = $qual;
                }
                else {
                    $this_rg_data[6] += ($qual - $this_rg_data[6]) / $this_rg_data[5];
                }
            }
            
            # avg insert size and keep track of 's' for later calculation of sd.
            # algorithm based on http://www.johndcook.com/standard_deviation.html
            if ($ps->is_mapped_paired($flag)) {
                $this_rg_data[7]++;
                
                $isize ||= 0;
                
                if ($mapq > 0) {
                    if ($isize > 0) { # avoids counting the isize twice for a pair, since one will be negative
                        $this_rg_data[8]++;
                        
                        if ($this_rg_data[8] == 1) {
                            $this_rg_data[9] = $isize;
                            $this_rg_data[10] = 0;
                        }
                        else {
                            my $old_mean = $this_rg_data[9];
                            $this_rg_data[9] += ($isize - $old_mean) / $this_rg_data[8];
                            $this_rg_data[10] += ($isize - $old_mean) * ($isize - $this_rg_data[9]);
                        }
                        
                        # also, median insert size. Couldn't find an accurate
                        # running algorithm, but just keeping a histogram is
                        # accurate and uses less than 1MB. We can use the same
                        # histogram later to calculate the MAD.
                        $this_rg_data[11]->{$isize}++;
                    }
                }
            }
            
            # for later calculation of mismatch %
            if (defined $nm && $nm ne '*') {
                $this_rg_data[12] += $read_length;
                $this_rg_data[13] += $nm;
            }
        }
        
        $readgroup_data{$rg} = \@this_rg_data;
    }
    
    # calculate the means etc.
    my %stats;
    my $math_util = VertRes::Utils::Math->new();
    foreach my $rg (sort keys %readgroup_data) {
        my @data = @{$readgroup_data{$rg}};
        
        # calculate/round stats
        my $avg_qual = defined $data[5] ? sprintf("%0.2f", $data[6]) : 0;
        my $avg_isize = defined $data[8] ? sprintf("%0.0f", $data[9]) : 0;
        my $sd_isize = $avg_isize ? sprintf("%0.2f", sqrt($data[10] / $data[8])) : 0;
        my $percent_mismatch = defined $data[13] ? sprintf("%0.2f", (100 / $data[12]) * $data[13]) : 0;
        
        my $median_isize = 0;
        my $mad = 0;
        if (defined $data[11]) {
            $median_isize = $math_util->histogram_median($data[11]);
            
            my %ads;
            while (my ($isize, $freq) = each %{$data[11]}) {
                my $ad = abs($median_isize - $isize);
                $ads{$ad} += $freq;
            }
            
            $mad = $math_util->histogram_median(\%ads);
        }
        
        my %rg_stats;
        $rg_stats{total_reads} = $data[0];
        $rg_stats{total_bases} = $data[1];
        $rg_stats{mapped_reads} = $data[2];
        $rg_stats{mapped_bases} = $data[3];
        $rg_stats{mapped_reads_paired_in_seq} = $data[4];
        $rg_stats{mapped_reads_properly_paired} = $data[7];
        $rg_stats{avg_qual} = $avg_qual;
        $rg_stats{avg_isize} = $avg_isize;
        $rg_stats{sd_isize} = $sd_isize;
        $rg_stats{percent_mismatch} = $percent_mismatch;
        $rg_stats{median_isize} = $median_isize;
        $rg_stats{mad} = $mad;
        $stats{$rg} = \%rg_stats;
    }
    
    return %stats;
}

=head2 make_unmapped_bam

 Title   : make_unmapped_bam
 Usage   : $obj->make_unmapped_bam('in.bam', 'out.bam');
 Function: Given a bam file, generates another bam file that contains only the
           unmapped reads.
 Returns : boolean (true on success)
 Args    : starting bam file, output name for bam file.

=cut

sub make_unmapped_bam {
    my ($self, $in_bam, $out_bam) = @_;
    
    # create a new sam file with just the header and unmapped reads
    my $in = VertRes::Wrapper::samtools->new(file => $in_bam,
                                             run_method => 'open',
                                             verbose => $self->verbose);
    my $in_fh = $in->view($in_bam, undef, h => 1);
    my $filtered_sam = $in_bam.'.unmapped_sam';
    $self->register_for_unlinking($filtered_sam);
    my $out = VertRes::IO->new(file => ">$filtered_sam");
    my $out_fh = $out->fh();
    my $sp = VertRes::Parser::sam->new();
    
    while (<$in_fh>) {
        if (/^\@/) {
            print $out_fh $_;
        }
        else {
            my (undef, $flag) = split;
            unless ($sp->is_mapped($flag)) {
                print $out_fh $_;
            }
        }
    }
    close($in_fh);
    close($out_fh);
    
    # convert to a bam file
    my $samtools = VertRes::Wrapper::samtools->new(verbose => $self->verbose);
    $samtools->view($filtered_sam, $out_bam, b => 1, S => 1);
    
    return $samtools->run_status() >= 1;
}

=head2 add_unmapped

 Title   : add_unmapped
 Usage   : $obj->add_unmapped($sam_file, @fastqs);
 Function: Append reads in the fastqs not already present in the sam file to the
           end of the sam file.
 Returns : boolean
 Args    : sam file name, fastq file name(s)

=cut

sub add_unmapped {
    my ($self, $sam, @fqs) = @_;
    
    my $paired = @fqs > 1;
    
    # store all the read names we already have in the sam. We can't use a hash
    # since that uses too much memory, and a tied hash is too slow. Just output
    # to disk and sort afterwards with unix sort
    my $snames_file = $sam.'.snames';
    my $snames_sorted_file = $snames_file.'.sorted';
    my $snames_sorted_fixed_file = $snames_sorted_file.'.fixed';
    my $num_mapped = 0;
    unless (-s $snames_sorted_fixed_file) {
        unless (-s $snames_sorted_file) {
            unless (-s $snames_file) {
                open(my $rfh, '>', $snames_file) || $self->throw("Could not write to $snames_file");
                open(my $sfh, $sam) || $self->throw("Could not open $sam");
                my $sp = VertRes::Parser::sam->new();
                while (<$sfh>) {
                    next if /^@/;
                    my ($qname, $flag) = split;
                    
                    unless ($qname =~ /\/([12])$/) {
                        my $read_num;
                        if ($sp->is_first($flag)) {
                            $read_num = 1;
                        }
                        elsif ($sp->is_second($flag)) {
                            $read_num = 2;
                        }
                        else {
                            # don't know what read this is; we'll call it read 1 and if both
                            # of these turn out to be present we'll add read 2 as well.
                            $read_num = 1;
                        }
                        
                        $qname =~ s/\/$//;
                        $qname .= "/$read_num";
                    }
                    
                    $num_mapped++;
                    print $rfh $qname, "\n";
                }
                close($sfh);
                close($rfh);
                
                # check rnames file isn't truncated
                my $actual_count = VertRes::IO->new(file => $snames_file)->num_lines;
                unless ($actual_count == $num_mapped) {
                    unlink($snames_file);
                    $self->throw("made an rnames file but it was truncated!");
                }
            }
            
            $num_mapped ||= VertRes::IO->new(file => $snames_file)->num_lines;
            
            system("sort $snames_file > $snames_sorted_file");
            my $actual_count = VertRes::IO->new(file => $snames_sorted_file)->num_lines;
            unless ($actual_count == $num_mapped) {
                unlink($snames_sorted_file);
                $self->throw("made an rnames_sorted_file file but it was truncated!");
            }
            
            unlink($snames_file);
        }
        
        # if we have the same read name twice in a row, change the second one
        # to be the second read of a pair
        open(my $ofh, '>', $snames_sorted_fixed_file) || $self->throw("Could not write to $snames_sorted_fixed_file");
        open(my $ifh, $snames_sorted_file) || $self->throw("Could not open $snames_sorted_file");
        my $previous_name = '';
        $num_mapped = 0;
        while (<$ifh>) {
            $num_mapped++;
            
            if ($_ eq $previous_name) {
                s/1\n$/2\n/;
            }
            
            print $ofh $_;
            
            $previous_name = $_;
        }
        close($ofh);
        close($ifh);
        
        # check rnames_sorted_fixed_file file isn't truncated
        my $actual_count = VertRes::IO->new(file => $snames_sorted_fixed_file)->num_lines;
        unless ($actual_count == $num_mapped) {
            unlink($snames_sorted_fixed_file);
            $self->throw("made an rnames_sorted_fixed_file file but it was truncated!");
        }
        
        unlink($snames_sorted_file);
    }
    $num_mapped ||= VertRes::IO->new(file => $snames_sorted_fixed_file)->num_lines;
    
    # now list out to another file all the read names in the fastqs and sort
    # that
    my $qnames_file = $sam.'.qnames';
    my $qnames_sorted_file = $snames_file.'.sorted';
    my $qlines = 0;
    unless (-s $qnames_sorted_file) {
        unless (-s $qnames_file) {
            open(my $qfh, '>', $qnames_file) || $self->throw("Could not write to $qnames_file");
            foreach my $fq (@fqs) {
                my $fqp = VertRes::Parser::fastq->new(file => $fq);
                my $rh = $fqp->result_holder();
                
                # loop through all the sequences in the fastq (without indexing
                # to save memory)
                while ($fqp->next_result(1)) {
                    $qlines++;
                    print $qfh $rh->[0], "\n";
                }
            }
            close($qfh);
            
            # check it's not truncated
            my $actual_count = VertRes::IO->new(file => $qnames_file)->num_lines;
            unless ($actual_count == $qlines) {
                unlink($qnames_file);
                $self->throw("made a qnames_file file but it was truncated!");
            }
        }
        
        $qlines ||= VertRes::IO->new(file => $qnames_file)->num_lines;
        
        system("sort $qnames_file > $qnames_sorted_file");
        my $actual_count = VertRes::IO->new(file => $qnames_sorted_file)->num_lines;
        unless ($actual_count == $qlines) {
            unlink($qnames_sorted_file);
            $self->throw("made a qnames_sorted_file file but it was truncated!");
        }
        
        unlink($qnames_file);
    }
    
    # do a simple diff to find the reads we're missing
    open(my $snamesfh, $snames_sorted_fixed_file) || $self->throw("Could not open $snames_sorted_fixed_file");
    open(my $qnamesfh, $qnames_sorted_file) || $self->throw("Could not open $qnames_sorted_file");
    my %missing_reads;
    while (my $sam_read = <$snamesfh>) {
        while (my $fastq_read = <$qnamesfh>) {
            if ($fastq_read ne $sam_read) {
                chomp($fastq_read);
                $missing_reads{$fastq_read} = 1;
            }
            else {
                last;
            }
        }
    }
    close($snamesfh);
    close($qnamesfh);
    
    # now for all the reads we're missing, append to original sam
    if (keys %missing_reads) {
        open(my $sfh, '>>', $sam) || $self->throw("Could not append to $sam");
        foreach my $fq (@fqs) {
            my $fqp = VertRes::Parser::fastq->new(file => $fq);
            my $rh = $fqp->result_holder();
            
            while ($fqp->next_result(1)) {
                my $id = $rh->[0];
                next unless exists $missing_reads{$id};
                my $seq = $rh->[1];
                my $qual = $rh->[2];
                
                # work out the SAM flags
                my ($mate_unmapped, $first, $second) = (0, 0, 0);
                if ($paired) {
                    my ($read_num) = $id =~ /\/([12])$/;
                    $read_num || $self->throw("Couldn't work out read number for $id");
                    if ($read_num == 1) {
                        $first = 1;
                    }
                    else {
                        $second = 1;
                    }
                    
                    my $mate_num = $read_num == 1 ? 2 : 1;
                    my $mate_id = $id;
                    $mate_id =~ s/$read_num$/$mate_num/;
                    if (exists $missing_reads{$mate_id}) {
                        $mate_unmapped = 1;
                    }
                }
                
                my $flags = $self->calculate_flag(self_unmapped => 1,
                                                  paired_tech => $paired,
                                                  mate_unmapped => $mate_unmapped,
                                                  '1st_in_pair' => $first,
                                                  '2nd_in_pair' => $second);
                
                print $sfh "$id\t$flags\t*\t0\t0\t*\t*\t0\t0\t$seq\t$qual\n";
            }
        }
        close($sfh);
    }
    
    unlink($snames_sorted_fixed_file);
    unlink($qnames_sorted_file);
    
    my $lines_now = VertRes::IO->new(file => $sam)->num_lines;
    my $unmapped_lines = keys %missing_reads;
    my $lines_expected = $num_mapped + $unmapped_lines;
    unless ($lines_now >= $lines_expected) {
        $self->warn("Tried appending $unmapped_lines unmapped records to the $num_mapped mapped records sam file, but ended up with only $lines_now lines!");
        return 0;
    }
    
    return 1;
}

=head2 calculate_flag

 Title   : calculate_flag
 Usage   : my $flag = $obj->calculate_flag(mapped => 1);
 Function: Make a sam flag that has the desired meaning. 
 Returns : int
 Args    : hash of desired meaning, where keys are one or more of the following
           and values are boolean:
           paired_tech
           paired_map
           self_unmapped
           mate_unmapped
           self_reverse
           mate_reverse
           '1st_in_pair'
           '2nd_in_pair'
           not_primary 
           failed_qc
           duplicate

=cut

sub calculate_flag {
    my ($self, %desired) = @_;
    
    # first, fix to make sane
    if ($desired{'self_unmapped'}) {
        foreach my $meaning ('duplicate', 'self_reverse', 'paired_map') {
            if ($desired{$meaning}) {
                $self->warn("'self_unmapped' was set, but so was '$meaning'; forcing $meaning off");
                $desired{$meaning} = 0;
            }
        }
    }
    if ($desired{'1st_in_pair'} && $desired{'2nd_in_pair'}) {
        $self->warn("Can't have both 1st_in_pair and 2nd_in_pair; forcing both off");
        $desired{'1st_in_pair'} = 0;
        $desired{'2nd_in_pair'} = 0;
    }
    if (! $desired{'paired_tech'}) {
        foreach my $meaning ('paired_map', 'mate_unmapped', 'mate_reverse', '1st_in_pair', '2nd_in_pair') {
            if ($desired{$meaning}) {
                $self->warn("'$meaning' was set, but 'paired_tech' was not; forcing paired_tech on");
                $desired{'paired_tech'} = 1;
            }
        }
    }
    
    my $result = 0;
    while (my ($flag, $val) = each %desired) {
        $val || next;
        exists $flags{$flag} || next;
        $result += hex($flags{$flag});
    }
    
    return $result;
}

=head2 check_bam_header

 Title   : check_bam_header
 Usage   : if ($obj->check_bam_header('a.bam',
                                    readgroup1 => { sample_name => 'NA000001' })
 Function: Prior to calling rewrite_bam_header() with the same arguments, just
           check that a rewrite is even necessary.
 Returns : boolean
 Args    : as per check_bam_header()

=cut

sub check_bam_header {
    my ($self, $bam, %rg_changes) = @_;
    
    keys %rg_changes || return 1;
    my %arg_to_tag = (sample_name => 'SM',
                      library => 'LB',
                      platform => 'PL',
                      centre => 'CN',
                      insert_size => 'PI',
                      project => 'DS');
    
    my $stin = VertRes::Wrapper::samtools->new(quiet => 1, run_method => 'open');
    my $bamfh = $stin->view($bam, undef, H => 1);
    my $made_changes = 0;
    while (<$bamfh>) {
        if (/^\@RG.+ID:([^\t]+)/) {
            my $rg = $1;
            if (exists $rg_changes{$rg}) {
                while (my ($arg, $value) = each %{$rg_changes{$rg}}) {
                    next unless defined $value;
                    my $tag = $arg_to_tag{$arg} || next;
                    
                    if (/\t$tag:([^\t\n]+)/) {
                        if ($1 ne $value) {
                            $made_changes = 1;
                            # all of <$bamfh> must be gone through, otherwise
                            # things will break, so we don't 'last' here.
                        }
                    }
                    else {
                        $made_changes = 1;
                        # no 'last'; see above comment
                    }
                }
            }
        }
    }
    close($bamfh);
    
    return $made_changes;
}

=head2 rewrite_bam_header

 Title   : rewrite_bam_header
 Usage   : $obj->rewrite_bam_header('a.bam',
                                    readgroup1 => { sample_name => 'NA000001' });
 Function: Corrects the @RG header lines in a bam header, eg. to account for
           any swaps (say, the sample changed since the bam was made).
 Returns : boolean (true on success, meaning a change was made and the bam has
           been confirmed OK, or no change was necessary and the bam was
           untouched)
 Args    : path to sam and a hash with keys as readgroup identifiers and values
           as hash refs with at least one of the following:
           sample_name => string, eg. NA00000;
           library => string
           platform => string
           centre => string
           insert_size => int
           project => string, the study id, eg. SRP000001 (overwrites DS)

=cut

sub rewrite_bam_header {
    my ($self, $bam, %rg_changes) = @_;
    
    keys %rg_changes || return 1;
    my %arg_to_tag = (sample_name => 'SM',
                      library => 'LB',
                      platform => 'PL',
                      centre => 'CN',
                      insert_size => 'PI',
                      project => 'DS');
    
    my $stin = VertRes::Wrapper::samtools->new(quiet => 1, run_method => 'open');
    my $stout = VertRes::Wrapper::samtools->new(quiet => 1, run_method => 'write_to');
    
    # output to bam
    my $temp_bam = $bam.'.rewrite_header.tmp.bam';
    $self->register_for_unlinking($temp_bam);
    my $sfh = $stout->view(undef, $temp_bam, b => 1, S => 1);
    $sfh || $self->throw("failed to get a filehandle for writing to '$temp_bam'");
    
    # parse just the header through perl to make changes
    my $expected_lines = 0;
    my $made_changes = 0;
    my $bamfh = $stin->view($bam, undef, H => 1);
    while (<$bamfh>) {
        $expected_lines++;
        
        if (/^\@RG.+ID:([^\t]+)/) {
            my $rg = $1;
            if (exists $rg_changes{$rg}) {
                while (my ($arg, $value) = each %{$rg_changes{$rg}}) {
                    next unless defined $value;
                    my $tag = $arg_to_tag{$arg} || next;
                    
                    if (/\t$tag:([^\t\n]+)/) {
                        if ($1 ne $value) {
                            $made_changes = 1;
                            s/\t$tag:[^\t\n]+/\t$tag:$value/;
                        }
                    }
                    else {
                        $made_changes = 1;
                        s/\n/\t$tag:$value\n/
                    }
                }
            }
        }
        
        print $sfh $_;
    }
    close($bamfh);
    
    unless ($made_changes) {
        close($sfh);
        unlink($temp_bam);
        return 1;
    }
    
    # then quickly (no regex) append all the bam records
    $bamfh = $stin->view($bam);
    while (<$bamfh>) {
        $expected_lines++;
        print $sfh $_;
    }
    close($bamfh);
    close($sfh);
    
    # check the new bam for truncation
    $bamfh = $stin->view($temp_bam, undef, h => 1);
    my $new_bam_lines = 0;
    while (<$bamfh>) {
        $new_bam_lines++;
    }
    close($bamfh);
    
    if ($new_bam_lines == $expected_lines) {
        unlink($bam);
        move($temp_bam, $bam) || $self->throw("Failed to move $temp_bam to $bam: $!");;
        return 1;
    }
    else {
        $self->warn("$temp_bam is bad (only $new_bam_lines lines vs $expected_lines), will unlink it");
        unlink($temp_bam);
        return 0;
    }
}

=head2 rewrite_bas_meta

 Title   : rewrite_bas_meta
 Usage   : $obj->rewrite_bas_meta('a.bam.bas',
                                  readgroup1 => { sample_name => 'NA000001' });
 Function: Corrects the meta columns 1-6 in a bas file, eg. to account for
           any swaps (say, the sample changed since the bas was made).
 Returns : boolean (true on success, meaning a change was made and the bas has
           been confirmed OK, or no change was necessary and the bas was
           untouched)
 Args    : path to bas and a hash with keys as readgroup identifiers and values
           as hash refs with at least one of the following:
           sample_name => string, eg. NA00000;
           library => string
           platform => string
           project => string, the study id, eg. SRP000001
           md5 => string (NB: if you've just used rewrite_bam_header on a bam
                          and are making the same changes to its bas, the bam
                          md5 will have changed so you should supply the new
                          md5 here)

=cut

sub rewrite_bas_meta {
    my ($self, $bas, %rg_changes) = @_;
    
    keys %rg_changes || return 1;
    my %arg_to_col = (sample_name => 3,
                      library => 5,
                      project => 2,
                      md5 => 1);
    # NAXXXXX.[chromN].technology.[center].algorithm.study_id.YYYY_MM.bam
    my %arg_to_filename_regex = (sample_name => qr{^[^\.]+(\.)},
                                 platform => qr{(?:ABI_SOLID|ILLUMINA|LS454|SLX|454|SOLID)(\.)},
                                 project => qr{[^\.]+(\.\d{4}_)});
    
    my $temp_bas = $bas.'.rewrite_bas_meta.tmp.bas';
    $self->register_for_unlinking($temp_bas);
    open(my $ofh, '>', $temp_bas) || $self->throw("Could not write to '$temp_bas'");
    
    open(my $ifh, $bas) || $self->throw("Could not open '$bas'");
    my $made_changes = 0;
    my $expected_lines = 0;
    while (<$ifh>) {
        $expected_lines++;
        
        my @cols = split("\t", $_);
        my $rg = $cols[6];
        
        if (exists $rg_changes{$rg}) {
            while (my ($arg, $value) = each %{$rg_changes{$rg}}) {
                my $col = $arg_to_col{$arg};
                if ($col) {
                    if ($cols[$col] ne $value) {
                        $made_changes = 1;
                        $cols[$col] = $value;
                    }
                }
                
                my $regex = $arg_to_filename_regex{$arg};
                if ($regex) {
                    my $orig = $cols[0];
                    $cols[0] =~ s/$regex/$value$1/;
                    if ($cols[0] ne $orig){
                        $made_changes = 1;
                    }
                }
            }
        }
        
        print $ofh join("\t", @cols);
    }
    close($ifh);
    close($ofh);
    
    unless ($made_changes) {
        unlink($temp_bas);
        return 1;
    }
    
    # check the new bas for truncation
    my $new_bas_lines = VertRes::IO->new(file => $temp_bas)->num_lines;
    
    if ($new_bas_lines == $expected_lines) {
        unlink($bas);
        move($temp_bas, $bas) || $self->throw("Failed to move $temp_bas to $bas: $!");;
        return 1;
    }
    else {
        $self->warn("$temp_bas is bad (only $new_bas_lines lines vs $expected_lines), will unlink it");
        unlink($temp_bas);
        return 0;
    }
}

1;
