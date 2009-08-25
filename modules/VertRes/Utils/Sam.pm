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
 Args    : n/a

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
           run_name => string, the platform unit, eg. 7563
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
    my $insert_size = $args{insert_size} || $parser->lane_info($lane, 'INSERT_SIZE') || $lane_info{insert_size};
    my $run_name = $args{run_name} || $parser->lane_info($lane, 'run_name') || $self->throw("Unable to get run_name for $lane");
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
    
    print $shfh "\@RG\tID:$lane\tPU:$run_name\tLB:$library\tSM:$individual";
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
    my $temp_dir = $io->tempdir();
    my $tmp_bam = $io->catfile($temp_dir, 'fixed_sorted.bam');
    
    # sam -> fixed, sorted bam
    my $wrapper = VertRes::Wrapper::samtools->new(verbose => $self->verbose, @args);
    $wrapper->sam_to_fixed_sorted_bam($in_sam, $tmp_bam);
    unless ($wrapper->run_status() >= 1) {
        $self->warn("Failed the initial sam->bam step");
        return 0;
    }
    
    # bam -> fillmd'd bam
    my $tmp2_bam = $io->catfile($temp_dir, 'fillmd.bam');
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
    
    my $lane_path = delete $args{lane_path};
    if ($lane_path) {
        my @fastq_info = @{HierarchyUtilities::getFastqInfo($lane_path)};
        if ($fastq_info[0]->[0] && ! $fastq_info[1]->[0] && ! $fastq_info[2]->[0]) {
            $command = 'rmdupse';
        }
    }
    if (delete $args{single_ended}) {
        $command = 'rmdupse';
    }
    
    my $wrapper = VertRes::Wrapper::samtools->new(verbose => $self->verbose, %args);
    $wrapper->$command($in_bam, $out_bam);
    
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
        return symlink($in_bams[0], $out_bam);
    }
    
    my $io = VertRes::IO->new();
    my $verbose = $self->verbose();
    my $wrapper = VertRes::Wrapper::picard->new(verbose => $verbose,
                                                quiet => $verbose ? 0 : 1,
                                                validation_stringency => 'silent',
                                                tmp_dir => $io->tempdir());
    
    $wrapper->merge_and_check($out_bam, \@in_bams);
    
    return $wrapper->run_status() >= 1;
}

=head2 stats

 Title   : stats
 Usage   : $obj->stats('in.bam', 'in2.bam', ...);
 Function: Generate both flagstat and bas files for the given bam files.
 Returns : boolean (true on success)
 Args    : list of bam files (for each, a .bamstat and .flagstat file will be
           created)

=cut

sub stats {
    my ($self, @in_bams) = @_;
    
    my $wrapper = VertRes::Wrapper::samtools->new(verbose => $self->verbose);
    my $all_ok = 1;
    foreach my $in_bam (@in_bams) {
        my $flagstat = $in_bam.'.flagstat';
        $wrapper->flagstat($in_bam, $flagstat);
        unless ($wrapper->run_status() >= 1) {
            $all_ok = 0;
            unlink($flagstat);
        }
        
        my $bas = $in_bam.'.bas';
        my $ok = $self->bas($in_bam, $bas);
        unless ($ok) {
            $all_ok = 0;
            unlink($bas);
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
           are accurate NM tags for each record), output filename

=cut

sub bas {
    my ($self, $in_bam, $out_bas) = @_;
    
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
    
    # print stats per-readgroup
    my %readgroup_data = $self->bam_statistics($in_bam);
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
 Usage   : $obj->bam_statistics('in.bam');
 Function: Calculate all the stats per-readgroup for a bam needed for the bas
           format.
 Returns : hash with keys as readgroup ids and values as hash refs. Those refs
           have the keys:
           dcc_filename, md5, study, sample, platform, library, total_bases,
           mapped_bases, total_reads, mapped_reads, mapped_reads_paired_in_seq,
           mapped_reads_properly_paired, percent_mismatch, avg_qual, avg_isize,
           sd_isize, median_isize, mad
 Args    : input bam file (must have been run via samtools fillmd so that there
           are accurate NM tags for each record), optionally a sequence.index
           file from the DCC if working on bams with poor headers

=cut

sub bam_statistics {
    my ($self, $bam, $seq_index) = @_;
    
    my $sip = VertRes::Parser::sequence_index->new(file => $seq_index) if $seq_index;
    
    my $hu = VertRes::Utils::Hierarchy->new(verbose => $self->verbose);
    my ($dcc_filename, $study, $sample, $technology, $previous_rg) = $hu->dcc_filename($bam);
    my $md5;
    my $md5_file = $bam.'.md5';
    if (-s $md5_file) {
        open(my $md5fh, $md5_file) || $self->throw("Couldn't open md5 file '$md5_file': $!");
        my $line = <$md5fh>;
        ($md5) = split(' ', $line);
    }
    else {
        my $dmd5 = Digest::MD5->new();
        open(FILE, $bam) or $self->throw("Couldn't open bam '$bam': $!");
        binmode(FILE);
        $dmd5->addfile(*FILE);
        $md5 = $dmd5->hexdigest;
    }
    
    # view the bam and accumulate all the raw stats in little memory
    my $stw = VertRes::Wrapper::samtools->new(quiet => 1);
    $stw->run_method('open');
    my $view_fh = $stw->view($bam, undef, h => 1);
    $view_fh || $self->throw("Failed to samtools view '$bam'");
    
    my $ps = VertRes::Parser::sam->new(fh => $view_fh);
    my $rh = $ps->result_holder;
    
    my $fqu = VertRes::Utils::FastQ->new();
    
    my %readgroup_data;
    while ($ps->next_result) {
        my $rg = $rh->{RG};
        unless ($rg) {
            $self->warn("$rh->{QNAME} had no RG tag, using previous RG tag '$previous_rg'");
            $rg = $previous_rg;
        }
        $previous_rg = $rg;
        
        $readgroup_data{$rg}->{total_reads}++;
        my $read_length = length($rh->{SEQ});
        $readgroup_data{$rg}->{total_bases} += $read_length;
        
        my $flag = $rh->{FLAG};
        if ($ps->is_mapped($flag)) {
            $readgroup_data{$rg}->{mapped_reads}++;
            $readgroup_data{$rg}->{mapped_bases} += $read_length;
            $readgroup_data{$rg}->{mapped_reads_paired_in_seq}++ if $ps->is_sequencing_paired($flag);
            
            # avg quality of mapped bases
            foreach my $qual ($fqu->qual_to_ints($rh->{QUAL})) {
                $readgroup_data{$rg}->{qual_count}++;
                
                if ($readgroup_data{$rg}->{qual_count} == 1) {
                    $readgroup_data{$rg}->{qual_mean} = $qual;
                }
                else {
                    $readgroup_data{$rg}->{qual_mean} += ($qual - $readgroup_data{$rg}->{qual_mean}) / $readgroup_data{$rg}->{qual_count};
                }
            }
            
            # avg insert size and keep track of 's' for later calculation of sd.
            # algorithm based on http://www.johndcook.com/standard_deviation.html
            if ($ps->is_mapped_paired($flag)) {
                $readgroup_data{$rg}->{mapped_reads_properly_paired}++;
                
                if ($rh->{MAPQ} > 0) {
                    my $isize = $rh->{ISIZE} || 0;
                    if ($isize > 0) { # avoids counting the isize twice for a pair, since one will be negative
                        $readgroup_data{$rg}->{isize_count}++;
                        
                        if ($readgroup_data{$rg}->{isize_count} == 1) {
                            $readgroup_data{$rg}->{isize_mean} = $isize;
                            $readgroup_data{$rg}->{isize_s} = 0;
                        }
                        else {
                            my $old_mean = $readgroup_data{$rg}->{isize_mean};
                            $readgroup_data{$rg}->{isize_mean} += ($isize - $old_mean) / $readgroup_data{$rg}->{isize_count};
                            $readgroup_data{$rg}->{isize_s} += ($isize - $old_mean) * ($isize - $readgroup_data{$rg}->{isize_mean});
                        }
                        
                        # also, median insert size. Couldn't find an accurate
                        # running algorithm, but just keeping a histogram is
                        # accurate and uses less than 1MB. We can use the same
                        # histogram later to calculate the MAD.
                        $readgroup_data{$rg}->{isize_median_histogram}->{$isize}++;
                    }
                }
            }
            
            # for later calculation of mismatch %
            if (defined $rh->{NM}) {
                $readgroup_data{$rg}->{total_bases_with_nm} += $read_length;
                $readgroup_data{$rg}->{edits} += $rh->{NM};
            }
        }
    }
    
    # calculate the means etc.
    my $math_util = VertRes::Utils::Math->new();
    foreach my $rg (sort keys %readgroup_data) {
        my %data = %{$readgroup_data{$rg}};
        
        # calculate/round stats
        my $avg_qual = defined $data{qual_count} ? sprintf("%0.2f", $data{qual_mean}) : 0;
        my $avg_isize = defined $data{isize_count} ? sprintf("%0.0f", $data{isize_mean}) : 0;
        my $sd_isize = $avg_isize ? sprintf("%0.2f", sqrt($data{isize_s} / $data{isize_count})) : 0;
        my $percent_mismatch = defined $data{edits} ? sprintf("%0.2f", (100 / $data{total_bases_with_nm}) * $data{edits}) : 0;
        
        my $median_isize = 0;
        my $mad = 0;
        if (defined $data{isize_median_histogram}) {
            $median_isize = $math_util->histogram_median($data{isize_median_histogram});
            
            my %ads;
            while (my ($isize, $freq) = each %{$data{isize_median_histogram}}) {
                my $ad = abs($median_isize - $isize);
                $ads{$ad} += $freq;
            }
            
            $mad = $math_util->histogram_median(\%ads);
        }
        
        delete $readgroup_data{$rg}->{qual_count};
        delete $readgroup_data{$rg}->{qual_mean};
        delete $readgroup_data{$rg}->{isize_count};
        delete $readgroup_data{$rg}->{isize_mean};
        delete $readgroup_data{$rg}->{isize_s};
        delete $readgroup_data{$rg}->{edits};
        delete $readgroup_data{$rg}->{total_bases_with_nm};
        delete $readgroup_data{$rg}->{isize_median_histogram};
        
        $readgroup_data{$rg}->{avg_qual} = $avg_qual;
        $readgroup_data{$rg}->{avg_isize} = $avg_isize;
        $readgroup_data{$rg}->{sd_isize} = $sd_isize;
        $readgroup_data{$rg}->{percent_mismatch} = $percent_mismatch;
        $readgroup_data{$rg}->{median_isize} = $median_isize;
        $readgroup_data{$rg}->{mad} = $mad;
        
        # add in header info
        $readgroup_data{$rg}->{dcc_filename} = $dcc_filename;
        $readgroup_data{$rg}->{study} = $study;
        $readgroup_data{$rg}->{sample} = $sample;
        $readgroup_data{$rg}->{platform} = $ps->readgroup_info($rg, 'PL') || 'unknown_platform';
        $readgroup_data{$rg}->{library} = $ps->readgroup_info($rg, 'LB') || 'unknown_library';
        
        if ($readgroup_data{$rg}->{library} eq 'unknown_library' && $sip) {
            $readgroup_data{$rg}->{platform} = $sip->lane_info($rg, 'INSTRUMENT_PLATFORM') || 'unknown_platform';
            readgroup_data{$rg}->{library} = $sip->lane_info($rg, 'LIBRARY_NAME') || 'unknown_library';
        }
        
        $readgroup_data{$rg}->{md5} = $md5;
    }
    
    return %readgroup_data;
}

=head2 make_unmapped_bam

 Title   : make_unmapped_bam
 Usage   : $obj->make_unmapped_bam('in.bam', 'out.bam');
 Function: Given a bam file, generates another bam file that contains only the
           unmapped reads.
 Returns : boolean (true on success)
 Args    : starting sam file, output name for bam file.

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
    
    while (<$in_fh>) {
        if (/^\@/) {
            print $out_fh $_;
        }
        else {
            my (undef, $flag) = split;
            # 0x0004 is the unmapped flag
            if ($flag & 0x0004) {
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

1;
