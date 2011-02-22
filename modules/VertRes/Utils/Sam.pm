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
petr.danecek@sanger
Martin Hunt: mh12@sanger.ac.uk

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
use VertRes::Wrapper::fastqcheck;
use HierarchyUtilities;
use VertRes::Parser::dict;
use VertRes::Parser::sequence_index;
use VertRes::Parser::sam;
use VertRes::Parser::bam;
use VertRes::Parser::dict;
use VertRes::Parser::fasta;
use VertRes::Utils::FastQ;
use VertRes::Utils::Math;
use VertRes::Utils::Hierarchy;
use Digest::MD5;
use VertRes::Parser::bamcheck;
use VertRes::Parser::fastqcheck;
use List::Util qw(min max sum);
use Test::Deep::NoTest;
use Graphs;
use Data::Dumper;

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

my %tech_to_platform = (SLX => 'ILLUMINA',
                        '454' => 'LS454',
                        SOLID => 'ABI_SOLID',
                        ILLUMINA => 'ILLUMINA',
                        'LS454' => 'LS454',
                        ABI_SOLID => 'ABI_SOLID',);

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

=head2 num_bam_records

 Title   : num_bam_records
 Usage   : my $num_records = $obj->num_bam_records($bam_file);
 Function: Find the number of records (reads) in a bam file.
 Returns : int
 Args    : bam filename, optionally these hash options:
           only => 'regex' only count sequences that match the regex.

=cut

sub num_bam_records {
    my ($self, $bam_file, %opts) = @_;
    my $pars = VertRes::Parser::bam->new(file => $bam_file);
    my $rh = $pars->result_holder;
    my $rname_regex;
    if ($opts{only}) {
        $pars->get_fields('RNAME');
        $rname_regex = $opts{only};
    }
    
    my $records = 0;
    while ($pars->next_result) {
        if ($rname_regex) {
            $records++ if $rh->{RNAME} =~ /$rname_regex/;
        }
        else {
            $records++;
        }
    }
    return $records;
}

=head2 num_bam_header_lines

 Title   : num_bam_header_lines
 Usage   : my $num_lines = $obj->num_bam_header_lines($bam_file);
 Function: Find the number of header lines in a bam file.
 Returns : int
 Args    : bam filename

=cut

sub num_bam_header_lines {
    my ($self, $bam_file) = @_;
    
    my $st = VertRes::Wrapper::samtools->new(quiet => 1, run_method => 'open');
    my $fh = $st->view($bam_file, undef, H => 1);
    my $header_lines = 0;
    while (<$fh>) {
        $header_lines++;
    }
    close($fh);
    
    return $header_lines;
}
=head2 num_bam_lines

 Title   : num_bam_lines
 Usage   : my $num_lines = $obj->num_bam_lines($bam_file);
 Function: Find the number of records (reads) + header lines in a bam file.
 Returns : int
 Args    : bam filename

=cut

sub num_bam_lines {
    my ($self, $bam_file) = @_;
    my $header_lines = $self->num_bam_header_lines($bam_file);
    my $records = $self->num_bam_records($bam_file);
    return $header_lines + $records;
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
           merge => {'regex' => 'prefix'} to put sequences with ids that match
                     the regex into a single file with the corresponding prefix.
                     ignore takes precendence.
           non_chr => boolean (default true; ignore takes precendence: does a
                               merge for sequences that don't look chromosomal
                               into a 'nonchrom' prefixed file; does not make
                               individual bams for each non-chromosomsal seq if
                               ignore has not been set)
           make_unmapped => boolean (default true: a file prefixed with
                                     'unmapped' is made containing all the
                                     unmapped reads. NB: only unmapped pairs or
                                     singletons appear in this file; the other
                                     split bam files will contain unmapped reads
                                     where that read's mate was mapped)
           all_unmapped => boolean (default false; when true, the described
                                    behaviour of make_unmapped above changes so
                                    that the unmapped file contains all unmapped
                                    reads, potentially duplicating reads in
                                    different split files)
           only => 'regex' to only makes splits for sequences that match the
                    regex. This changes the default of make_unmapped to false,
                    but you can turn it back on explicitly. non_chr is disabled.
           output_dir => 'path' to specify where the split bams are created;
                         default is the same dir as the input bam
           pretend => boolean (if true, don't actually do anything, just return
                               what files would be made)

=cut

sub split_bam_by_sequence {
    my ($self, $bam, %opts) = @_;
    
    my $basename = basename($bam);
    my $output_dir;
    if (defined $opts{output_dir}) {
        $output_dir = $opts{output_dir};
        $self->throw("output_dir '$output_dir' does not exist") unless -d $output_dir;
    }
    else {
        $output_dir = $bam;
        $output_dir =~ s/$basename$//;
    }
    unless (defined $opts{make_unmapped}) {
        $opts{make_unmapped} = $opts{only} ? 0 : 1;
    }
    unless (defined $opts{non_chr}) {
        $opts{non_chr} = 1;
    }
    if ($opts{non_chr}) {
        $opts{merge}->{'^(?:N[TC]_\d+|GL\d+)'} = 'nonchrom';
        unless (defined $opts{ignore}) {
            $opts{ignore} = '^(?:N[TC]_\d+|GL\d+)';
        }
    }
    if ($opts{only}) {
        $opts{non_chr} = 0;
    }
    
    # find out what sequences there are
    my $pb = VertRes::Parser::bam->new(file => $bam);
    my %all_sequences = $pb->sequence_info();
    
    # work out the output filenames for each sequence
    my %seq_to_bam;
    my %output_bams;
    foreach my $seq (keys %all_sequences) {
        my @prefixes;
        if (defined $opts{merge}) {
            while (my ($regex, $this_prefix) = each %{$opts{merge}}) {
                if ($seq =~ /$regex/) {
                    push(@prefixes, $this_prefix);
                }
            }
        }
        unless (scalar(@prefixes) || ($opts{ignore} && $seq =~ /$opts{ignore}/)) {
            @prefixes = ('chrom'.$seq);
        }
        
        if ($opts{only}) {
            next unless $seq =~ /$opts{only}/;
        }
        
        foreach my $prefix (@prefixes) {
            my $out_bam = File::Spec->catfile($output_dir, $prefix.'.'.$basename);
            $output_bams{$out_bam} = 1;
            next if -s $out_bam;
            push(@{$seq_to_bam{$seq}}, $out_bam);
        }
    }
    
    my @get_fields = ('RNAME');
    my ($unmapped_bam, $skip_mate_mapped);
    if ($opts{make_unmapped}) {
        $unmapped_bam = File::Spec->catfile($output_dir, 'unmapped.'.$basename);
        $output_bams{$unmapped_bam} = 1;
        push(@get_fields, 'FLAG');
        $skip_mate_mapped = $opts{all_unmapped} ? 0 : 1;
        if (-s $unmapped_bam) {
            undef $unmapped_bam;
        }
    }
    
    if ($opts{pretend}) {
        my @outs = sort keys %output_bams;
        return @outs;
    }
    
    # stream through the bam, outputting all our desired files
    $pb->get_fields(@get_fields);
    my $rh = $pb->result_holder;
    my %counts;
    while ($pb->next_result) {
        my $out_bams = $seq_to_bam{$rh->{RNAME}};
        if ($out_bams) {
            foreach my $out_bam (@{$out_bams}) {
                $pb->write_result($out_bam.'.unchecked');
                $counts{$out_bam}++;
            }
        }
        
        if ($unmapped_bam) {
            my $flag = $rh->{FLAG};
            unless ($pb->is_mapped($flag)) {
                if ($skip_mate_mapped && $pb->is_sequencing_paired($flag)) {
                    next if $pb->is_mate_mapped($flag);
                }
                $pb->write_result($unmapped_bam.'.unchecked');
                $counts{$unmapped_bam}++;
            }
        }
    }
    $pb->close;
    
    # check all the bams we created
    my @outs;
    foreach my $out_bam (sort keys %output_bams) {
        if (-s $out_bam) {
            push(@outs, $out_bam);
            next;
        }
        
        my $unchecked = $out_bam.'.unchecked';
        # the input bam might not have had reads mapped to every sequence in the
        # header, so we might not have created all the bams expected. Create
        # a header-only bam in that case:
        unless (-s $unchecked) {
            $pb = VertRes::Parser::bam->new(file => $bam);
            $pb->next_result;
            $pb->_create_no_record_output_bam($unchecked);
            $pb->close;
        }
        
        my $actual_reads = $self->num_bam_records($unchecked);
        my $expected_reads = $counts{$out_bam} || 0;
        
        if ($expected_reads == $actual_reads) {
            move($unchecked, $out_bam) || $self->throw("Couldn't move $unchecked to $out_bam");
            push(@outs, $out_bam);
        }
        else {
            $self->warn("When attempting $bam -> $out_bam split, got $actual_reads reads instead of $expected_reads; will delete it");
            unlink($unchecked);
        }
    }
    
    return @outs;
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
           platform => string (Must be sanger or DCC standard naming: SLX|SOLEXA
                               454|LS454 SOLID|ABI_SOLID)
           centre => string
           insert_size => int
           study => string, the study id, eg. SRP000001
 
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
    $platform = $tech_to_platform{$platform} || $self->throw("Bad platform '$platform'");
    my $centre = $args{centre} || $parser->lane_info($lane, 'CENTER_NAME');
    my $project = $args{study} || $args{project} || $parser->lane_info($lane, 'STUDY_ID');
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
        my $sp = $rh->{SP} || '';
        $sp = "\t$sp" if $sp;
        my $local_ref_name = $rh->{AS} || $ref_name;
        my $local_ur = $rh->{UR} || $ref_fa;
        print $shfh "\@SQ\tSN:$rh->{SN}\tLN:$rh->{LN}\tAS:$local_ref_name\tUR:$local_ur\tM5:$rh->{M5}$sp\n";
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
    
    my $tmp3_bam = "$out_bam.copied";
    unlink($tmp3_bam);
    move($tmp2_bam, $tmp3_bam) || $self->throw("Failed to move $tmp2_bam to $tmp3_bam: $!");
    
    # check it
    $wrapper->run_method('open');
    my $fh = $wrapper->view($tmp3_bam, undef, h => 1);
    my $bam_count = 0;
    while (<$fh>) {
        $bam_count++;
    }
    close($fh);
    $io->file($in_sam);
    my $sam_count = $io->num_lines();
    if ($bam_count >= $sam_count) {
        move($tmp3_bam, $out_bam) || $self->throw("Failed to move $tmp2_bam to $out_bam: $!");
        return 1;
    }
    else {
        $self->warn("$tmp3_bam is bad (only $bam_count lines vs $sam_count)");
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

=head2 markdup

 Title   : markdup
 Usage   : $obj->markdup('in.bam', 'dupsmarked.bam');
 Function: Marks duplicates in a bam file, automatically taking into account
           if these are single ended or paired reads.
 Returns : boolean (true on success)
 Args    : starting bam file, output name for bam file

           Optionally, args to pass to VertRes::Wrapper::picard (eg.
           quiet => 1)

=cut

sub markdup {
    my ($self, $in_bam, $out_bam, %args) = @_;
    
    my $verbose = $self->verbose();
    
    # picard needs a tmp dir, but we don't use /tmp because it's likely to fill
    # up
    my (undef, $path) = fileparse($out_bam);
    my $fsu = VertRes::Utils::FileSystem->new();
    my $tmp_dir = $fsu->tempdir('_markdup_tmp_XXXXXX', DIR => $path);
    
    my $wrapper = VertRes::Wrapper::picard->new(verbose => $verbose,
                                                quiet => $verbose ? 0 : 1,
                                                $self->{java_memory} ? (java_memory => $self->{java_memory}) : (),
                                                tmp_dir => $tmp_dir,
                                                %args);
    
    $wrapper->markdup($in_bam, $out_bam);
    
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
    
    my ($basename, $path) = fileparse($out_bam);
    if (@in_bams == 1) {
        # we want a relative symlink, so use the path of out_bam and make the
        # symlink relative to that
        my $rel_path = File::Spec->abs2rel($in_bams[0], $path);
        return symlink($rel_path, $out_bam);
    }
    
    # picard needs a tmp dir for merging, but we don't use /tmp because it's
    # likely to fill up
    my $fsu = VertRes::Utils::FileSystem->new();
    my $tmp_dir = $fsu->tempdir('_merge_tmp_XXXXXX', DIR => $path);
    
    my $verbose = $self->verbose();
    my $wrapper = VertRes::Wrapper::picard->new(verbose => $verbose,
                                                quiet => $verbose ? 0 : 1,
                                                $self->{java_memory} ? (java_memory => $self->{java_memory}) : (),
                                                validation_stringency => 'silent',
                                                tmp_dir => $tmp_dir);
    
    $wrapper->merge_and_check($out_bam, \@in_bams);
    
    return $wrapper->run_status() >= 1;
}

=head2 stats

 Title   : stats
 Usage   : $obj->stats('YYYYMMDD', 'in.bam', 'in2.bam', ...);
 Function: Generate both flagstat and bas files for the given bam files.
 Returns : boolean (true on success)
 Args    : YYYYMMDD format date string to signify the release date, list of bam
           files (for each, a .bas and .flagstat file will be created)

=cut

sub stats {
    my ($self, $release_date, @in_bams) = @_;
    $release_date =~ /^\d{8}$/ || $self->throw("bad release date '$release_date'");
    
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
        my $ok = $self->bas($in_bam, $release_date, $bas.'.tmp');
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

=head2 index_bams

 Title   : index_bams
 Usage   : $obj->index_bams(files=>['file.bam',...]);
           $obj->index_bams(fofn=>'bam_files.list');
 Function: Index one or more BAM files, optionally reading the file names from a file.
 Returns : boolean (true on success)
 Args    : bam file names or file with a list of bam files

=cut

sub index_bams 
{
    my ($self, %args) = @_;
    my @bams;
    if ( exists($args{files}) ) { @bams = @{$args{files}} }
    if ( exists($args{fofn}) )
    {
        open(my $fh,'<',$args{fofn}) or $self->throw("$args{fofn}: $!");
        while (my $line=<$fh>)
        {
            chomp($line);
            push @bams,$line;
        }
        close($fh);
    }
    my $samtools = VertRes::Wrapper::samtools->new(verbose => $self->verbose, quiet => 1);
    for my $bam (@bams)
    {
        $samtools->run_method('system');
        $samtools->index($bam, $bam.'.bai');
        $samtools->run_status >= 1 || $self->throw("Failed to create $bam.bai");
    }
    return 1;
}


=head2 bas

 Title   : bas
 Usage   : $obj->bas('in.bam', 'YYYYMMDD', 'out.bas', 'sequence.index');
 Function: Generate a 'bas' file. These are bam statistic files that provide
           various stats about each readgroup in a bam.
 Returns : boolean (true on success)
 Args    : input bam file (must have been run via samtools fillmd so that there
           are accurate NM tags for each record), YYYYMMDD string describing the
           date of the release, output filename, sequence.index file (optional
           if your bam headers are good and you don't care about column 1 of the
           bas file, which is the DCC filename)

=cut

sub bas {
    my ($self, $in_bam, $release_date, $out_bas, $seq_index) = @_;
    $release_date =~ /^\d{8}$/ || $self->throw("bad release date '$release_date'");
    
    my $working_bas = $out_bas.'.working';
    open(my $bas_fh, '>', $working_bas) || $self->throw("Couldn't write to '$working_bas': $!");
    
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
                             'insert_size_median_absolute_deviation',
                             '#_duplicate_reads'), "\n";
    $expected_lines++;
    
    # get the stats for each read group
    my %readgroup_data = $self->bam_statistics($in_bam);
    
    # get the meta data
    my $hu = VertRes::Utils::Hierarchy->new(verbose => $self->verbose);
    my $dcc_filename = $hu->dcc_filename($in_bam, $release_date, $seq_index);
    
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
    
    my $pb = VertRes::Parser::bam->new(file => $in_bam);
    
    # add in the meta data
    while (my ($rg, $data) = each %readgroup_data) {
        $readgroup_data{$rg}->{dcc_filename} = $dcc_filename;
        $readgroup_data{$rg}->{study} = $pb->readgroup_info($rg, 'DS') || 'unknown_study';
        $readgroup_data{$rg}->{sample} = $pb->readgroup_info($rg, 'SM') || 'unknown_sample';
        $readgroup_data{$rg}->{platform} = $pb->readgroup_info($rg, 'PL') || 'unknown_platform';
        $readgroup_data{$rg}->{library} = $pb->readgroup_info($rg, 'LB') || 'unknown_library';
        
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
        
        if (defined $tech_to_platform{$readgroup_data{$rg}->{platform}}) {
            $readgroup_data{$rg}->{platform} = $tech_to_platform{$readgroup_data{$rg}->{platform}};
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
                                 $data{mad},
                                 $data{duplicate_reads} || 0), "\n";
        $expected_lines++;
    }
    close($bas_fh);
    
    my $io = VertRes::IO->new(file => $working_bas);
    my $actual_lines = $io->num_lines;
    if ($actual_lines == $expected_lines) {
        move($working_bas, $out_bas) || $self->throw("Could not move $working_bas to $out_bas");
        return 1;
    }
    else {
        unlink($working_bas);
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
    my $pb = VertRes::Parser::bam->new(file => $bam);
    $pb->get_fields('SEQ_LENGTH', 'MAPPED_SEQ_LENGTH', 'FLAG', 'QUAL', 'MAPQ', 'ISIZE', 'RG', 'NM');
    my $rh = $pb->result_holder;
    
    my $fqu = VertRes::Utils::FastQ->new();
    
    my %readgroup_data;
    my $previous_rg = 'unknown_readgroup';
    while ($pb->next_result) {
        my $rg = $rh->{RG};
        my $flag = $rh->{FLAG};
        
        unless ($rg) {
            $self->warn("a read had no RG tag, using previous RG tag '$previous_rg'");
            $rg = $previous_rg;
        }
        $previous_rg = $rg;
        
        my @this_rg_data = @{$readgroup_data{$rg} || []};
        $this_rg_data[0]++;
        $this_rg_data[1] += $rh->{SEQ_LENGTH};
        
        if ($pb->is_mapped($flag)) {
            $this_rg_data[2]++;
            $this_rg_data[3] += $rh->{MAPPED_SEQ_LENGTH};
            $this_rg_data[4]++ if $pb->is_sequencing_paired($flag);
            
            # avg quality of mapped bases
            foreach my $qual ($fqu->qual_to_ints($rh->{QUAL})) {
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
            if ($pb->is_mapped_paired($flag)) {
                $this_rg_data[7]++;
                
                my $isize = $rh->{ISIZE} || 0;
                
                if ($rh->{MAPQ} > 0) {
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
            my $nm = $rh->{NM};
            if (defined $nm && $nm ne '*') {
                $this_rg_data[12] += $rh->{SEQ_LENGTH};
                $this_rg_data[13] += $nm;
            }
        }
        
        if ($pb->is_duplicate($flag)) {
            $this_rg_data[14]++;
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
        $rg_stats{duplicate_reads} = $data[14] || 0;
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
 Args    : starting bam file, output name for bam file. Optional boolean value
           which if true will not output unmapped reads if that read's mate was
           mapped

=cut

sub make_unmapped_bam {
    my ($self, $in_bam, $out_bam, $skip_mate_mapped) = @_;
    
    my $pb = VertRes::Parser::bam->new(file => $in_bam);
    $pb->flag_selector(self_unmapped => 1);
    if ($skip_mate_mapped) {
        $pb->get_fields('FLAG');
    }
    my $rh = $pb->result_holder;
    
    my $working_bam = $out_bam.'.working';
    my $count = 0;
    while ($pb->next_result) {
        if ($skip_mate_mapped) {
            my $flag = $rh->{FLAG};
            if ($pb->is_sequencing_paired($flag)) {
                next if $pb->is_mate_mapped($flag);
            }
        }
        $pb->write_result($working_bam);
        $count++;
    }
    $pb->close;
    
    my $actual = $self->num_bam_records($working_bam);
    
    if ($actual == $count) {
        move($working_bam, $out_bam) || $self->throw("Failed to move $working_bam to $out_bam");
        return 1;
    }
    else {
        $self->warn("Created an output bam, but it had $actual instead of $count records; will unlink it");
        unlink($working_bam);
        return 0;
    }
}

=head2 add_unmapped

 Title   : add_unmapped
 Usage   : $obj->add_unmapped($sam_file, @fastqs);
 Function: Append reads in the fastqs not already present in the sam file to the
           end of the sam file.
 Returns : boolean
 Args    : sam file name, fastq file name(s) (NB: must be supplied in the order
           first of a pair, second of a pair, if supplying a pair of fastqs.
           Also, these fastqs are expected to have read names ending in /1 and
           /2 respectively)

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
    my $snames_uniq_file = $snames_file.'.uniq';
    my $s_lines = 0;
    unless (-s $snames_uniq_file) {
        unless (-s $snames_sorted_fixed_file) {
            my $dubious = 0;
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
                                $dubious = 1;
                            }
                            
                            $qname =~ s/\/$//;
                            $qname .= "/$read_num";
                        }
                        
                        $s_lines++;
                        print $rfh $qname, "\n";
                    }
                    close($sfh);
                    close($rfh);
                    
                    # check rnames file isn't truncated
                    my $actual_count = VertRes::IO->new(file => $snames_file)->num_lines;
                    unless ($actual_count == $s_lines) {
                        unlink($snames_file);
                        $self->throw("made an rnames file but it was truncated!");
                    }
                }
                
                $s_lines ||= VertRes::IO->new(file => $snames_file)->num_lines;
                
                system("sort $snames_file > $snames_sorted_file");
                my $actual_count = VertRes::IO->new(file => $snames_sorted_file)->num_lines;
                unless ($actual_count == $s_lines) {
                    unlink($snames_sorted_file);
                    $self->throw("made an rnames_sorted_file file but it was truncated!");
                }
                
                unlink($snames_file);
            }
            
            if ($dubious) {
                # if we have the same read name twice in a row, change the second one
                # to be the second read of a pair
                open(my $ofh, '>', $snames_sorted_fixed_file) || $self->throw("Could not write to $snames_sorted_fixed_file");
                open(my $ifh, $snames_sorted_file) || $self->throw("Could not open $snames_sorted_file");
                my $previous_name = '';
                $s_lines = 0;
                while (<$ifh>) {
                    $s_lines++;
                    
                    if ($_ eq $previous_name) {
                        if (/1\n$/) {
                            s/1\n$/2\n/;
                        }
                        else {
                            s/2\n$/1\n/;
                        }
                    }
                    
                    print $ofh $_;
                    
                    $previous_name = $_;
                }
                close($ofh);
                close($ifh);
                
                # check rnames_sorted_fixed_file file isn't truncated
                my $actual_count = VertRes::IO->new(file => $snames_sorted_fixed_file)->num_lines;
                unless ($actual_count == $s_lines) {
                    unlink($snames_sorted_fixed_file);
                    $self->throw("made an rnames_sorted_fixed_file file but it was truncated!");
                }
                
                unlink($snames_sorted_file);
            }
            else {
                move($snames_sorted_file, $snames_sorted_fixed_file) || $self->throw("Could not move $snames_sorted_file to $snames_sorted_fixed_file");
            }
        }
        
        # for some reason we can end up with multiple copies of some read names;
        # get rid of them with uniq. Unfortunately we have no way of knowning
        # how many lines to expect, so can't check for truncation :(
        my $failed = system("uniq -u $snames_sorted_fixed_file > $snames_uniq_file");
        if ($failed) {
            unlink($snames_uniq_file);
            $self->throw("failed to uniq $snames_sorted_fixed_file");
        }
        unlink($snames_sorted_fixed_file);
    }
    $s_lines = VertRes::IO->new(file => $snames_uniq_file)->num_lines;
    
    # now list out to another file all the read names in the fastqs
    my $qnames_file = $sam.'.qnames';
    my $q_lines = 0;
    unless (-s $qnames_file) {
        open(my $qfh, '>', $qnames_file) || $self->throw("Could not write to $qnames_file");
        foreach my $fq (@fqs) {
            my $fqp = VertRes::Parser::fastq->new(file => $fq);
            my $rh = $fqp->result_holder();
            
            # loop through all the sequences in the fastq (without indexing
            # to save memory)
            while ($fqp->next_result(1)) {
                $q_lines++;
                print $qfh $rh->[0], "\n";
            }
        }
        close($qfh);
        
        # check it's not truncated
        my $actual_count = VertRes::IO->new(file => $qnames_file)->num_lines;
        unless ($actual_count == $q_lines) {
            unlink($qnames_file);
            $self->throw("made a qnames_file file but it was truncated!");
        }
    }
    $q_lines ||= VertRes::IO->new(file => $qnames_file)->num_lines;
    
    my $unmapped = $q_lines - $s_lines;
    unless ($unmapped) {
        unlink($snames_uniq_file);
        unlink($qnames_file);
        $self->warn("no reads were unmapped, nothing to do!");
        return 1;
    }
    
    # Now cat the qnames and snames files, sort and use uniq on that to find
    # the names that weren't in both. (unix diff, or an attempt at finding a
    # diff in perl uses too much memory.)
    my $sq_file = $sam.'.sqnames';
    my $sq_sorted_file = $sq_file.'.sorted';
    my $sq_uniq_file = $sq_sorted_file.'.uniq';
    unless (-s $sq_uniq_file) {
        unless (-s $sq_sorted_file) {
            my $total_lines = $s_lines + $q_lines;
            unless (-s $sq_file) {
                system("cat $snames_uniq_file $qnames_file > $sq_file");
                my $actual_count = VertRes::IO->new(file => $sq_file)->num_lines;
                unless ($actual_count == $total_lines) {
                    unlink($sq_file);
                    $self->throw("made an sq_file file but it was truncated!");
                }
            }
            
            system("sort $sq_file > $sq_sorted_file");
            my $actual_count = VertRes::IO->new(file => $sq_sorted_file)->num_lines;
            unless ($actual_count == $total_lines) {
                unlink($sq_sorted_file);
                $self->throw("made an sq_sorted_file file but it was truncated!");
            }
            
            unlink($sq_file);
        }
        
        my $failed = system("uniq -u $sq_sorted_file > $sq_uniq_file");
        if ($failed) {
            unlink($snames_uniq_file);
            $self->throw("failed to uniq $snames_sorted_fixed_file");
        }
        
        my $actual_count = VertRes::IO->new(file => $sq_uniq_file)->num_lines;
        unless ($actual_count == $unmapped) {
            unlink($sq_uniq_file);
            $self->throw("made a sq_uniq_file file but it was truncated! ($actual_count lines vs $unmapped unmapped, from $q_lines q - $s_lines s");
        }
        
        unlink($sq_sorted_file);
    }
    
    # now for all the reads we're missing, append to original sam
    my @fqps;
    foreach my $fq (@fqs) {
        my $fqp = VertRes::Parser::fastq->new(file => $fq);
        my $rh = $fqp->result_holder();
        push(@fqps, [$fqp, $rh]);
    }
    open(my $sfh, '>>', $sam) || $self->throw("Could not append to $sam");
    open (my $ufh, $sq_uniq_file) || $self->throw("Could not open $sq_uniq_file");
    
    if (-s $sam < 999999999) {
        # we can hash the unmapped and go through the fastqs sequentially, and
        # we don't care about the order of sequences in the fastqs
        my %unmapped;
        while (<$ufh>) {
            chomp;
            $unmapped{$_} = 1;
        }
        
        foreach my $fqp_data (@fqps) {
            my ($fqp, $rh) = @{$fqp_data};
            
            while ($fqp->next_result(1)) {
                my $id = $rh->[0];
                next unless exists $unmapped{$id};
                
                my $seq = $rh->[1];
                my $qual = $rh->[2];
                
                my ($read_num) = $id =~ /\/([12])$/;
                my ($first, $second) = (0, 0, 0);
                my $mate_id = $id;
                if ($paired) {
                    if ($read_num == 2) {
                        $second = 1;
                        $mate_id =~ s/2$/1/;
                    }
                    else {
                        $first = 1;
                        $mate_id =~ s/1$/2/;
                    }
                }
                
                # work out the SAM flags
                my $flags = $self->calculate_flag(self_unmapped => 1,
                                                  paired_tech => $paired,
                                                  mate_unmapped => exists $unmapped{$mate_id},
                                                  '1st_in_pair' => $first,
                                                  '2nd_in_pair' => $second);
                
                print $sfh "$id\t$flags\t*\t0\t0\t*\t*\t0\t0\t$seq\t$qual\n";
            }
        }
    }
    else {
        # we must assume that our list of unmapped read names in $ufh is in the
        # same sort order as sequences in the fastqs, so for each unmapped read
        # name...
        my $previous_root = '';
        while (<$ufh>) {
            chomp;
            my $read_name = $_;
            my ($read_num) = $read_name =~ /\/([12])$/;
            my ($fqp, $rh);
            if ($read_num) {
                ($fqp, $rh) = @{$fqps[$read_num - 1]};
            }
            else {
                ($fqp, $rh) = @{$fqps[0]};
            }
            
            my ($mate_unmapped, $first, $second) = (0, 0, 0);
            my $root_name = $read_name;
            $root_name =~ s/\/[12]$//;
            if ($paired) {
                if ($read_num == 2) {
                    $second = 1;
                    
                    if ($previous_root eq $root_name) {
                        $mate_unmapped = 1;
                    }
                }
                else {
                    $first = 1;
                    
                    # see if the next line matches root name
                    my $tell = tell($ufh);
                    my $next_line = <$ufh>;
                    seek($ufh, $tell, 0);
                    my ($next_root) = $next_line =~ /^(.+)\/[12]\n$/;
                    if ($next_root eq $root_name) {
                        $mate_unmapped = 1;
                    }
                }
            }
            $previous_root = $root_name;
            
            # ... loop through the fastq until we find the corresponding sequence
            $fqp->_seek_first_result;
            my $found = 0;
            while ($fqp->next_result(1)) {
                my $id = $rh->[0];
                next unless $id eq $read_name;
                $found = 1;
                my $seq = $rh->[1];
                my $qual = $rh->[2];
                
                # work out the SAM flags
                my $flags = $self->calculate_flag(self_unmapped => 1,
                                                  paired_tech => $paired,
                                                  mate_unmapped => $mate_unmapped,
                                                  '1st_in_pair' => $first,
                                                  '2nd_in_pair' => $second);
                
                print $sfh "$id\t$flags\t*\t0\t0\t*\t*\t0\t0\t$seq\t$qual\n";
                
                last;
            }
            
            unless ($found) {
                $self->throw("Didn't find the unmapped read $read_name in the fastqs!");
            }
        }
    }
    close($sfh);
    close($ufh);
    
    unlink($snames_uniq_file);
    unlink($qnames_file);
    unlink($sq_uniq_file);
    
    my $lines_now = VertRes::IO->new(file => $sam)->num_lines;
    unless ($lines_now >= $q_lines) {
        $self->warn("Tried appending $unmapped unmapped records to the $s_lines mapped records sam file, but ended up with only $lines_now lines!");
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

=head2 header_rewrite_required

 Title   : header_rewrite_required
 Usage   : if ($obj->header_rewrite_required('a.bam', 'RG', 
                                    readgroup1 => { sample_name => 'NA000001' })
 Function: Prior to calling change_header_lines() with the same arguments, just
           check that a rewrite is even necessary.
 Returns : boolean - true if changes will be made, false otherwise.
 Args    : as per change_header_lines()

=cut

sub header_rewrite_required {
    my ($self, $bam, %changes) = @_;
    
    keys %changes || return 1;
    my %arg_to_tag = (RG => {sample_name => 'SM',
                      	library => 'LB',
                      	platform => 'PL',
                      	centre => 'CN',
                      	insert_size => 'PI',
                      	study => 'DS',
                      	project => 'DS'},
                      PG => {program => 'PN',
                      	command => 'CL',
                      	version => 'VN'},
                      SQ => {seq_name => 'SN',
                      	ref_length => 'LN',
                      	assembler => 'AS',
                      	md5 => 'M5',
                      	species => 'SP',
                      	uri => 'UR'});

	my $stin = VertRes::Wrapper::samtools->new(quiet => 1, run_method => 'open');

    my $remove_unique = 0;
    if (exists $changes{'PG'}{'remove_unique'}) {
    	$remove_unique = $changes{'PG'}{'remove_unique'};
    }

    my $dict;
    if (exists $changes{'SQ'}{'from_dict'}) {
		$dict = $changes{'SQ'}{'from_dict'};
    }
    
    # Compare the SQ lines of the bam and the dict file
	if ($dict) {
		my $bfh = $stin->view($bam, undef, H => 1);
    	$bfh || $self->throw("Could not read header from '$bam'");
    	my @bheader;
    	while (<$bfh>) {
    		chomp;
    		next unless /\@SQ/;
    		push @bheader, $_;
    	}
    	close $bfh;
    	
        open my $dfh, "<$dict" || $self->throw("Could not open dictionary, $dict");
    	my @dheader;
    	while (<$dfh>) {
    		chomp;
    		next unless /\@SQ/;
    		push @dheader, $_;
    	}
    	close $bfh;
    	
    	my $dict_header = join "\n", @dheader;
    	my $bam_header = join "\n", @bheader;
    	
    	return 1 unless ($dict_header eq $bam_header);
	}
	
    if (defined $changes{platform}) {
        $changes{platform} = $tech_to_platform{$changes{platform}} || $self->throw("Bad platform '$changes{platform}'");
    }
    
        
    # Parse original header and grab the lane id
    my $lane;
    if ($remove_unique) {
		my $pars = VertRes::Parser::sam->new(file => $bam);
		my %read_info = $pars->readgroup_info();
		my @lane_ids = keys %read_info;
		$self->throw("$bam does not have a single lane id.") unless (scalar @lane_ids == 1);
		$lane = $lane_ids[0];
    }
    
	# Write new header
    my $made_changes = 0;
	my $bamfh = $stin->view($bam, undef, H => 1);
    $bamfh || $self->throw("Could not read header from '$bam'");
	my $sq_from_dict_done = 0;
    while (<$bamfh>) {
    	foreach my $line_tag (keys %changes) {
    		my $id_tag = $line_tag eq 'SQ' ? 'SN' : 'ID';
	        if (/^\@$line_tag.+$id_tag:([^\t]+)/) {
	            my $id = $1;
	            if (exists $changes{$line_tag}{$id}) {
	                while (my ($arg, $value) = each %{$changes{$line_tag}{$id}}) {
	                    next unless defined $value;
	                    my $tag = $arg_to_tag{$line_tag}{$arg} || next;
	                    
	                    if (/\t$tag:([^\t\n]+)/) {
	                        if ($1 ne $value) {
	                            $made_changes = 1;
	                        }
	                    }
	                    else {
	                        $made_changes = 1;
	                    }
	                }
	            }

	            if (exists $changes{$line_tag}{ all }) {
	                while (my ($arg, $value) = each %{$changes{$line_tag}{all}}) {
	                    next unless defined $value;
	                    my $tag = $arg_to_tag{$line_tag}{$arg} || next;
	                    
	                    if (/\t$tag:([^\t\n]+)/) {
	                        if ($1 ne $value) {
	                            $made_changes = 1;
	                        }
	                    }
	                    else {
	                        $made_changes = 1;
	                    }
	                }
	            }
	            
	            if ($remove_unique) {
					if (/\tCL:([^\t]+)/) {
			        	if ($1 =~ /\@|$lane/) {
							$made_changes = 1;
						}
			        }
				}
			}
        }
    }
    close($bamfh);

    return $made_changes;
}


=head2 replace_bam_header

 Title   : replace_bam_header
 Usage   : $obj->replace_bam_header('a.bam', 'sam.header');
 Function: Replaces the header of a given bam with the given header in sam format.
 Returns : boolean (true on success, meaning a change was made and the bam has
           been confirmed OK, or no change was necessary and the bam was
           untouched)
 Args    : path to bam, path to new header.

=cut

sub replace_bam_header {
    my ($self, $bam, $new_header) = @_;
    
    # Setup temporary output bam
    my $temp_bam = $bam.'.reheader.tmp.bam';
    $self->register_for_unlinking($temp_bam);
	
	# Run samtools reheader
	my $sw = VertRes::Wrapper::samtools->new(run_method => 'system');
	$sw->reheader($new_header, $bam, $temp_bam);
		
    # Check for trunction in the number of bam records (headers may have different number of lines)
    my $old_bam_records = num_bam_records($bam);
    my $new_bam_records = num_bam_records($temp_bam);
	    
    if ( $old_bam_records == $new_bam_records ) {
        unlink($bam);
        move($temp_bam, $bam) || $self->throw("Failed to move $temp_bam to $bam: $!");
        return 1;
    }
    else {
        $self->warn("$temp_bam is bad, will unlink it");
        unlink($temp_bam);
        return 0;
    }
}

=head2 change_header_lines

 Title   : change_header_lines
 Usage   : $obj->change_header_lines('a.bam', (RG => readgroup1 => { sample_name => 'NA000001' }));
           $obj->change_header_lines('a.bam', (PG => 'PG ID' => { version => '1.0.34' }));
           $obj->change_header_lines('a.bam', (PG => { remove_unique => 1}));
           $obj->change_header_lines('a.bam', (SQ => { from_dict => '/path/to/dict/'}));
           $obj->change_header_lines('a.bam', (SQ => { all => { uri => 'ftp://something.awesome.com/huzzah.fasta'}}));
           $obj->change_header_lines('a.bam', (PG => {'PG ID' => { command => {'foo' => 'bar'} })); # NOT YET IMPLEMENTED
 Function: Changes @RG or @PG or @SQ header lines as specified. Option 'remove_unique => 1' for @PG 
           tags will remove any labels which contain lane specific information or which 
           contain random elements generated by GATK (which look like label=something@random).
           Option 'from_dict => 'path/to/dict' for @SQ will load @SQ as specified in a dict file.
           Using all => {tag => value} will do the replacement for all not, sot just a specified
           readgroup or id.
           Only intended to be run on lane level bams. 
           ** Deprecates functionality of rewrite_bam_header and standardise_pg_header_lines. 
 Returns : boolean (true on success, meaning a change was made and the bam has
           been confirmed OK, or no change was necessary and the bam was
           untouched)
 Args    : path to bam, string 'RG' or 'PG' according to which types of header lines you 
           want to edit and a hash of changes to be made. Keys are readgroup or program 
           identifiers and values are hash refs with at least one of the following:
           For @RG changes:
           	sample_name => string, eg. 'NA00000';
           	library => string
           	platform => string
           	centre => string
           	insert_size => int
           	study => string, the study id, eg. 'SRP000001' (overwrites DS)
           For @PG changes:
           	version => string, eg. '0.5.5';
           	command => string, eg. '-q 0.9 out=some_file
           	program => string, eg. 'bwa' (overwrites PN tag)
           For @PG, there is the added option 'remove_unique => 1' which will 
           perform as specified above.
           For @SQ changes:
           	seq_name => string
           	ref_length => int
           	assembler => string
           	md5 => string
           	species => string
           	uri => string
           For @SQ, there is the added option 'dict => /path/to/dict' which will 
           perform as specified above.

=cut

sub change_header_lines {
    my ($self, $bam, %changes) = @_;
    
    keys %changes || return 1;
    my %arg_to_tag = (RG => {sample_name => 'SM',
                      	library => 'LB',
                      	platform => 'PL',
                      	centre => 'CN',
                      	insert_size => 'PI',
                      	study => 'DS',
                      	project => 'DS'},
                      PG => {program => 'PN',
                      	command => 'CL',
                      	version => 'VN'},
                      SQ => {seq_name => 'SN',
                      	ref_length => 'LN',
                      	assembler => 'AS',
                      	md5 => 'M5',
                      	species => 'SP',
                      	uri => 'UR'});

	my $stin = VertRes::Wrapper::samtools->new(quiet => 1, run_method => 'open');

    my $remove_unique = 0;
    if (exists $changes{'PG'}{'remove_unique'}) {
    	$remove_unique = $changes{'PG'}{'remove_unique'};
    }

    my $dict;
    if (exists $changes{'SQ'}{'from_dict'}) {
		$dict = $changes{'SQ'}{'from_dict'};
    }
    
	# Compare the SQ lines of the bam and the dict file
	if ($dict) {
		my $bfh = $stin->view($bam, undef, H => 1);
    	$bfh || $self->throw("Could not read header from '$bam'");
    	my @bheader;
    	while (<$bfh>) {
    		chomp;
    		next unless /\@SQ/;
    		push @bheader, $_;
    	}
    	close $bfh;
    	
        open my $dfh, "<$dict" || $self->throw("Could not open dictionary, $dict");
    	my @dheader;
    	while (<$dfh>) {
    		chomp;
    		next unless /\@SQ/;
    		push @dheader, $_;
    	}
    	close $bfh;
    	
    	my $dict_header = join "\n", @dheader;
    	my $bam_header = join "\n", @bheader;
    	
    	$dict = "" if ($dict_header eq $bam_header);
	}
    
    if (defined $changes{platform}) {
        $changes{platform} = $tech_to_platform{$changes{platform}} || $self->throw("Bad platform '$changes{platform}'");
    }
    
    my $new_header = $bam.'.new_header';
    $self->register_for_unlinking($new_header);
    
    # Parse original header and grab the lane id
    my $lane;
    if ($remove_unique) {
		my $pars = VertRes::Parser::sam->new(file => $bam);
		my %read_info = $pars->readgroup_info();
		my @lane_ids = keys %read_info;
		$self->throw("$bam does not have a single lane id.") unless (scalar @lane_ids == 1);
		$lane = $lane_ids[0];
    }
    
	open my $hfh, ">$new_header";
    $hfh || $self->throw("failed to get a filehandle for writing to '$new_header'");

	# Write new header
    my $made_changes = 0;
	my $bamfh = $stin->view($bam, undef, H => 1);
    $bamfh || $self->throw("Could not read header from '$bam'");
	my $sq_from_dict_done = 0;
    while (<$bamfh>) {
    	
		if (/^\@SQ/ && $dict) {
			$made_changes = 1;
        	unless ( $sq_from_dict_done ) {
        		open my $dictfh, "<$dict" || $self->throw("Could not open dictionary, $dict");
        		while (<$dictfh>) {
        			next unless /^\@SQ/;
					print $hfh $_;
        		}
		
        		close $dictfh;
        		$sq_from_dict_done = 1;
        	}
        	next;
        }

    	foreach my $line_tag (keys %changes) {
    		my $id_tag = $line_tag eq 'SQ' ? 'SN' : 'ID';
	        if (/^\@$line_tag.+$id_tag:([^\t]+)/) {
	            my $id = $1;
	            if (exists $changes{$line_tag}{$id}) {
	                while (my ($arg, $value) = each %{$changes{$line_tag}{$id}}) {
	                    next unless defined $value;
	                    my $tag = $arg_to_tag{$line_tag}{$arg} || next;
	                    
	                    if (/\t$tag:([^\t\n]+)/) {
	                        if ($1 ne $value) {
	                            $made_changes = 1;
	                            if ($value eq '') {
									s/\t$tag:[^\t\n]+//;
	                            } else {
	                            	s/\t$tag:[^\t\n]+/\t$tag:$value/;
	                            }
	                        }
	                    }
	                    else {
	                        $made_changes = 1;
	                        s/\n/\t$tag:$value\n/;
	                    }
	                }
	            }

	            if (exists $changes{$line_tag}{ all }) {
	                while (my ($arg, $value) = each %{$changes{$line_tag}{all}}) {
	                    next unless defined $value;
	                    my $tag = $arg_to_tag{$line_tag}{$arg} || next;
	                    
	                    if (/\t$tag:([^\t\n]+)/) {
	                        if ($1 ne $value) {
	                            $made_changes = 1;
	                            if ($value eq '') {
									s/\t$tag:[^\t\n]+//;
	                            } else {
	                            	s/\t$tag:[^\t\n]+/\t$tag:$value/;
	                            }
	                        }
	                    }
	                    else {
	                        $made_changes = 1;
	                        s/\n/\t$tag:$value\n/;
	                    }
	                }
	            }
	            
	            if ($remove_unique) {
					if (/\tCL:([^\t]+)/) {
			        	if ($1 =~ /\@|$lane/) {
							$made_changes = 1;
							s/[^\s\:]+=[^\s]+$lane[^\s]+//g;
							s/[^\s\:]+=[^\s]+\@[^\s]+//g;
							s/:\s+/:/;
							s/[^\S\n\t]+/ /g;
						}
			        }
				}
			}
        }
        print $hfh $_;
    }
    close($bamfh);
    close($hfh);
    
    # Exit with true if no changes were necessary
    unless ($made_changes) {
        unlink($new_header);
        return 1;
    }
    
    # Otherwise do the replacement and return its status
	my $status = $self->replace_bam_header($bam, $new_header);
	unlink($new_header);
	
	return $status;
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
                      study => 'DS',
                      project => 'DS');
    
    if (defined $rg_changes{platform}) {
        $rg_changes{platform} = $tech_to_platform{$rg_changes{platform}} || $self->throw("Bad platform '$rg_changes{platform}'");
    }
    
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
 Args    : path to bam and a hash with keys as readgroup identifiers and values
           as hash refs with at least one of the following:
           sample_name => string, eg. NA00000;
           library => string
           platform => string
           centre => string
           insert_size => int
           study => string, the study id, eg. SRP000001 (overwrites DS)

=cut

sub rewrite_bam_header {
    my ($self, $bam, %rg_changes) = @_;
    
    keys %rg_changes || return 1;
    my %arg_to_tag = (sample_name => 'SM',
                      library => 'LB',
                      platform => 'PL',
                      centre => 'CN',
                      insert_size => 'PI',
                      study => 'DS',
                      project => 'DS');
    
    if (defined $rg_changes{platform}) {
        $rg_changes{platform} = $tech_to_platform{$rg_changes{platform}} || $self->throw("Bad platform '$rg_changes{platform}'");
    }
    
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
    my $new_bam_lines = $self->num_bam_lines($temp_bam);
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

=head2 standardise_pg_header_lines

 Title   : standardise_pg_header_lines
 Usage   : $obj->standardise_pg_header_lines('a.bam',
                                'PG ID' => { VN => 'CL' });
 Function: For a given software and version, will change the CL to the one
           supplied, so that you can standardise PG lines that are supposed to
           be 'the same'.
 Returns : boolean (true on success, meaning a change was made and the bam has
           been confirmed OK, or no change was necessary and the bam was
           untouched)
 Args    : path to bam and a hash with keys as PG identifiers and values
           as hash refs where keys are version numbers and values are the
           CL command line you want to set, eg.
           'GATK TableRecalibration' => { '1.0.3016' => 'pQ=5, maxQ=40' }

=cut

sub standardise_pg_header_lines {
    my ($self, $bam, %pg_changes) = @_;
    
    keys %pg_changes || return 1;
    
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
        
        if (/^\@PG.+ID:([^\t]+)/) {
            my $pg = $1;
            if (exists $pg_changes{$pg}) {
                while (my ($vn, $cl) = each %{$pg_changes{$pg}}) {
                    next unless defined $cl;
                    
                    if (/\tVN:([^\t\n]+)/) {
                        if ($1 eq $vn) {
                            if (/\tCL:([^\t\n]+)/) {
                                if ($1 ne $cl) {
                                    $made_changes = 1;
                                    s/\tCL:[^\t\n]+/\tCL:$cl/;
                                }
                            }
                            else {
                                $made_changes = 1;
                                s/\n/\tCL:$cl\n/
                            }
                        }
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
    my $new_bam_lines = $self->num_bam_lines($temp_bam);
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
           study => string, the study id, eg. SRP000001
           md5 => string (NB: if you've just used rewrite_bam_header on a bam
                          and are making the same changes to its bas, the bam
                          md5 will have changed so you should supply the new
                          md5 here)
           filename => string OR ['/path/to/bam', 'release_date',
                                  '/path/to/sequence.index']
                                 (recalculates what the DCC filename of the bam
                                  should be)

=cut

sub rewrite_bas_meta {
    my ($self, $bas, %rg_changes) = @_;
    
    keys %rg_changes || return 1;
    my %arg_to_col = (sample_name => 3,
                      library => 5,
                      platform => 4,
                      study => 2,
                      project => 2,
                      md5 => 1,
                      filename => 0);
    # NAXXXXX.[chromN].technology.[center].algorithm.study_id.YYYYMMDD.bam
    my %arg_to_filename_regex = (sample_name => qr{^[^\.]+(\.)},
                                 platform => qr{(?:ABI_SOLID|ILLUMINA|LS454|SLX|454|SOLID)(\.)},
                                 project => qr{[^\.]+(\.\d{8})});
    
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
                if ($arg eq 'filename' && ref($value) && ref($value) eq 'ARRAY') {
                    $value = VertRes::Utils::Hierarchy->new->dcc_filename(@{$value});
                }
                
                my $col = $arg_to_col{$arg};
                if (defined $col) {
                    if ($cols[$col] ne $value) {
                        $made_changes = 1;
                        $cols[$col] = $value;
                    }
                }
                
                unless (exists $rg_changes{$rg}->{filename}) {
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

=head2 extract_intervals_from_bam

 Title   : extract_intervals_from_bam
 Usage   : $obj->extract_intervals_from_bam('in.bam', 'intervals_file', 'out.bam');
 Function: Takes a bam and makes a new bam file containing only those reads
           which overlap required intervals.  At least 1 base of a read needs to
           overlap an interval for it to be included.
           Input bam file must be sorted by coordinate.
 Returns : boolean (true on success)
 Args    : in.bam = input bam file
           intervals_file = file specifying the positions to keep.
             One line per interval, tab delimited: chromosome start end
             Order of intervals doesn't matter.  Any overlapping intervals
             will be merged into one interval, e.g.
             10 100 200
             10 120 170
             10 190 300
             would be merged into one interval 100-300 on chromosome 10.
           out.bam = name of output bam file (will be sorted by coordinate)
           opts = hash of optional options:
             max_read_length => int (default 200; length of longest read)

=cut

sub extract_intervals_from_bam {
    my ($self, $bam_in, $intervals_file, $bam_out, %opts) = @_;
    my $tmp_out = $bam_out . ".tmp";
    my %intervals;  # chromosome => [[start1,end1], [start2, end2], ...]
    my $lines_out_counter = 0;  # number of lines written to output file
    my $samtools = VertRes::Wrapper::samtools->new(verbose => $self->verbose, quiet => 1);
    my $bam_parser = VertRes::Parser::bam->new(file => $bam_in);
    my $result_holder = $bam_parser->result_holder();
    $bam_parser->flag_selector(self_unmapped => 0);
    $bam_parser->get_fields('CIGAR', 'QNAME', 'FLAG');
    %intervals = %{$self->_intervals_file2hash($intervals_file)};
    $opts{max_read_length} = 200 unless $opts{max_read_length};

    # we must have a bam index file
    $samtools->run_method('system');
    my $bai = $bam_in.'.bai';
    my $created_bai = 0;
    unless (-s $bai) {
        $samtools->index($bam_in, $bai);
        $samtools->run_status >= 1 || $self->throw("Failed to create $bai");
        $created_bai = 1;
    }
    
    # to keep the order of the output bam the same as the input, we need to
    # order the reference sequence names from intervals file to be
    # same as in the input bam header
    $samtools = VertRes::Wrapper::samtools->new(verbose => $self->verbose, quiet => 1, run_method => 'open');
    my $fh = $samtools->view($bam_in, undef, H => 1);
    my @ordered_ref_seqs;

    while (my $line = <$fh>) {
        if ($line =~ m/^\@SQ\tSN:(\w+)/){
            push @ordered_ref_seqs, $1 if $intervals{$1};
        }
    }

    close $fh;

    # extract intervals from input bam
    foreach my $chr (@ordered_ref_seqs) {
        my %reads_written;
        my $pos = 0;

        # make/append to the output bam file
        foreach my $interval (@{$intervals{$chr}}) {
            if ($interval->[0] - $pos > $opts{max_read_length}) {
                %reads_written = ();
            }

            $pos = $interval->[1];
            $bam_parser->region("$chr:$interval->[0]-$interval->[1]");
            while ($bam_parser->next_result) {
                next if $result_holder->{CIGAR} eq "*";

                unless ($reads_written{$result_holder->{QNAME} . $result_holder->{FLAG}}){
                    $bam_parser->write_result("$tmp_out");
                    $lines_out_counter++;
                    $reads_written{$result_holder->{QNAME} . $result_holder->{FLAG}} = 1;
                }
            }
        }
    }

    $bam_parser->close();
    unlink $bai if $created_bai;

    # check the right number of lines got written
    $lines_out_counter += $self->num_bam_header_lines("$bam_in");
    my $actual_lines = $self->num_bam_lines("$tmp_out");
    
    if ($actual_lines == $lines_out_counter) {
        rename "$tmp_out", $bam_out;
        return 1;
    }
    else {
        $self->warn("$tmp_out is bad, deleteing it (only $actual_lines lines vs $lines_out_counter)");
        unlink "$tmp_out";
        return 0;
    }
}

=head2 tag_strip

 Title   : tag_strip
 Usage   : $obj->tag_strip('in.bam', 'out.bam', qw(OQ XM XG XO));
 Function: Strips tags from bam records. By default all tags are stripped.
 Returns : boolean (true on success, meaning the output bam was successfully
           made)
 Args    : paths to input and output bams, list of tags to remove.

=cut

sub tag_strip {
    my ($self, $in_bam, $out_bam, @to_strip) = @_;
    @to_strip > 0 || $self->throw("You must supply tags to strip");
    
    my $pb = VertRes::Parser::bam->new(file => $in_bam);
    $pb->ignore_tags_on_write(@to_strip);
    my $tmp_out = $out_bam.'.running';
    
    my $lines = 0;
    while ($pb->next_result) {
        $lines++;
        $pb->write_result($tmp_out);
    }
    $pb->close;
    
    # check for truncation
    my $actual_lines = $self->num_bam_records($tmp_out);
    
    if ($actual_lines == $lines) {
        move($tmp_out, $out_bam) || $self->throw("Failed to move $tmp_out to $out_bam");
        return 1;
    }
    else {
        $self->warn("$tmp_out is bad (only $actual_lines lines vs $lines)");
        unlink $tmp_out;
        return 0;
    }
}

=head2 bam2fastq

 Title   : bam2fastq
 Usage   : $obj->bam2fastq('in.bam', 'fastq_base');
 Function: Converts bam to fastq and runs fastqcheck on the result. 
           bamcheck will be run, if a bamcheck file does not already exist.
 Returns : boolean (true on success, meaning the output bam was successfully
           made)
 Args    : paths to input and output bams, list of tags to remove.

=cut

sub bam2fastq {
    my ($self, $bam, $fastq_base) = @_;

	# read bam info from bamcheck file - create if necessary
	my $bamcheck = $bam.'.bc';
	unless (-s $bamcheck) {
        Utils::CMD(qq[bamcheck $bam > $bamcheck.tmp]);
        move("$bamcheck.tmp", $bamcheck) || $self->throw("Could not rename '$bamcheck.tmp' to '$bamcheck'\n");
	}
	my $bc_parser = VertRes::Parser::bamcheck->new(file => $bamcheck);
    my $total_reads = $bc_parser->get('sequences');
    my $total_bases = $bc_parser->get('total_length');
	my $is_paired = $bc_parser->get('is_paired');
	$self->throw("failed to parse $bamcheck\n") unless ($total_reads && $total_bases);
	
    my $fsu = VertRes::Utils::FileSystem->new();
    my (undef, $path) = fileparse($bam);
	my @out_fastq;
	if ($is_paired) {
		push @out_fastq, $fsu->catfile($path, "${fastq_base}_1.fastq");
		push @out_fastq, $fsu->catfile($path, "${fastq_base}_2.fastq");
	} else {
		push @out_fastq, "${fastq_base}.fastq";
	}
	
	my @tmp_fastq = map("$_.tmp.fastq", @out_fastq);
    
	# picard needs a tmp dir, but we don't use /tmp because it's likely to fill up
    my $tmp_dir = $fsu->tempdir('_bam2fastq_tmp_XXXXXX', DIR => $path);
    
    my $verbose = $self->verbose();
    my $picard = VertRes::Wrapper::picard->new(verbose => $verbose,
                                                quiet => $verbose ? 0 : 1,
                                                $self->{java_memory} ? (java_memory => $self->{java_memory}) : (),
                                                validation_stringency => 'silent',
                                                tmp_dir => $tmp_dir);
    
    $picard->SamToFastq($bam, @tmp_fastq);
	$self->throw("SamToFastq of $bam failed\n") unless $picard->run_status >= 1;
	
	# run fastqcheck on all of the fastq files made record total length and number of reads
	my $num_bases = 0;
	my $num_sequences = 0;
	foreach my $precheck_fq (@tmp_fastq) {
		my $fqc_file = $precheck_fq;
		$fqc_file =~ s/\.tmp\.fastq$/.fastqcheck/;
		unless (-s $fqc_file) {
		    # make the fqc file
		    my $fqc = VertRes::Wrapper::fastqcheck->new();
		    $fqc->run($precheck_fq, $fqc_file.'.temp');
		    $self->throw("fastqcheck on $precheck_fq failed - try again?") unless $fqc->run_status >= 1;
		    $self->throw("fastqcheck failed to make the file $fqc_file.temp") unless -s $fqc_file.'.temp';
	    
		    # check it is parsable
		    my $parser = VertRes::Parser::fastqcheck->new(file => $fqc_file.'.temp');
		    my $num_seq = $parser->num_sequences();
		    if ($num_seq) {
		        move($fqc_file.'.temp', $fqc_file) || $self->throw("failed to rename '$fqc_file.temp' to '$fqc_file'");
		    }
		    else {
				unlink $fqc_file.'.temp';
		        $self->throw("fastqcheck file '$fqc_file.temp' was created, but doesn't seem valid\n");
		    }
		}
	    
	    if (-s $fqc_file) {
		    my $parser = VertRes::Parser::fastqcheck->new(file => $fqc_file);
		    $num_sequences += $parser->num_sequences();
		    $num_bases += $parser->total_length();		    
		}
	}
	
    # check fastqcheck output against expectation
    my $ok = 0;
	if ($num_sequences == $total_reads) {
	    $ok++;
	} else {
	    $self->warn("fastqcheck reports number of reads as $num_sequences, not the expected $total_reads\n");
	}
	if ($num_bases == $total_bases) {
	    $ok++;
	} else {
	    $self->warn("fastqcheck reports number of bases as $num_bases, not the expected $total_bases\n");
	}
	
	# if everything checks out rename the tmp fastqs, otherwise unlink created files
    if ($ok == 2) {
		foreach my $tmp_fq (@tmp_fastq) {
			my $fq = $tmp_fq;
			$fq =~ s/\.tmp\.fastq$//;
	        move($tmp_fq, $fq) || $self->throw("Could not rename '$tmp_fq' to '$fq'\n");
		}
    } else {
		foreach my $tmp_fq (@tmp_fastq) {
			unlink $tmp_fq;
		}
        $self->throw("fastqcheck of fastq file doesn't match expectations\n");
    }
}

=head2 filter_readgroups

 Title   : filter_readgroups
 Usage   : $obj->filter_readgroups('in.bam','out.bam',include=>[{SM=>'HG01522',PL=>'ILLUMINA'}]);
 Function: Given a bam file, generates another bam file that contains only the
           requested read groups. Currently only 'include' logic implemented; it is trivial to add 'exclude' though.
           Note: The header is not modified.
 Returns : 1 on success or 0 when the filters would not exclude any reads
 Args    : starting bam file, output name for bam file and filter specification

=cut

sub filter_readgroups
{
    my ($self, $in_bam, $out_bam, %opts) = @_;

    if ( !exists($opts{include}) ) { $self->throw("Missing the 'include' parameter.\n"); }
    if ( ! -e $in_bam ) { $self->throw("No such file: $in_bam\n"); }

    my $pars = VertRes::Parser::bam->new(file=>$in_bam);
    my $tmp  = $out_bam.'.part';

    # Init the read groups: loop through all RG records of the BAM file and see if they pass
    #   any of the filters.
    #
    my %include;
    my %rg_info = $pars->readgroup_info();
    while (my ($rg,$info) = each %rg_info)
    {
        # Will the RG pass any of these filters?
        for my $inc (@{$opts{include}})
        {
            # All requested tags must match
            my $match = 0;
            while (my ($tag,$value) = each %$inc)
            {
                if ( exists($$info{$tag}) && $$info{$tag} eq $value ) { $match++; }
            }
            if ( $match == scalar keys %$inc ) 
            {
                # All tags match the criteria, include this RG
                $include{$rg} = 1;
                last;
            }
        }
    }

    if ( !scalar keys %include ) { $self->throw("The filter is too strict, no read group passes the criteria.\n"); }
    if ( scalar keys %include == scalar keys %rg_info ) { $self->warn("The filter matches all read groups, use cp instead..."); return 0; }

    $pars->get_fields('RG');
    my $rh = $pars->result_holder();
    
    while ($pars->next_result()) 
    {
        if (exists $include{$rh->{RG}}) 
        {
            $pars->write_result($tmp);
        }
    }

    rename($tmp,$out_bam) or $self->throw("rename $tmp $out_bam: $!");
}

=head2 bam_exome_qc_stats

 Title   : bam_exome_qc_stats
 Usage   : $obj->bam_exome_qc_stats('in.bam', %options);
 Function: calculates QC statistics from a BAM file of exome reads.
 Returns : reference to hash of statistics, contents of which are: 
               readlen                  = length of first read (we assume all reads same length)
               bait_bases               = number of bases in the baits (after merging)
               bait_bases_mapped        = number of bases mapped onto a bait
               bait_near_bases          = number of bases near to, but not in, a bait
               bait_near_bases_mapped   = number of bases mapped near to, but not on, a bait
               baits                    = total number of baits (after merging)
               bait_design_efficiency   = target_bases / bait_bases.  (1 => perfectly efficient,
                                          0.5 => half of baited bases not on target)
               mean_bait_coverage       = bait_bases_mapped / bait_bases
               mean_target_coverage     = target_bases_mapped / target_bases
               bait_coverage_sd         = standard deviation of coverage distribution of bait bases
               target_coverage_sd       = standard deviation of coverage distribution of target bases
               off_bait_bases           = number of bases not mapped on or near to a bait
               bases_mapped             = number of bases mapped (to anything)
               bases_mapped_reads       = total length of mapped reads
               clip_bases               = number of clipped bases (S or H in cigar string)
               mapped_as_pair           = number of reads mapped as a proper pair (0x0002 in flag)
               num_mismatches           = number of mismatches (total of NM:i:.. values)
               pct_mismatches           = 100 * num_mismatches / bases_mapped_reads
               reads_paired             = number of paired reads (0x0001 in flag)
               raw_reads                = number of records in the input bam file
               reads_mapped             = number of reads mapped (0x0004 not in flag)
               rmdup_reads_mapped       = number of reads mapped excluding duplicates
                                          (i.e. don't count reads with 0x0400 in flag)
               rmdup_bases_mapped       = as previous, but bases
               target_bases             = total number of target bases (after merging)
               target_bases_mapped      = number of bases mapped on a target
               target_near_bases        = total number of bases near to, but not on, a target
               target_near_bases_mapped = total number of bases mapped near to, but not on a target
               targets                  = total number of targets (after merging)
               raw_bases                = total number of bases of all reads in bam file
               reads_on_target          = number of reads mapped to a target (at least 1 base
                                          mapped onto a target)
               reads_on_target_near     = as above, but mapped near to, not on, a target
               reads_on_bait            = as reads_on_target, but for baits
               reads_on_bait_near       = as above, but mapped near to, not on, a bait
               low_cvg_targets          = % targets that didn't reach a minimum coverage over
                                          any base.  Minimum coverage is set by low_cvg in
                                          options hash.
               mean_insert_size         = mean insert size of all read pairs
               median_insert_size       = median insert size of all read pairs
               insert_size_sd           = standard deviation of insert size

               gc_hist_unmapped_1 = array of length readlen + 1, i'th element is the number
                                    of unmapped reads with i bases equal to G or C.
               gc_hist_unmapped_2 = same as gc_hist_unmapped_1, but for 2nd of pair
               gc_hist_mapped_1   = same as gc_hist_unmapped, but for mapped reads, 1st of pair
               gc_hist_mapped_2   = same as gc_hist_unmapped, but for mapped reads, 2nd of pair
  
               qual_scores_1      = an array of hashes.  i'th hash stores quality values of the
                                    (i+1)'th base of each read which is first of a pair.
                                    hash: quality => count.
                                    e.g. qual_box_data[3]{30} = 42 means that the 4th base of
                                    42 reads had quality score 30.
                                    Empty if no reads in the input bam have flag 0x0040.
               qual_scores_2      = same as qual_box_data_1, but reads which are second of a pair,
                                    i.e. those with 0x0080 in the flag.
               qual_scores_up     = same as for qual_box_data_1, but reads which are not paired,
                                    i.e. those without 0x0001 in the flag.
               bait_cvg_hist      = hash of coverage => number of bait bases with that coverage
               target_cvg_hist    = hash of coverage => number of target bases with that coverage
               gc_vs_bait_cvg     = hash to see how coverage varies with GC content.
                                    %GC (nearest int) => {mean coverage (nearest int) => count}
               gc_vs_target_cvg   = as for gc_vs_bait_cvg, but for targets instead of baits
               insert_size_hist   = hash of insert size -> count of read pairs
               bait_cumulative_coverage = hash of coverage => number of bait bases achieving
               target_cumulative_coverage = as previous, but targets instead of baits
               target_bases_1X   = fraction of target bases with >= 1X coverage
               target_bases_2X   = fraction of target bases with >= 2X coverage
               target_bases_5X   = fraction of target bases with >= 5X coverage
               target_bases_10X  = fraction of target bases with >= 10X coverage
               target_bases_20X  = fraction of target bases with >= 20X coverage
               target_bases_50X  = fraction of target bases with >= 50X coverage
               target_bases_100X = fraction of target bases with >= 100X coverage

 Args    : Options in a hash:
             MANDATORY: either all four of:
                          - bait_interval, target_interval, ref_fa, ref_fai
                        or this:
                          - load_intervals_dump_file
                        must be given in the options hash.   Info about baits/targets and reference
                        is calulcated prior to parsing the BAM file.  This info can be dumped
                        to a file to save calculating it more than once, by running first with
                        -dump_intervals, then do the QC with -load_intervals_dump_file.

             bait_interval   => filename of baits intervals file, of the form 1 bait per line, tab
                                separated:
                                reference_seq start end
                                Order doesn't matter, and baits can overlap.  This will be dealt with.
             target_interval => filename of targets intervals file, same format as bait_interval.
             ref_fa          => filename of reference fasta file
             ref_fai         => fai file of reference fasta file
             dump_intervals  => name of file to dump info to.  If this is used, the BAM is NOT QC'd.
             load_intervals_dump_file => name of file made when dump_intervals was used.  This will
                                         load all reference/bait/target info from the dump file, which
                                         is faster than creating it from scratch
             bam             => name of bam file to be QC'd.
                                MANDATORY, unless the dump_intervals option is used (in which case
                                the BAM file is ignored)
             low_cvg         => int (default 2).  Used to make the statistic low_cvg_targets.
                                Any target with all its bases' coverage less than this value will
                                be counted as low coverage.  e.g. default means that a target
                                must have coverage of at least 2 at one base to not count as a
                                low coverage target.
             near_length     => int (default 100)
                                Distance to count as 'near' to a bait or a target
             max_insert      => int (default 2000) Ignore insert sizes > this number.  Sometimes
                                (depending on the mapper) BAM files can have unrealisticly large
                                insert sizes.  This option removes the outliers, so calculated
                                mean insert size is more realistic.

=cut

sub bam_exome_qc_stats {
    my ($self, %opts) = @_;
    my $intervals;
    my $ref_lengths; 
    my $current_rname;
    $opts{near_length} = 100 unless $opts{near_length};
    $opts{low_cvg} = 2 unless $opts{low_cvg};
    $opts{max_insert} = 1000 unless $opts{max_insert};
    my $first_sam_record = 1;

    if ($opts{fake_it}) {  # top-secret(!) option for quick pipeline debugging
        my $s = do $opts{fake_it} or $self->throw("Error loading data from file $opts{fake_it}");
        print STDERR "VertRes::Utils::Sam->bam_exome_qc_stats using fake file $opts{fake_it}\n";
        return $s;
    }

    if (!($opts{dump_intervals}) and !($opts{bam})) {
        $self->throw("If dump_intervals not used, then a BAM file must be given");
    }


    my %stats = ('bait_bases', 0,
                 'bait_bases_mapped', 0,
                 'bait_near_bases_mapped', 0,
                 'bait_near_bases', 0,
                 'bases_mapped', 0,
                 'bases_mapped_reads', 0,
                 'clip_bases', 0,
                 'low_cvg_targets', 0,
                 'mapped_as_pair', 0,
                 'num_mismatches', 0,
                 'off_bait_bases', 0,
                 'raw_bases', 0,
                 'raw_reads', 0,
                 'reads_mapped', 0,
                 'reads_on_target', 0,
                 'reads_on_target_near', 0,
                 'reads_on_bait', 0,
                 'reads_on_bait_near', 0,
                 'reads_paired', 0,
                 'rmdup_reads_mapped', 0,
                 'rmdup_bases_mapped', 0,
                 'target_bases', 0,
                 'target_bases_mapped', 0,
                 'target_near_bases', 0,
                 'target_near_bases_mapped', 0);

    my %pileup; # for the current reference sequence in the bam file ($current_rname), this
                # will store pileup info of the target bases, near target bases, bait bases and
                # near bait bases.  One hash for each.  reference position => number of reads
                # covering that position.  'near' does *not* include on, i.e. a given position
                # cannot be both near to a target and in a target (similarly for baits).

    # Need to get the target/bait info.  Either it was precomputed and dumped
    # to a file we need to load, or we need to calculate it all now.
    # (Calulcating the GC of each interval is slow)
    if ($opts{load_intervals_dump_file}) {
        $intervals = do $opts{load_intervals_dump_file} or $self->throw("Error loading data from load_intervals_dump_file file $opts{load_intervals_dump_file}");
    }
    # parse the interval files and get gc from reference fasta
    elsif ($opts{bait_interval} and $opts{target_interval} and $opts{ref_fa} and $opts{ref_fai}) {
        $ref_lengths = Utils::fai_chromosome_lengths($opts{ref_fai});
        $intervals->{target} = $self->_intervals_file2hash($opts{target_interval});
        $intervals->{target_near} = $self->_intervals_file2hash($opts{target_interval}, ('expand_intervals' => $opts{near_length}));
        $intervals->{bait} = $self->_intervals_file2hash($opts{bait_interval});
        $intervals->{bait_near} = $self->_intervals_file2hash($opts{bait_interval}, ('expand_intervals' => $opts{near_length}));
        

        foreach my $bait_or_target qw(bait target) {
            my $pars = VertRes::Parser::fasta->new(file => $opts{ref_fa});
            my $result_holder = $pars->result_holder();

            while ($pars->next_result(1)) {
                my $id = \$result_holder->[0];
                my $seq = \$result_holder->[1];
                exists $intervals->{$bait_or_target}{$$id} || next;

                foreach my $interval (@{$intervals->{$bait_or_target}{$$id}}) {
                    my $length = $interval->[1] - $interval->[0] + 1;
                    # coords of fasta sequence are zero-based, interval coords are 1-based
                    my $gc_count = substr($$seq, $interval->[0] - 1, $length) =~ tr/cgCG//;
                    push @{$interval}, int(100 * $gc_count / $length);
                } 
            }
        }

        if ($opts{dump_intervals}) {
            open my $fh, '>', $opts{dump_intervals} or self->throw("Error opening file $opts{dump_intervals}");
            print $fh Dumper($intervals);
            close $fh;
            return %stats;
        }
    }
    else {
        $self->throw("Error. must use either intervals_dump_file, or give all four of ref_fa, ref_fai, bait_interval and target_interval\n");
    }

    print STDERR "Intervals etc sorted. Parse BAM file...\n" if $opts{verbose};

    # Everything is initialised.  It's time to parse the bam file to get the stats
    my $bam_parser = VertRes::Parser::bam->new(file => $opts{bam});
    my $result_holder = $bam_parser->result_holder();
    $bam_parser->get_fields('FLAG', 'RNAME', 'POS', 'CIGAR', 'ISIZE', 'SEQ', 'QUAL', 'SEQ_LENGTH', 'MAPPED_SEQ_LENGTH', 'NM');

    while ($bam_parser->next_result()) {
        my $flag = $result_holder->{FLAG};
        my $rname = $result_holder->{RNAME};
        my $pos = $result_holder->{POS};
        my $cigar = $result_holder->{CIGAR};
        my $isize = $result_holder->{ISIZE};
        my $seq = $result_holder->{SEQ};
        my $qual = $result_holder->{QUAL};
        my $seq_length = $result_holder->{SEQ_LENGTH};
        my $mapped_seq_length = $result_holder->{MAPPED_SEQ_LENGTH};
        my $num_mismatches = $result_holder->{NM};
        my $gc_array = $stats{gc_hist_unmapped_1};
        my $qual_array = $stats{qual_scores_up};
        
        # if first line of SAM, initlalize hashes/lists etc
        if ($first_sam_record) {
            $stats{readlen} = length $qual;
            $stats{paired} = $bam_parser->is_sequencing_paired($flag);
            $stats{gc_hist_unmapped_1} = [(0) x ($stats{readlen} + 1)];
            $stats{gc_hist_unmapped_2} = [(0) x ($stats{readlen} + 1)];
            $stats{gc_hist_mapped_1} = [(0) x ($stats{readlen} + 1)];
            $stats{gc_hist_mapped_2} = [(0) x ($stats{readlen} + 1)];
            
            foreach my $a qw(gc_vs_bait_cvg gc_vs_target_cvg) {
                $stats{$a} = ();
                foreach (0 .. 100) {
                    push @{$stats{$a}}, {};
                }
            }

            foreach my $a qw(qual_scores_1 qual_scores_2 qual_scores_up) {
                $stats{$a} = ();
                foreach (1 .. $stats{readlen}) {
                    push @{$stats{$a}}, {};
                }
            }

            $first_sam_record = 0;
        }

        $stats{raw_reads}++; 
        $stats{raw_bases} += length $qual;
        $stats{reads_paired}++ if $bam_parser->is_sequencing_paired($flag);
        $stats{num_mismatches} += $num_mismatches unless $num_mismatches eq '*';
        
        if ($bam_parser->is_mapped($flag)) {
            $stats{reads_mapped}++;
            $stats{bases_mapped} += $mapped_seq_length;
            $stats{mapped_as_pair}++ if $bam_parser->is_mapped_paired($flag);
            $stats{insert_size_hist}{$isize}++ if ($bam_parser->is_sequencing_paired($flag) and $isize > 0 and $isize < $opts{max_insert});
            $stats{bases_mapped_reads}+= $seq_length;
            
            unless ($bam_parser->is_duplicate($flag)) {
                $stats{rmdup_reads_mapped}++;
                $stats{rmdup_bases_mapped} += $mapped_seq_length;
            }
            
            if ($bam_parser->is_first($flag)) {
                $gc_array = $stats{gc_hist_mapped_1};
            }
            else {
                $gc_array = $stats{gc_hist_mapped_2}; 
            }
        }
        elsif ($bam_parser->is_second($flag)) {
            $gc_array = $stats{gc_hist_unmapped_2};
        }
            
        $gc_array->[$seq =~ tr/cgCG/cgCG/]++;

        # update quality scores arrays
        $qual_array = $stats{qual_scores_1} if $bam_parser->is_first($flag);
        $qual_array = $stats{qual_scores_2} if $bam_parser->is_second($flag);
        my @qual_ints = VertRes::Utils::FastQ->qual_to_ints($qual); 
        @qual_ints = reverse @qual_ints if $bam_parser->is_reverse_strand($flag);
        foreach my $i (0 .. $#qual_ints){
            $qual_array->[$i]{$qual_ints[$i]}++;
        }

        print STDERR "$stats{raw_reads} reads processed\n" if ($opts{verbose} && $stats{raw_reads} % 1000000 == 0);

        # nothing left to do if the read is unmapped...
        next unless $bam_parser->is_mapped($flag);

        # if current reference sequence doesn't have any baits or targets...
        unless (exists $intervals->{target}{$rname} or exists $intervals->{bait}{$rname}) {
            $stats{off_bait_bases} += $mapped_seq_length; 
            $stats{off_target_bases} += $mapped_seq_length;
            next;
        }

        # have we just reached a new reference sequence?
        if ((!defined $current_rname) or $current_rname ne $rname) {
            # sort out stats etc for the ref sequence we just finished
            if (defined $current_rname) {
                $self->_update_bam_exome_qc_stats($intervals, \%pileup, \%stats, $current_rname, \%opts);
                print STDERR "Chromosome $current_rname stats done\n" if $opts{verbose};
            }

            $current_rname = $rname;

            # initialize the pileup hashes for the new reference sequence, which
            # means setting the value of each position to be zero
            foreach (keys %{$intervals}) {
                $pileup{$_} = {};
            }  
            foreach my $interval_type ('bait', 'target') {
                foreach my $interval_array ($intervals->{$interval_type}{$current_rname}) {
                    foreach my $interval (@{$interval_array}) {
                        foreach my $i ($interval->[0] .. $interval->[1]) {
                            $pileup{$interval_type}{$i} = 0;
                        }
                    }
                }
            }

            # That's baits/targets initialized.  Now need to do the near bases.
            # We don't want a 'near base' to be actually in a target/bait, since doing so
            # would mean we're recording the same information twice, using extra memory
            foreach my $interval_type ('bait', 'target') {
                foreach my $interval_array ($intervals->{$interval_type . '_near'}{$current_rname}) {
                    foreach my $interval (@{$interval_array}) {
                        foreach my $i ($interval->[0] .. $interval->[1]) {
                            $pileup{$interval_type . '_near'}{$i} = 0 unless exists $pileup{$interval_type}{$i};
                        }
                    }
                }
            }
        }

        # only gather bait/target stats for reads which are not duplicates
        next if $bam_parser->is_duplicate($flag);

        # update the number of reads hitting targets/baits and the target/bait coverage
        my $ref_pos = $pos;  
        my @numbers = split /[A-Z]/, $cigar;
        my @letters = split /[0-9]+/, $cigar;
        shift @letters;
        my %read_hits;
        
        foreach my $letter_index (0 ..  $#letters) {
            if ($letters[$letter_index] eq 'M') {

                foreach my $i ($ref_pos .. $ref_pos + $numbers[$letter_index] - 1) {
                    foreach my $key ('bait', 'target') {
                        if (exists $pileup{$key}{$i}) {
                            $pileup{$key}{$i}++;
                            $read_hits{$key} = 1;
                            $read_hits{$key . '_near'} = 1;
                        }
                        elsif (exists $pileup{$key . '_near'}{$i}) {
                            $pileup{$key . '_near'}{$i}++;
                            $read_hits{$key . '_near'} = 1;
                        }
                    }

                    $stats{off_bait_bases}++ unless (exists $pileup{bait}{$i} or exists $pileup{bait_near}{$i});
                    $stats{off_target_bases}++ unless (exists $pileup{target}{$i} or exists $pileup{target_near}{$i});
                } 

                $ref_pos += $numbers[$letter_index];
            }
            elsif ($letters[$letter_index] eq 'D' or $letters[$letter_index] eq 'N') {
                $ref_pos += $numbers[$letter_index];
            }
            elsif ($letters[$letter_index] eq 'S' or $letters[$letter_index] eq 'H') {
                $stats{clip_bases} += $numbers[$letter_index];
            }
        }

        delete $read_hits{bait_near} if $read_hits{bait};
        delete $read_hits{target_near} if $read_hits{target};

        foreach (keys %read_hits) {
            $stats{'reads_on_' . $_}++;
        }
    } 

    $bam_parser->close();
    print STDERR "Finished parsing BAM file\n" if $opts{verbose};

    $stats{'bait_gc_hist'} = [(0) x 101 ];
    $stats{'target_gc_hist'} = [(0) x 101];

    # add bait/target gc histogram to qc dump hash
    for my $b_or_t ('bait', 'target'){
        foreach my $chr (keys %{$intervals->{$b_or_t}}) {
            foreach my $a (@{$intervals->{$b_or_t}{$chr}}) {
                $stats{$b_or_t . '_gc_hist'}[$a->[2]]++;
            }
        }
    }

    # sort out stats for the last ref sequence in the bam file
    $self->_update_bam_exome_qc_stats($intervals, \%pileup, \%stats, $current_rname, \%opts);
    undef %pileup;
    undef $intervals;

    # sort out other global stats
    $stats{bait_design_efficiency} = $stats{bait_bases} != 0 ? $stats{target_bases} / $stats{bait_bases} : 'NA';
    $stats{pct_mismatches} = 100 * $stats{num_mismatches} / $stats{bases_mapped_reads};
    my %qual_stats = VertRes::Utils::Math->new()->histogram_stats($stats{insert_size_hist});
    $stats{median_insert_size} = $qual_stats{q2};
    $stats{mean_insert_size} = $qual_stats{mean};
    $stats{insert_size_sd} = $qual_stats{standard_deviation};
    $stats{mean_bait_coverage} = $stats{bait_bases_mapped} / $stats{bait_bases};
    $stats{mean_target_coverage} = $stats{target_bases_mapped} / $stats{target_bases};

    # calculate cumulative coverage of baits and targets
    foreach my $bait_or_target qw(bait target) {
        my %cov_stats = VertRes::Utils::Math->new()->histogram_stats($stats{$bait_or_target . '_cvg_hist'});
        $stats{$bait_or_target . '_coverage_sd'} = $cov_stats{standard_deviation};

        my $base_count = 0;
        my $total_bases = $stats{$bait_or_target . '_bases'};

        foreach my $cov (reverse sort {$a <=> $b} keys %{$stats{$bait_or_target . '_cvg_hist'}}) {
            $base_count += $stats{$bait_or_target . '_cvg_hist'}{$cov};
            $stats{$bait_or_target . '_cumulative_coverage'}{$cov} = $base_count;
            $stats{$bait_or_target . '_cumulative_coverage_pct'}{$cov} = $base_count / $stats{$bait_or_target . '_bases'};
        }
    }

    # work out % bases coverage above certain values
    foreach my $bait_or_target qw(bait target) { 
        my @vals = (1, 2, 5, 10, 20, 50, 100);
        foreach (@vals) {
            $stats{$bait_or_target . '_bases_' . $_ . 'X'} = 0;
        }

        foreach my $cov (sort {$a <=> $b} keys %{$stats{$bait_or_target . '_cumulative_coverage_pct'}}) {
            while ((0 < scalar @vals) and $cov >= $vals[0]) {
                $stats{$bait_or_target . '_bases_' . $cov . 'X'} = $stats{$bait_or_target . '_cumulative_coverage_pct'}{$cov};
                shift @vals;
            }
            last if (0 == scalar @vals);
        }
    }

    return \%stats;
}

=head2 bam_exome_qc_make_plots

 Title   : bam_exome_qc_make_plots
 Usage   : $obj->bam_exome_qc_make_plots(\%stats, 'outfiles_prefix', 'plot_type');
 Function: Makes exome QC plots using data gathered by bam_exome_qc_stats
 Returns : nothing (makes plots)
 Args    : reference to hash (returned by bam_exome_qc_stats)
           prefix of output plot files to be made
           plot type (must be recognised by R, e.g. 'pdf', 'png', ...)

=cut

sub bam_exome_qc_make_plots {
    my ($self, $stats, $outfiles_prefix, $plot_type) = @_;

    # Plot: mean coverage distribution of baits/targets
    # This info is in the gc vs coverage info, so extract it then plot
    my %coverage;
    my %plot_data;

    foreach my $bait_or_target qw(bait target) {
        foreach my $hash (@{$stats->{'gc_vs_' . $bait_or_target . '_cvg'}}) {
            while (my ($cov, $freq) = each(%{$hash})) {
                $coverage{$bait_or_target}{$cov} += $freq;
            }
        }
        my %cov_stats = VertRes::Utils::Math->new()->histogram_stats($coverage{$bait_or_target});
        
        # plot 4 standard devations away from the mean
        my $xmin = max (0, $cov_stats{mean} - 4 * $cov_stats{standard_deviation});
        my $xmax = $cov_stats{mean} + 4 * $cov_stats{standard_deviation};

        $plot_data{$bait_or_target}{x} = [];
        $plot_data{$bait_or_target}{y} = [];

        foreach my $i ($xmin .. $xmax) {
            if (exists $coverage{$bait_or_target}{$i}) {
                push @{$plot_data{$bait_or_target}{x}}, $i;
                push @{$plot_data{$bait_or_target}{y}}, $coverage{$bait_or_target}{$i};
            }
        }
    }
    my @data = ({xvals => $plot_data{bait}{x}, yvals => $plot_data{bait}{y}, legend => 'baits'},
                {xvals => $plot_data{target}{x}, yvals => $plot_data{target}{y}, legend => 'targets'},);

    Graphs::plot_stats({outfile=>"$outfiles_prefix.mean_coverage.$plot_type",
                        title => "Mean coverage distribution of baits and targets",
                        desc_xvals => 'Mean coverage',
                        desc_yvals => 'Frequency',
                        data => \@data});

    # Coverage plots
    %plot_data = ();

    foreach my $bait_or_target qw(bait target) {
        my %cov_stats = VertRes::Utils::Math->new()->histogram_stats($stats->{$bait_or_target . '_cvg_hist'});
        my $xmin = max (0, $cov_stats{mean} - 4 * $cov_stats{standard_deviation});
        my $xmax = $cov_stats{mean} + 4 * $cov_stats{standard_deviation};

        $plot_data{$bait_or_target}{x_per_base} = [];
        $plot_data{$bait_or_target}{y_per_base} = [];
        $plot_data{$bait_or_target}{x_norm} = [];
        $plot_data{$bait_or_target}{y_norm} = [];
        $plot_data{$bait_or_target}{x_cumulative} = [];
        $plot_data{$bait_or_target}{y_cumulative} = [];

        # get the per-base coverage data
        foreach my $i ($xmin .. $xmax) {
            if (exists $stats->{$bait_or_target . '_cvg_hist'}{$i}) {
                push @{$plot_data{$bait_or_target}{x_per_base}}, $i;
                push @{$plot_data{$bait_or_target}{y_per_base}}, $stats->{$bait_or_target . '_cvg_hist'}{$i};
            }
        }

        foreach my $cov (sort {$a <=> $b} keys %{$stats->{$bait_or_target . '_cumulative_coverage_pct'}}) {
            push @{$plot_data{$bait_or_target}{x_cumulative}}, $cov;
            push @{$plot_data{$bait_or_target}{y_cumulative}}, $stats->{$bait_or_target . '_cumulative_coverage_pct'}{$cov};
            last if $cov >= 2 * $stats->{mean_target_coverage};
        }

        # work out the normalised coverage plot data
        foreach my $cov (sort {$a <=> $b} keys %{$stats->{$bait_or_target . '_cumulative_coverage'}}) {
            my $bases_fraction = $stats->{$bait_or_target . '_cumulative_coverage'}{$cov} / $stats->{$bait_or_target . '_bases'};
            my $normalised_cov = $cov / $stats->{'mean_' . $bait_or_target . '_coverage'};
            push @{$plot_data{$bait_or_target}{x_norm}}, $normalised_cov;
            push @{$plot_data{$bait_or_target}{y_norm}}, $bases_fraction;
            last if $normalised_cov > 1;
        }
    }

    @data = ({xvals => $plot_data{bait}{x_per_base}, yvals => $plot_data{bait}{y_per_base}, legend => 'baits'},
                {xvals => $plot_data{target}{x_per_base}, yvals => $plot_data{target}{y_per_base}, legend => 'targets'},);

    Graphs::plot_stats({outfile=>"$outfiles_prefix.coverage_per_base.$plot_type",
                        title => "Bait and target coverage per base",
                        desc_xvals => 'Coverage',
                        desc_yvals => 'Frequency',
                        data => \@data});

    @data = ({xvals => $plot_data{bait}{x_norm}, yvals => $plot_data{bait}{y_norm}, legend => 'baits'},
                {xvals => $plot_data{target}{x_norm}, yvals => $plot_data{target}{y_norm}, legend => 'targets'},);

    Graphs::plot_stats({outfile=>"$outfiles_prefix.normalised_coverage.$plot_type",
                        title => "Bait and target normalised coverage",
                        desc_xvals => 'Normalised coverage',
                        desc_yvals => 'Fraction of bases',
                        r_plot => 'ylim=c(0,1)', 
                        data => \@data});

    @data = ({xvals => $plot_data{bait}{x_cumulative}, yvals => $plot_data{bait}{y_cumulative}, legend => 'baits'},
                {xvals => $plot_data{target}{x_cumulative}, yvals => $plot_data{target}{y_cumulative}, legend => 'targets'},);

    Graphs::plot_stats({outfile=>"$outfiles_prefix.cumulative_coverage.$plot_type",
                        title => "Bait and target cumulative coverage",
                        desc_xvals => 'Coverage',
                        desc_yvals => 'Fraction of bases',
                        r_plot => 'ylim=c(0,1)', 
                        data => \@data});

    # Plot: GC of baits/targets vs mean, quartiles coverage
    foreach (qw(bait target)) {
        Graphs::plot_histograms_distributions({outfile=>"$outfiles_prefix.$_" . "_gc_vs_cvg.$plot_type",
                                               title => "Coverage vs GC plot of $_" . 's',
                                               desc_xvals => 'GC (%)',
                                               desc_yvals => 'Mapped depth',
                                               xdata => [(0..100)],
                                               ydata => $stats->{"gc_vs_$_" . '_cvg'}});
        Graphs::plot_histograms_distributions({outfile=>"$outfiles_prefix.$_" . "_gc_vs_cvg.scaled.$plot_type",
                                               title => "Coverage vs GC plot of $_" . 's',
                                               x_scale => 1, 
                                               x_scale_values => [(30,40,50)],
                                               desc_xvals => "Percentile of $_ sequence ordered by GC content",
                                               desc_xvals_top => '%GC',
                                               desc_yvals => 'Mapped depth',
                                               xdata => [(0..100)],
                                               ydata => $stats->{"gc_vs_$_" . '_cvg'}});
    }

    # Plot: Quality scores by cycle
    foreach (qw(1 2 up)) {

        my %key2string = (1, 'first of pair', 2, 'second of pair', 'up', 'unpaired');
        Graphs::plot_histograms_distributions({outfile=>"$outfiles_prefix.quality_scores_$_.$plot_type",
                                               xdata => [(0 .. $stats->{readlen} - 1)],
                                               ydata => $stats->{'qual_scores_' . $_},
                                               title => "Quality scores distribution, $key2string{$_}",
                                               desc_xvals => 'Cycle',
                                               desc_yvals => 'Quality score',
                                               y_min => 0,
                                               y_max => 41});
    }

    # Plot: reads/targets/baits GC content
    my @read_xvals = (0..$stats->{readlen});
    foreach(@read_xvals){$_ = 100 * $_ / $stats->{readlen}}

    my @zero_to_100 = (0..100);
    

    foreach my $type ('unmapped', 'mapped') {
        my @yvals_1 = @{$stats->{"gc_hist_$type" . '_1'}};
        my @yvals_2 = @{$stats->{"gc_hist_$type" . '_2'}};
        my @data = ({xvals => \@read_xvals, yvals => \@yvals_1, legend => '_1'},
                    {xvals => \@read_xvals, yvals => \@yvals_2, legend => '_2'},
                    {xvals => \@zero_to_100, yvals => $stats->{bait_gc_hist}, legend => 'bait'},
                    {xvals => \@zero_to_100, yvals => $stats->{target_gc_hist}, legend => 'target'},);
        Graphs::plot_stats({outfile=>"$outfiles_prefix.gc_$type.$plot_type",
                            normalize => 1,
                            title => "GC Plot $type reads",
                            desc_xvals => '%GC',
                            desc_yvals => 'Frequency',
                            data => \@data});
    }

    # Plot: insert size
    my %insert_stats = VertRes::Utils::Math->new()->histogram_stats($stats->{insert_size_hist});

    # plot 4 standard devations away from the mean
    my $xmin = max (0, $insert_stats{mean} - 4 * $insert_stats{standard_deviation});
    my $xmax = $insert_stats{mean} + 4 * $insert_stats{standard_deviation};
    my @xvals = ();
    my @yvals = ();

    foreach my $i ($xmin .. $xmax) {
        if (exists $stats->{insert_size_hist}{$i}) {
            push @xvals, $i;
            push @yvals, $stats->{insert_size_hist}{$i};
        }
    }

    Graphs::plot_stats({outfile=>"$outfiles_prefix.insert_size.$plot_type",
                        title => "Insert Size Distribution",
                        desc_xvals => 'Insert Size',
                        desc_yvals => 'Number of read pairs',
                        data => [{xvals => \@xvals, yvals =>\@yvals}]});
}

sub _update_bam_exome_qc_stats {
    my ($self, $intervals, $pileup, $stats, $ref_id, $opts) = @_;
    $stats->{per_ref_sequence}{$ref_id} = {};
    my $ref_id_hash = $stats->{per_ref_sequence}{$ref_id};

    # update stats specific to the given reference sequence
    $ref_id_hash->{bait_bases} = scalar keys %{$pileup->{bait}};
    $ref_id_hash->{bait_near_bases} = scalar keys %{$pileup->{bait_near}};
    $ref_id_hash->{target_bases} = scalar keys %{$pileup->{target}};
    $ref_id_hash->{target_near_bases} = scalar keys %{$pileup->{target_near}};
    $ref_id_hash->{targets} = scalar @{$intervals->{target}->{$ref_id}};
    $ref_id_hash->{baits} = scalar @{$intervals->{bait}->{$ref_id}};

    $ref_id_hash->{bait_bases_mapped} = sum 0, values %{$pileup->{bait}};
    $ref_id_hash->{bait_near_bases_mapped} = sum 0, values %{$pileup->{bait_near}};
    $ref_id_hash->{target_bases_mapped} = sum 0, values %{$pileup->{target}};
    $ref_id_hash->{target_near_bases_mapped} = sum 0, values %{$pileup->{target_near}};

    $ref_id_hash->{target_coverage_mean} = $ref_id_hash->{target_bases_mapped} / $ref_id_hash->{target_bases};
    $ref_id_hash->{bait_coverage_mean} = $ref_id_hash->{bait_bases_mapped} / $ref_id_hash->{bait_bases};

    # update global stats (combination of all reference sequences)
    foreach (qw(bait_bases bait_near_bases target_bases target_near_bases bait_bases_mapped bait_near_bases_mapped target_bases_mapped target_near_bases_mapped targets baits)) {
        $stats->{$_} += $ref_id_hash->{$_};
    }

    # update the per base coverage data
    foreach my $bait_or_target (qw(bait target)) {
        foreach my $cov (values %{$pileup->{$bait_or_target}}) {
            $stats->{$bait_or_target . '_cvg_hist'}{$cov}++;
        }
    }

    # for each bait and target, we want to know:
    # 1) the mean coverage and the %GC
    # 2) whether the min coverage was achieved over the whole target.
    foreach my $bait_or_target qw(bait target) {
        foreach my $interval (@{$intervals->{$bait_or_target}->{$ref_id}}) {
            my $cov_sum = int( 0.5 + sum @{$pileup->{$bait_or_target}}{($interval->[0] .. $interval->[1])} );
            my $mean_cov = int(0.5 + $cov_sum / ($interval->[1] - $interval->[0] + 1) );
            my $gc = int(0.5 + $interval->[2]);        
            $stats->{'gc_vs_' . $bait_or_target . '_cvg'}->[$gc]->{$mean_cov}++;
            my $low_cvg = 1;

            foreach my $i ($interval->[0] .. $interval->[1]) {
                if ($pileup->{$bait_or_target}{$i} >= $opts->{low_cvg}) {
                    $low_cvg = 0;
                    last;
                }                
            } 
            $stats->{'low_cvg_' . $bait_or_target . 's'}++ if $low_cvg;
        }
    }
}

# reads a file of intervals, each line tab delimited:
# reference_id start end
# where start/end coordinates are 1-based
# and returns a hash of the form
# reference_id => ((start1, end1), (start2, end2) ... )
# where each array is sorted according to start position, e.g.
# Order of intervals in input file doesn't matter.  Overlapping
# or adjcent intervals are merged together. e.g.
# 10 100 200
# 10 120 170
# 10 190 250
# 10 251 300
# would be merged into one interval 100-300 on chromosome 10.
#
# Arg1: file of intervals
# Arg2: optional, hash of options.  Currently only option is
# expand_intervals => int. Default = 0.  Adds N bases to the
#                     start and end of each interval.
sub _intervals_file2hash {
    my ($self, $intervals_file, %opts) = @_;
    my %intervals;
    
    $opts{expand_intervals} = 0 unless $opts{expand_intervals};


    # get the all intervals into a hash
    open my $fh, $intervals_file or $self->throw("Cannot open $intervals_file: $!");
    
    while (my $line = <$fh>) { 
        chomp $line;
        my ($ref_id, $start, $end) = split /\t/, $line;
        $intervals{$ref_id} = () unless (exists $intervals{$ref_id});

        $start = max(0, $start - $opts{expand_intervals});
        $end += $opts{expand_intervals};

        push @{$intervals{$ref_id}}, [$start, $end];
    }
    
    close $fh;

    # sort each array and merge overlaps
    foreach my $ref_id (keys %intervals) {
        my $list = $intervals{$ref_id};
        
        # sort by positions, first by coord 1, then by coord 2,  e.g.
        # [1,2] < [1,3] < [2,1] < [2,10]
        @$list = sort {$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @$list;
        
        # merge overalapping intervals.  Also merge adjacent ones,
        # e.g. [1,2], [3,4] becomes [1,4]
        my $i = 0;
        
        while ($i < scalar @$list - 1) {
            # if current interval overlaps next interval, or is adjacent to it
            if ($list->[$i][1] >= $list->[$i+1][0] - 1) {
                $list->[$i + 1] = [$list->[$i][0], max($list->[$i][1], $list->[$i + 1][1])];
                splice @$list, $i, 1;
            }
            else {
                $i++;
            }
        }
    }

    return \%intervals;
}

1;
