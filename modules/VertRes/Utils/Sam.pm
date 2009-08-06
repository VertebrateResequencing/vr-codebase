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
use VertRes::IO;
use VertRes::Wrapper::samtools;
use VertRes::Wrapper::picard;
use HierarchyUtilities;
use VertRes::Parser::sequence_index;
use SamTools;

use base qw(VertRes::Base);


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
 Function: Adds more meta infomation to the header of a sam file.
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
 
           -and-

           either just lane_path => '/path/to/lane' (to calculate the below
           using HierarchyUtilities::lane_info)
           -or-
           all of these key => value pairs:
           lane => run_id, eg. SRR00000
           ref_fa => '/full/path/to/reference.fasta'
           ref_fai => '/full/path/to/reference.fasta.fai'
           ref_name => standard name for the reference, eg. NCBI36
           ref_md5 => the md5sum string for the reference fasta file

=cut

sub add_sam_header {
    my ($self, $raw_sam_file, %args) = @_;
    -s $raw_sam_file || $self->throw("sam file '$raw_sam_file' doesn't exist, can't add header to it!");
    
    # gather information
    my $lane_path = $args{lane_path};
    my %lane_info = %{HierarchyUtilities::lane_info($lane_path)} if $lane_path;
    my $lane = $args{lane} || $lane_info{lane} || $self->throw("lane must be supplied");
    my $ref_fa = $args{ref_fa} || $lane_info{fa_ref} || $self->throw("the reference fasta must be supplied");
    my $ref_fai = $args{ref_fai} || $lane_info{fai_ref} || $self->throw("the reference fai must be supplied");
    my $ref_name = $args{ref_name} || $lane_info{ref_name} || $self->throw("the reference assembly name must be supplied");
    my $ref_md5 = $args{ref_md5} || $lane_info{md5_ref} || $self->throw("the reference md5 must be supplied");
    
    my $seq_index = $args{sequence_index} || '';
    my $parser = VertRes::Parser::sequence_index->new(file => $seq_index) if $seq_index;
    my $individual = $args{sample_name} || $parser->lane_info($lane, 'sample_name') || $lane_info{sample} || $self->throw("Unable to get sample_name for $lane from sequence index '$seq_index' or other supplied args");
    my $library = $args{library} || $parser->lane_info($lane, 'LIBRARY_NAME') || $lane_info{library};
    my $insert_size = $args{insert_size} || $parser->lane_info($lane, 'INSERT_SIZE') || $lane_info{insert_size};
    my $run_name = $args{run_name} || $parser->lane_info($lane, 'run_name') || $self->throw("Unable to get run_name for $lane");
    my $platform = $args{platform} || $parser->lane_info($lane, 'INSTRUMENT_PLATFORM');
    my $centre = $args{centre} || $parser->lane_info($lane, 'CENTER_NAME');
    
    # write the sam header
    my $headed_sam = $raw_sam_file.'.withheader';
    my $header_lines = 0;
    open(my $shfh, '>', $headed_sam) or $self->throw("Cannot create $headed_sam: $!");
    print $shfh "\@HD\tVN:1.0\tSO:coordinate\n";
    $header_lines++;
    
    open(my $faifh, $ref_fai) or $self->throw("Can't open ref fai file: $ref_fai $!");
    while (<$faifh>) {
        chomp;
        my @s = split;
        print $shfh "\@SQ\tSN:$s[0]\tLN:$s[1]\tAS:$ref_name\tM5:$ref_md5\tUR:file:$ref_fa\n";
        $header_lines++;
    }
    close($faifh);
    
    print $shfh "\@RG\tID:$lane\tPU:$run_name\tLB:$library\tSM:$individual";
    print $shfh "\tPI:$insert_size" if (defined $insert_size);
    print $shfh "\tCN:$centre" if (defined $centre);
    print $shfh "\tPL:$platform" if (defined $platform);
    
    print $shfh "\n";
    $header_lines++;
    
    # combine the header with the raw sam file
    open(my $rsfh, '<', $raw_sam_file) or $self->throw("Couldn't open raw sam file '$raw_sam_file'");
    while (<$rsfh>) {
        print $shfh $_;
    }
    close($rsfh);
    close($shfh);
    
    # check and mv
    my $io = VertRes::IO->new(file => $raw_sam_file);
    my $raw_lines = $io->num_lines();
    my $expected_lines = $raw_lines + $header_lines;
    $io->file($headed_sam);
    my $actual_lines = $io->num_lines();
    
    if ($expected_lines == $actual_lines) {
        move($headed_sam, $raw_sam_file) || $self->throw("Failed to move $headed_sam to $raw_sam_file: $!");
        return 1;
    }
    else {
        unlink($headed_sam);
        $self->throw("Failed to prepend sam header to sam file '$raw_sam_file'");
    }
}

=head2 sam_to_fixed_sorted_bam

 Title   : sam_to_fixed_sorted_bam
 Usage   : $obj->sam_to_fixed_sorted_bam('in.sam', 'sorted.bam', 'ref.fai');
 Function: Converts a sam file to a mate-fixed and sorted bam file.
 Returns : boolean (true on success)
 Args    : starting sam file, output name for bam file. The reference.fai is
           only needed if add_sam_header() hasn't already been called on the
           in.sam. Optionally, args to pass to VertRes::Wrapper::samtools (eg.
           quiet => 1).

=cut

sub sam_to_fixed_sorted_bam {
    my ($self, $in_sam, $out_bam, $ref_fai, @args) = @_;
    
    my $wrapper = VertRes::Wrapper::samtools->new(verbose => $self->verbose, @args);
    $wrapper->sam_to_fixed_sorted_bam($in_sam, $out_bam, $ref_fai);
    
    return $wrapper->run_status() >= 1;
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
    my $wrapper = VertRes::Wrapper::picard->new(quiet => 1,
                                                validation_stringency => 'silent',
                                                tmp_dir => $io->tempdir());
    
    $wrapper->MergeSamFiles($out_bam, @in_bams);
    
    return $wrapper->run_status() >= 1;
}

=head2 stats

 Title   : stats
 Usage   : $obj->stats('in.bam', 'in2.bam', ...);
 Function: Generate both flagstat and bamstat files for the given bam files
 Returns : boolean (true on success)
 Args    : list of bam files (for each, a .bamstat and .flagstat file will be
           created)

=cut

sub stats {
    my ($self, @in_bams) = @_;
    
    my $wrapper = VertRes::Wrapper::samtools->new(verbose => $self->verbose);
    my $all_ok = 1;
    foreach my $in_bam (@in_bams) {
        #  $cmd .= qq[ $SAMTOOLS index $lane_dir/rmdup.bam;]; ??? old mapping code used to index the rmdup.bam (only)
        $wrapper->flagstat($in_bam, $in_bam.'.flagstat');
        $all_ok = 0 unless $wrapper->run_status() >= 1;
        
        my $bamstat = $in_bam.'.bamstat';
        SamTools::bam_stat($in_bam, $bamstat);
        $all_ok = 0 unless -s $bamstat;
    }
    
    return $all_ok;
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

1;
