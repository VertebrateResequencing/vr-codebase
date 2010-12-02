=head1 NAME

VertRes::Wrapper::picard - wrapper for picard tools

=head1 SYNOPSIS

use VertRes::Wrapper::picard;

my $wrapper = VertRes::Wrapper::picard->new(validation_stringency => 'silent');

# Now you can call a method corrseponding to a picard tools jar.
# All methods take the filenames as a list to start with, followed by a hash
# of any additional options the picard jar understands.
$wrapper->MarkDuplicates($in_bam, $out_bam);
$wrapper->MergeSamFiles($out_bam, @bams_to_mergs);

# not yet wrapped:
$wrapper->CreateSequenceDictionary();
$wrapper->ValidateSamFile();
$wrapper->ViewSam();

=head1 DESCRIPTION

Runs picard tools in a nice way.
NB: you need to reserve at least 5GB to run this, or the java engine might fail
to start!

* not all picard tools have been wrapped yet*

For default "exe" path assumes you have the env variable PICARD pointing to the
directory containing the picard .jar files etc.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Wrapper::picard;

use strict;
use warnings;
use File::Copy;
use VertRes::Utils::FileSystem;
use VertRes::Wrapper::samtools;
use VertRes::Utils::Sam;
use VertRes::Parser::sam;

use base qw(VertRes::Wrapper::WrapperI);

my $fsu = VertRes::Utils::FileSystem->new();
our $DEFAULT_PICARD_DIR = $ENV{PICARD} || die "PICARD environment variable not set\n";

=head2 new

 Title   : new
 Usage   : my $wrapper = VertRes::Wrapper::picard->new();
 Function: Create a VertRes::Wrapper::picard object.
 Returns : VertRes::Wrapper::picard object
 Args    : quiet   => boolean
           exe     => string (full path to the location of the picard tools jar
                              files; a TEAM145 default exists)
           java_memory => int (the amount of memory in MB to give java; default
                               5000)
           validation_stringency => STRICT|LENIENT|SILENT (silent by default,
                                    overriden if VALIDATION_STRINGENCY is set
                                    directly in any other method call)
           max_records_in_ram => int (default 100x java_memory)
           tmp_dir => /tmp/dir (VertRes::Utils::FileSystem->tempdir by default,
                                overriden if TMP_DIR is set directly in any
                                other method call. Any supplied directory is
                                used as the root for a new directory that will
                                be auto-deleted)

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    $self->{picard_dir} = $self->exe() || $DEFAULT_PICARD_DIR;
    my $java_mem = delete $self->{java_memory} || 5000;
    $self->{base_exe} = "java -Xmx${java_mem}m -jar ";
    
    # our bsub jobs will get killed if we don't select high-mem machines
    $self->bsub_options(M => ($java_mem * 1000), R => "'select[mem>$java_mem] rusage[mem=$java_mem]'");
    
    my $stringency = delete $self->{validation_stringency} || 'silent';
    $self->{_default_validation_stringency} = uc($stringency);
    
    my $temp_dir = delete $self->{tmp_dir};
    $temp_dir = $fsu->tempdir($temp_dir ? (DIR => $temp_dir) : ());
    $self->{_default_tmp_dir} = $temp_dir;
    
    $self->{_default_max_records_in_ram} = delete $self->{max_records_in_ram} || ($java_mem * 100);
    
    return $self;
}

sub _handle_common_params {
    my ($self, $params) = @_;
    
    foreach my $key (qw(I INPUT O OUTPUT)) {
        delete $params->{$key};
    }
    
    unless (defined $params->{VALIDATION_STRINGENCY}) {
        $params->{VALIDATION_STRINGENCY} = $self->{_default_validation_stringency};
    }
    
    unless (defined $params->{MAX_RECORDS_IN_RAM}) {
        $params->{MAX_RECORDS_IN_RAM} = $self->{_default_max_records_in_ram};
    }
    
    unless (defined $params->{TMP_DIR}) {
        $params->{TMP_DIR} = $self->{_default_tmp_dir};
    }
    else {
        $params->{TMP_DIR} = $fsu->tempdir(DIR => $params->{TMP_DIR});
    }
}

=head2 MergeSamFiles

 Title   : MergeSamFiles
 Usage   : $wrapper->MergeSamFiles('out.bam', @in_bams, %options);
 Function: MergeSamFiles... (can merge both sam and bam files)
 Returns : n/a
 Args    : list of file paths (output bam, input bams), followed by a hash of
           options understood by MergeSamFiles, eg. VALIDATION_STRINGENCY =>
           'SILENT'. (case matters: must be uppercase)

=cut

sub MergeSamFiles {
    my ($self, $out_bam, @args) = @_;
    
    $self->exe($self->{base_exe}.$fsu->catfile($self->{picard_dir}, 'MergeSamFiles.jar'));
    
    $self->switches([]);
    $self->params([qw(SORT_ORDER SO ASSUME_SORTED AS TMP_DIR VERBOSITY QUIET
                      VALIDATION_STRINGENCY COMPRESSION_LEVEL OPTIONS_FILE
                      MAX_RECORDS_IN_RAM)]);
    
    my (@in_bams, @params);
    foreach my $arg (@args) {
        if (-e $arg && ! @params) {
            push(@in_bams, " I=$arg");
        }
        else {
            push(@params, $arg);
        }
    }
    my @file_args = (" O=$out_bam", @in_bams);
    
    my %params = @params;
    $self->_handle_common_params(\%params);
    
    $self->register_output_file_to_check($out_bam);
    $self->_set_params_and_switches_from_args(%params);
    
    return $self->run(@file_args);
}

=head2 merge_and_check

 Title   : merge_and_check
 Usage   : $wrapper->merge_and_check('out.bam', \@in_bams, %options);
 Function: Merges multiple sam/bam files together and checks the output merged
           file isn't truncated.
 Returns : n/a
 Args    : output bam, array ref of input bams, options hash as understood by
           MergeSamFiles()

=cut

sub merge_and_check {
    my ($self, $out_bam, $in_bams, %opts) = @_;
    
    my $orig_run_method = $self->run_method;
    
    # do the merge
    $self->run_method('system');
    my $tmp_out = $out_bam;
    $tmp_out =~ s/\.bam$/.tmp.bam/; # won't work unless the final suffix is .bam
    $self->MergeSamFiles($tmp_out, @{$in_bams}, %opts);
    $self->throw("failed during the merge step, giving up for now") unless $self->run_status >= 1;
    
    # find out our expectation
    my $su = VertRes::Utils::Sam->new();
    my $bam_count = 0;
    foreach my $bam_file (@{$in_bams}) {
        if ($bam_file =~ /\.sam$/) {
            open(my $fh, $bam_file) || $self->throw("Could not open sam file '$bam_file'");
            while (<$fh>) {
                $bam_count++;
            }
            close($fh);
        }
        else {
            $bam_count += $su->num_bam_records($bam_file);
        }
    }
    
    # find out how many lines are in the merged bam file
    my $merge_count = 0;
    if ($out_bam =~ /\.sam$/) {
        open(my $fh, $tmp_out) || $self->throw("Could not open sam file '$tmp_out'");
        while (<$fh>) {
            $merge_count++;
        }
        close($fh);
    }
    else {
        $merge_count = $su->num_bam_records($tmp_out);
    }
    
    # check for truncation
    if ($merge_count == $bam_count) {
        move($tmp_out, $out_bam) || $self->throw("Failed to move $tmp_out to $out_bam: $!");
        $self->_set_run_status(2);
    }
    else {
        $self->warn("$tmp_out is bad (only $merge_count lines vs $bam_count), will unlink it");
        $self->_set_run_status(-1);
        unlink("$tmp_out");
    }
    
    $self->run_method($orig_run_method);
    return;
}

=head2 FixMateInformation

 Title   : FixMateInformation
 Usage   : $wrapper->FixMateInformation($in_bam, $out_bam, %options);
 Function: Can take a bam of any sort order and make a coordinate-sorted bam
           with the mate information fixed.
 Returns : n/a
 Args    : list of file paths (input bam, output bam), followed by a hash of
           options understood by FixMateInformation, eg. VALIDATION_STRINGENCY
           => 'SILENT'. (case matters: must be uppercase)

=cut

sub FixMateInformation {
    my ($self, $in_bam, $out_bam, %args) = @_;
    
    #java -Djava.io.tmpdir=/path/to/tmpdir \  [this argument recommended when dealing with large input]
    #   -jar FixMateInformation.jar \
    #   INPUT=<input1.bam> \
    #   OUTPUT=<fixedBam.bam> \
    #   SO=coordinate \
    #   VALIDATION_STRINGENCY=SILENT
    
    $self->exe($self->{base_exe}.$fsu->catfile($self->{picard_dir}, 'FixMateInformation.jar'));
    
    $self->switches([]);
    $self->params([qw(TMP_DIR VERBOSITY QUIET VALIDATION_STRINGENCY
                   COMPRESSION_LEVEL SORT_ORDER MAX_RECORDS_IN_RAM)]);
    
    my @file_args = (" I=$in_bam", " O=$out_bam");
    $self->_handle_common_params(\%args);
    unless (defined $args{SO} || defined $args{SORT_ORDER}) {
        $args{SORT_ORDER} = 'coordinate';
    }
    
    $self->register_output_file_to_check($out_bam);
    $self->_set_params_and_switches_from_args(%args);
    
    return $self->run(@file_args);
}

=head2 MarkDuplicates

 Title   : MarkDuplicates
 Usage   : $wrapper->MarkDuplicates($in_bam, $out_bam, %options);
 Function: MarkDuplicates...
 Returns : n/a
 Args    : list of file paths (input bam, output bam), followed by a hash of
           options understood by MarkDuplicates, eg. VALIDATION_STRINGENCY =>
           'SILENT'. (case matters: must be uppercase)
           By default, METRICS_FILE (M) goes to /dev/null

=cut

sub MarkDuplicates {
    my ($self, $in_bam, $out_bam, %args) = @_;
    
    $self->exe($self->{base_exe}.$fsu->catfile($self->{picard_dir}, 'MarkDuplicates.jar'));
    
    $self->switches([]);
    $self->params([qw(METRICS_FILE M TMP_DIR VERBOSITY QUIET
                      VALIDATION_STRINGENCY COMPRESSION_LEVEL OPTIONS_FILE
                      MAX_RECORDS_IN_RAM)]);
    
    my @file_args = (" I=$in_bam", " O=$out_bam");
    $self->_handle_common_params(\%args);
    unless (defined $args{M} || defined $args{METRICS_FILE}) {
        $args{M} = '/dev/null';
    }
    
    $self->register_output_file_to_check($out_bam);
    $self->_set_params_and_switches_from_args(%args);
    
    return $self->run(@file_args);
}

=head2 markdup

 Title   : markdup
 Usage   : $wrapper->markdup($in_bam, $out_bam, %options);
 Function: Calls MarkDuplicates() and checks the output bam is sane.
 Returns : n/a
 Args    : list of file paths (input bam, output bam), followed by a hash of
           options understood by MarkDuplicates().

=cut

sub markdup {
    my ($self, $in_bam, $out_bam, %args) = @_;
    
    # MarkDuplicates
    my $tmp_bam = $out_bam.'.tmp.bam';
    my $orig_run_method = $self->run_method;
    $self->run_method('system');
    $self->MarkDuplicates($in_bam, $tmp_bam, %args);
    $self->run_method($orig_run_method);
    $self->throw("failed during the MarkDuplicates step, giving up") unless $self->run_status >= 1;
    
    # check the output isn't truncated, or unlink it
    my $su = VertRes::Utils::Sam->new();
    my $bam_count = $su->num_bam_lines($tmp_bam);
    my $expected_count = $su->num_bam_lines($in_bam);
    if ($bam_count >= $expected_count) {
        move($tmp_bam, $out_bam) || $self->throw("Failed to move $tmp_bam to $out_bam: $!");
        $self->_set_run_status(2);
    }
    else {
        $self->warn("$tmp_bam is bad ($bam_count lines vs $expected_count), will unlink it");
        $self->_set_run_status(-1);
        unlink("$tmp_bam");
    }
    
    return;
}

=head2 rmdup

 Title   : rmdup
 Usage   : $wrapper->rmdup($in_bam, $out_bam, %options);
 Function: rmdup creates an output bam file with duplicate reads physically
           removed, based on the output of MarkDuplicates().
 Returns : n/a
 Args    : list of file paths (input bam, output bam), followed by a hash of
           options understood by MarkDuplicates().

=cut

sub rmdup {
    my ($self, $in_bam, $out_bam, %args) = @_;
    
    # MarkDuplicates
    my $temp_dir = $fsu->tempdir;
    my $marked_bam = $fsu->catfile($temp_dir, 'marked.bam');
    my $orig_run_method = $self->run_method;
    $self->run_method('system');
    $self->MarkDuplicates($in_bam, $marked_bam, %args);
    $self->run_method($orig_run_method);
    $self->throw("failed during the MarkDuplicates step, giving up") unless $self->run_status >= 1;
    
    # use samtools to run through the bam, filter out dups, and create a new bam
    my $st = VertRes::Wrapper::samtools->new(quiet => 1);
    $st->run_method('open');
    my $fh = $st->view($marked_bam, undef, h => 1);
    $self->throw("failed during the first view step, giving up") unless $st->run_status >= 1;
    $fh || $self->throw("failed to get a filehandle from the view step, giving up");
    
    my $sp = VertRes::Parser::sam->new();
    
    my $total_lines = 0;
    my $dup_lines = 0;
    my $filtered_sam = $fsu->catfile($temp_dir, 'filtered.sam');
    open(my $fsfh, '>', $filtered_sam) || $self->throw("Could not write to $filtered_sam: $!");
    while (<$fh>) {
        if (/^@/) {
            print $fsfh $_;
            next;
        }
        
        $total_lines++;
        
        my (undef, $flag) = split(qr/\t/, $_);
        if ($sp->is_duplicate($flag)) {
            $dup_lines++;
        }
        else {
            print $fsfh $_;
        }
    }
    close($fh);
    close($fsfh);
    
    $st->run_method('system');
    my $tmp_bam = $out_bam.'_tmp';
    $st->view($filtered_sam, $tmp_bam, h => 1, S => 1, b => 1);
    $self->throw("failed during the second view step, giving up") unless $st->run_status >= 1;
    
    # check the output isn't truncated, or unlink it
    $st->run_method('open');
    undef($fh);
    $fh = $st->view($tmp_bam, undef, h => 1);
    my $bam_count = 0;
    while (<$fh>) {
        $bam_count++;
    }
    close($fh);
    my $expected_count = $total_lines - $dup_lines;
    if ($bam_count >= $expected_count) { # >= because it might have 1 or 2 extra header lines
        move($tmp_bam, $out_bam) || $self->throw("Failed to move $tmp_bam to $out_bam: $!");
        $self->_set_run_status(2);
    }
    else {
        $self->warn("$tmp_bam.bam is bad ($bam_count lines vs $expected_count), will unlink it");
        $self->_set_run_status(-1);
        unlink("$tmp_bam.bam");
    }
    
    return;
}

sub _pre_run {
    my $self = shift;
    $self->_set_params_string(join => '=');
    return @_;
}

sub run {
    my $self = shift;
    
    # refuses to be quiet, so force the issue
    if ($self->quiet) {
        my $run_method = $self->run_method;
        $self->run_method('open');
        my ($fh) = $self->_run(@_);
        while (<$fh>) {
            next;
        }
        close($fh);
        $self->_post_run();
        $self->run_method($run_method);
    }
    else {
        return $self->SUPER::run(@_);
    }
}

1;
