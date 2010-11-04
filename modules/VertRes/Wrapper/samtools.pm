=head1 NAME

VertRes::Wrapper::samtools - wrapper for samtools

=head1 SYNOPSIS

use VertRes::Wrapper::samtools;

my $wrapper = VertRes::Wrapper::samtools->new();

# Now you can call a method corrseponding to a samtools command. (except for
# import, which must be called as sam_import()).
# All methods take the filenames as a list to start with, followed by a hash
# of any options the samtools command understands.

$wrapper->sam_import('in.ref_list', 'in.sam', 'out.bam');
$wrapper->view('in.[bs]am', 'output.file', regions => [], %options);
# eg. to convert sam to bam:
$wrapper->view('in.sam', 'out.bam', b => 1, t => 'ref.fai');
$wrapper->index('input.bam', 'output.bai');
# etc.

# There are also some special convienience methods that do multi-step
# operations:
$wrapper->sam_to_fixed_sorted_bam('in.sam', 'out.bam', 'ref.fai');

=head1 DESCRIPTION

Runs samtools in a nice way.

*** not all samtools commands have been wrapped yet...

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Wrapper::samtools;

use strict;
use warnings;
use File::Copy;
use VertRes::IO;
use VertRes::Utils::FileSystem;
BEGIN { eval 'use VertRes::Utils::Sam'; }

use base qw(VertRes::Wrapper::WrapperI);


=head2 new

 Title   : new
 Usage   : my $wrapper = VertRes::Wrapper::samtools->new();
 Function: Create a VertRes::Wrapper::samtools object.
 Returns : VertRes::Wrapper::samtools object
 Args    : quiet   => boolean
           exe     => string (default 'samtools' - the first one on your path
                              will get used)

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    $self->{base_exe} = $self->exe() || 'samtools';
    
    return $self;
}

=head2 sam_import

 Title   : sam_import
 Usage   : $wrapper->sam_import('in.ref_list', 'in.sam', 'out.bam');
 Function: import...
 Returns : n/a
 Args    : list of file paths

=cut

sub sam_import {
    my ($self, $in_ref_list, $in_sam, $out_bam) = @_;
    
    $self->exe($self->{base_exe}.' import');
    
    $self->switches([]);
    $self->params([]);
    
    $self->register_output_file_to_check($out_bam);
    
    return $self->run($in_ref_list, $in_sam, $out_bam);
}

=head2 view

 Title   : view
 Usage   : $wrapper->view('in.[bs]am', 'output.file', regions => [], %options);
 Function: view...
 Returns : n/a
 Args    : list of file paths, options as a hash including regions as a key
           and the value as an array ref. If using the 'open' run_method, the
           second arg should be undef. Special option append => 1 can be used
           to append to output file instead of overwriting it; only makes sense
           when outputting sam format.

=cut

sub view {
    my ($self, $in_sam, $out_file, %options) = @_;
    
    $self->exe($self->{base_exe}.' view');
    
    my @file_and_regions = ($in_sam);
    if (defined $options{regions} && ref($options{regions}) eq 'ARRAY') {
        my $regions = delete $options{regions};
        push(@file_and_regions, @{$regions});
    }
    
    if ($out_file) {
        my $append = delete $options{append};
        my $redirect = $append ? ' >> ' : ' > ';
        push(@file_and_regions, $redirect.$out_file);
        $self->register_output_file_to_check($out_file);
    }
    
    $self->switches([qw(b h H S)]);
    $self->params([qw(t o f F q)]);
    $self->_set_params_and_switches_from_args(%options);
    
    
    return $self->run(@file_and_regions);
}

=head2 sort

 Title   : sort
 Usage   : $wrapper->sort('in.bam', 'out_prefix', %options);
 Function: sort...
 Returns : n/a
 Args    : list of file paths, options as a hash

=cut

sub sort {
    my ($self, $in_bam, $out_prefix, %options) = @_;
    
    $self->exe($self->{base_exe}.' sort');
    
    $self->switches([qw(n)]);
    $self->params([qw(m)]);
    $self->_set_params_and_switches_from_args(%options);
    
    $self->register_output_file_to_check($out_prefix.'.bam');
    
    return $self->run($in_bam, $out_prefix);
}

=head2 merge

 Title   : merge
 Usage   : $wrapper->merge('out.bam', \@in_bams, %options);
 Function: merge... 
 Returns : n/a
 Args    : paths to out.bam and array ref of input bams, options as a hash

=cut

sub merge {
    my ($self, $out_bam, $in_bams, %options) = @_;
    
    $self->exe($self->{base_exe}.' merge');
    
    $self->switches([qw(n)]);
    $self->params([]);
    $self->_set_params_and_switches_from_args(%options);
    
    $self->register_output_file_to_check($out_bam);
    
    return $self->run($out_bam, @{$in_bams});
}

=head2 merge_and_check

 Title   : merge_and_check
 Usage   : $wrapper->merge_and_check('out.bam', \@in_bams, %options);
 Function: Merges multiple bam files together and checks the output merged bam
           isn't truncated.
 Returns : n/a
 Args    : paths to out.bam and array ref of input bams, options as a hash

=cut

sub merge_and_check {
    my ($self, $out_bam, $in_bams, %options) = @_;
    
    my $orig_run_method = $self->run_method;
    
    # do the merge
    $self->run_method('system');
    $self->merge($out_bam.'.tmp', $in_bams, %options);
    $self->throw("failed during the merge step, giving up for now") unless $self->run_status >= 1;
    
    # find out our expectation
    my $su = VertRes::Utils::Sam->new;
    my $bam_count = 0;
    foreach my $bam_file (@{$in_bams}) {
        $bam_count += $su->num_bam_records($bam_file);
    }
    
    # find out how many lines are in the merged bam file
    my $merge_count = $su->num_bam_records($out_bam.'.tmp');
    
    # check for truncation
    if ($merge_count >= $bam_count) {
        move("$out_bam.tmp", $out_bam) || $self->throw("Failed to move $out_bam.tmp to $out_bam: $!");
        $self->_set_run_status(2);
    }
    else {
        $self->warn("$out_bam.tmp is bad (only $merge_count lines vs $bam_count), will unlink it");
        $self->_set_run_status(-1);
        unlink("$out_bam.tmp");
    }
    
    $self->run_method($orig_run_method);
    return;
}

=head2 pileup

 Title   : pileup
 Usage   : $wrapper->pileup('in.[sb]am', 'out.file', %options);
 Function: pileup...
 Returns : n/a
 Args    : list of file paths, options as a hash. If using the 'open' run_method
           you should set the second arg to undef.

=cut

sub pileup {
    my ($self, $in_file, $out_file, %options) = @_;
    
    $self->exe($self->{base_exe}.' pileup');
    
    $self->switches([qw(s i c g S a 2)]);
    $self->params([qw(m M t l f T N r G I)]);
    $self->_set_params_and_switches_from_args(%options);
    
    my @files = ($in_file);
    if ($out_file) {
        push(@files, ' > '.$out_file);
        $self->register_output_file_to_check($out_file);
    }
    
    return $self->run(@files);
}

=head2 fillmd

 Title   : fillmd
 Usage   : $wrapper->fillmd('in.bam', 'ref.fa', 'out.sam', %options);
 Function: fillmd...
 Returns : n/a
 Args    : list of file paths, options as a hash. If using the 'open' run_method
           you should set the third arg to undef.

=cut

sub fillmd {
    my ($self, $in_file, $ref, $out_file, %options) = @_;
    
    $self->exe($self->{base_exe}.' fillmd');
    
    $self->switches([qw(e u b S r)]);
    $self->params([]);
    $self->_set_params_and_switches_from_args(%options);
    
    my @files = ($in_file, $ref);
    if ($out_file) {
        push(@files, ' > '.$out_file);
        $self->register_output_file_to_check($out_file);
    }
    
    return $self->run(@files);
}
sub calmd; *calmd = \&fillmd; 

=head2 calmd_and_check

 Title   : calmd_and_check
 Usage   : $wrapper->calmd_and_check('in.bam', 'ref.fa', 'out.bam', %ops);
 Function: Runs calmd on a bam and creates an output bam (as opposed to the
           normal sam), checking the output bam is complete.
 Returns : n/a
 Args    : in and out bam paths, optional args as understood by calmd, like
           r => 1.

=cut

sub calmd_and_check {
    my ($self, $in_bam, $ref, $out_bam, %options) = @_;
    
    my $orig_run_method = $self->run_method;
    
    # do the calmd
    $self->run_method('open');
    my $tmp_bam = $out_bam.'.tmp';
    my $fh = $self->fillmd($in_bam, $ref, undef, %options);
    $self->throw("failed during the calmd step, giving up for now") unless $self->run_status >= 1;
    
    # convert to bam
    $self->run_method('open_to');
    $self->view($fh, $tmp_bam, S => 1, b => 1);
    $self->throw("failed during the view step, giving up for now") unless $self->run_status >= 1;
    
    # find out our expectation
    my $su = VertRes::Utils::Sam->new;
    my $bam_count = $su->num_bam_records($in_bam);
    
    # find out how many lines are in the merged bam file
    my $merge_count = $su->num_bam_records($tmp_bam);
    
    # check for truncation
    if ($merge_count >= $bam_count) {
        move($tmp_bam, $out_bam) || $self->throw("Failed to move $tmp_bam to $out_bam: $!");
        $self->_set_run_status(2);
    }
    else {
        $self->warn("$tmp_bam is bad (only $merge_count lines vs $bam_count), will unlink it");
        $self->_set_run_status(-1);
        unlink($tmp_bam);
    }
    
    $self->run_method($orig_run_method);
    return;
}

=head2 index

 Title   : index
 Usage   : $wrapper->index('in.bam', 'out.index');
 Function: index...
 Returns : n/a
 Args    : list of file paths

=cut

sub index {
    my ($self, $in_bam, $out_index) = @_;
    
    $self->exe($self->{base_exe}.' index');
    
    $self->switches([]);
    $self->params([]);
    
    $self->register_output_file_to_check($out_index);
    
    return $self->run($in_bam, $out_index);
}

=head2 faidx

 Title   : faidx
 Usage   : $wrapper->faidx('in.fa');
           my ($seq1, $seq2) = $wrapper->faidx('ref.fa', '1:1-5', '2:1-5');
 Function: faidx...
 Returns : no regions: n/a (creates 'in.fa.fai')
           regions: list of strings (reference bases, undef if region invalid)
 Args    : path to input fasta (only, to generate .fai file), optionally one or
           more region strings to get back the reference bases in that region.
           region specifications are strings like '1' for all of chr1, '1:20'
           for chr1 position 20 onward, and '1:20-21' to get the 2 bases at
           positions 20 and 21. Positions are 1-based.

=cut

sub faidx {
    my ($self, $in_fa, @regions) = @_;
    
    $self->exe($self->{base_exe}.' faidx');
    
    $self->switches([]);
    $self->params([]);
    
    
    my @files = ($in_fa);
    if (@regions) {
        push(@files, @regions);
        my $orig_method = $self->run_method;
        $self->run_method('open');
        
        my $fh = $self->run(@files);
        my @seqs;
        my $seq;
        while (<$fh>) {
            chomp;
            if (/^>/) {
                push(@seqs, $seq) if $seq;
                $seq = '';
                next;
            }
            $seq .= $_;
        }
        close($fh);
        push(@seqs, $seq) if $seq;
        
        $self->run_method($orig_method);
        return @seqs;
    }
    else {
        $self->register_output_file_to_check($in_fa.'.fai');
        return $self->run(@files);
    }
}

=head2 fixmate

 Title   : fixmate
 Usage   : $wrapper->fixmate('in.namesorted.bam', 'out.namesorted.bam');
 Function: fixmate...
 Returns : n/a
 Args    : list of file paths

=cut

sub fixmate {
    my ($self, $in_bam, $out_bam) = @_;
    
    $self->exe($self->{base_exe}.' fixmate');
    
    $self->switches([]);
    $self->params([]);
    
    $self->register_output_file_to_check($out_bam) unless $out_bam eq '-';
    
    return $self->run($in_bam, $out_bam);
}

=head2 rmdup

 Title   : rmdup
 Usage   : $wrapper->rmdup('sorted_in.bam', 'out.bam', %options);
 Function: rmdup...
 Returns : n/a
 Args    : paths to input and output bams, optionally options understood by
           rmdup, most significantly s => 1 to set single-ended mode

=cut

sub rmdup {
    my ($self, $in_bam, $out_bam, %options) = @_;
    
    $self->exe($self->{base_exe}.' rmdup');
    
    $self->switches([qw(s S)]);
    $self->params([]);
    $self->_set_params_and_switches_from_args(%options);
    
    $self->register_output_file_to_check($out_bam);
    
    return $self->run($in_bam, $out_bam);
}

=head2 flagstat

 Title   : flagstat
 Usage   : $wrapper->flagstat('in.bam', 'out.bam.flagstat');
 Function: flagstat...
 Returns : n/a
 Args    : paths to input and output bams. output can be excluded if you're
           piping out with run_method('open')

=cut

sub flagstat {
    my ($self, $in_bam, $out_bam) = @_;
    
    $self->exe($self->{base_exe}.' flagstat');
    
    $self->switches([]);
    $self->params([]);
    
    $self->register_output_file_to_check($out_bam) if $out_bam;
    my @args = ($in_bam);
    push(@args, ' > '.$out_bam) if $out_bam;
    
    return $self->run(@args);
}

=head2 sam_to_fixed_sorted_bam

 Title   : sam_to_fixed_sorted_bam
 Usage   : $wrapper->sam_to_fixed_sorted_bam('in.sam', 'out.bam', 'ref.fai');
 Function: Given a sam file, converts it to a name-sorted bam file, runs
           fixmate on it, and re-sorts it by position. Checks the output file
           is not truncated.
 Returns : n/a
 Args    : list of file paths. The last, ref.fai, can be left out if the sam
           file has a header with @SQ lines.

=cut

sub sam_to_fixed_sorted_bam {
    my ($self, $in_sam, $out_bam, $ref_fai) = @_;
    
    my $orig_run_method = $self->run_method;
    
    $self->run_method('open');
    my $fh = $self->view($in_sam, undef, b => 1, S => 1, $ref_fai ? (t => $ref_fai) : ());
    $self->throw("failed during the view step, giving up for now") unless $self->run_status >= 1;
    $fh || $self->throw("failed to get a filehandle from the view step, giving up for now");
    
    my $io = VertRes::IO->new();
    my $fsu = VertRes::Utils::FileSystem->new();
    my $temp_dir = $fsu->tempdir();
    my $name_sorted_file = $fsu->catfile($temp_dir, 'namesorted');
    
    $self->run_method('open_to');
    $self->sort($fh, $name_sorted_file, n => 1);
    $self->throw("failed during the name sort step, giving up for now") unless $self->run_status >= 1;
    
    $self->run_method('open');
    undef($fh);
    $fh = $self->fixmate($name_sorted_file.'.bam', '-');
    $self->throw("failed during the fixmate step, giving up for now") unless $self->run_status >= 1;
    $fh || $self->throw("failed to get a filehandle from the fixmate step, giving up for now");
    
    $self->run_method('open_to');
    $out_bam =~ s/\.bam$//;
    my $tmp_bam = $out_bam.'_tmp';
    $self->sort($fh, $tmp_bam, n => 0);
    
    # check the output isn't truncated, or unlink it
    $self->run_method('open');
    undef($fh);
    $fh = $self->view($tmp_bam.'.bam', undef, h => 1);
    my $bam_count = 0;
    while (<$fh>) {
        $bam_count++;
    }
    close($fh);
    $io->file($in_sam);
    my $sam_count = $io->num_lines();
    if ($bam_count >= $sam_count) {
        move("$tmp_bam.bam", "$out_bam.bam") || $self->throw("Failed to move $tmp_bam.bam to $out_bam.bam: $!");;
        $self->_set_run_status(2);
    }
    else {
        $self->warn("$tmp_bam.bam is bad (only $bam_count lines vs $sam_count), will unlink it");
        $self->_set_run_status(-1);
        unlink("$tmp_bam.bam");
    }
    
    $self->run_method($orig_run_method);
    return;
}

1;
