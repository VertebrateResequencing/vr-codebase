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
use VertRes::IO;

use base qw(VertRes::Wrapper::WrapperI);


=head2 new

 Title   : new
 Usage   : my $wrapper = VertRes::Wrapper::samtools->new();
 Function: Create a VertRes::Wrapper::samtools object.
 Returns : VertRes::Wrapper::samtools object
 Args    : quiet   => boolean

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args, exe => 'samtools');
    
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
    
    $self->exe('samtools import');
    
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
           second arg should be undef.

=cut

sub view {
    my ($self, $in_sam, $out_file, %options) = @_;
    
    $self->exe('samtools view');
    
    my @file_and_regions = ($in_sam);
    if (defined $options{regions} && ref($options{regions}) eq 'ARRAY') {
        my $regions = delete $options{regions};
        push(@file_and_regions, @{$regions});
    }
    
    if ($out_file) {
        push(@file_and_regions, ' > '.$out_file);
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
    
    $self->exe('samtools sort');
    
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
    
    $self->exe('samtools merge');
    
    $self->switches([qw(n)]);
    $self->params([]);
    $self->_set_params_and_switches_from_args(%options);
    
    $self->register_output_file_to_check($out_bam);
    
    return $self->run($out_bam, @{$in_bams});
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
    
    $self->exe('samtools pileup');
    
    $self->switches([qw(s i c g)]);
    $self->params([qw(m t l f T N r G I)]);
    $self->_set_params_and_switches_from_args(%options);
    
    my @files = ($in_file);
    if ($out_file) {
        push(@files, ' > '.$out_file);
        $self->register_output_file_to_check($out_file);
    }
    
    return $self->run(@files);
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
    
    $self->exe('samtools index');
    
    $self->switches([]);
    $self->params([]);
    
    $self->register_output_file_to_check($out_index);
    
    return $self->run($in_bam, $out_index);
}

=head2 faidx

 Title   : faidx
 Usage   : $wrapper->faidx('in.fa');
 Function: faidx...
 Returns : n/a (creates 'in.fa.fai')
 Args    : paths to input fasta

=cut

sub faidx {
    my ($self, $in_fa) = @_;
    
    $self->exe('samtools faidx');
    
    $self->switches([]);
    $self->params([]);
    
    $self->register_output_file_to_check($in_fa.'.fai');
    
    return $self->run($in_fa);
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
    
    $self->exe('samtools fixmate');
    
    $self->switches([]);
    $self->params([]);
    
    $self->register_output_file_to_check($out_bam) unless $out_bam eq '-';
    
    return $self->run($in_bam, $out_bam);
}

=head2 rmdup

 Title   : rmdup
 Usage   : $wrapper->rmdup('sorted_in.bam', 'out.bam');
 Function: rmdup...
 Returns : n/a
 Args    : paths to input and output bams

=cut

sub rmdup {
    my ($self, $in_bam, $out_bam) = @_;
    
    $self->exe('samtools rmdup');
    
    $self->switches([]);
    $self->params([]);
    
    $self->register_output_file_to_check($out_bam);
    
    return $self->run($in_bam, $out_bam);
}

=head2 rmdupse

 Title   : rmdupse
 Usage   : $wrapper->rmdupse('sorted_in.bam', 'out.bam');
 Function: rmdupse...
 Returns : n/a
 Args    : paths to input and output bams

=cut

sub rmdupse {
    my ($self, $in_bam, $out_bam) = @_;
    
    $self->exe('samtools rmdupse');
    
    $self->switches([]);
    $self->params([]);
    
    $self->register_output_file_to_check($out_bam);
    
    return $self->run($in_bam, $out_bam);
}

=head2 flagstat

 Title   : flagstat
 Usage   : $wrapper->flagstat('in.bam', 'out.bam.flagstat');
 Function: flagstat...
 Returns : n/a
 Args    : paths to input and output bams

=cut

sub flagstat {
    my ($self, $in_bam, $out_bam) = @_;
    
    $self->exe('samtools flagstat');
    
    $self->switches([]);
    $self->params([]);
    
    $self->register_output_file_to_check($out_bam);
    
    return $self->run($in_bam, ' > '.$out_bam);
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
    my $temp_dir = $io->tempdir();
    my $name_sorted_file = $io->catfile($temp_dir, 'namesorted');
    
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
        system("mv $tmp_bam.bam $out_bam.bam");
    }
    else {
        $self->warn("$tmp_bam.bam is bad (only $bam_count lines vs $sam_count), will unlink it");
        $self->run_status(-1);
        unlink("$tmp_bam.bam");
    }
    
    $self->run_method($orig_run_method);
    return;
}

=head2 run

 Title   : run
 Usage   : Do not call directly: use one of the other methods like index()
           instead.
 Function: Run your chosen samtools command on the supplied file(s).
 Returns : n/a
 Args    : paths to input/output files

=cut

sub _pre_run {
    my ($self, @files) = @_;
    
    #*** what actually happens if open is used??...
    #my $run_method = $self->run_method();
    #if ($run_method eq 'open') {
    #    $self->warn("the open run method isn't compatible with this wrapper, switching to system");
    #    $self->run_method('system');
    #}
    
    return @files;
}

1;
