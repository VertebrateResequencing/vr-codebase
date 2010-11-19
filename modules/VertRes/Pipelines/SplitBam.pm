=head1 NAME

VertRes::Pipelines::SplitBam - pipeline for splitting bams into chromosome files

=head1 SYNOPSIS

# *** To do: take care of 454 bams made with ssaha2 and our cigar->sam code.
#
# Make a file of bam filenames which you'd like split into chromosomes.
# The output for each bam will be one bam per chromosome.
# /path/to/foo.bam will have output files put in directory /path/to/foo.split
#
# Make a conf file with...
# Optional settings also go here.
#
# Example splitbam.conf:
root    => '/abs/path/to/output/dir',  
module  => 'VertRes::Pipelines::SplitBam',
prefix  => '_',
bams_fofn => 'bams.fofn'
simultaneous_splits => 5,
data => {
    pretend => 0,
    split_bam_by_sequence_opts => {
        only => '^[0-9]',
        check => 1},
    }
}
#
# 'root' can be anything, just there to keep run-pipeline happy
#
# Option which can go in the data {} section are:

# split_bam_by_sequence_opts => hash of options to be passed into
#   VertRes::Utils::Sam::split_bam_by_sequence
#    For possible options, see the usage of this subroutine.
#    Exceptions:
#     1) output_dir will be ignored; the output directory
#        corresponding to each bam /path/to/foo.bam to be split is
#        /path/to/foo.split
#     2) pretend ignored; if you want this, use pretend option
#        in the data section, not here.
#
# pretend => boolean (default false; if set to true, will print the
#                     bam files which will be made, but not make them)
#
# make a pipeline file:
echo "__SPLITBAM__ splitbam.conf" > splitbam.pipeline

# run the pipeline:
run-pipeline -c splitbam.pipeline -v

# (and make sure it keeps running by adding that last to a regular cron job)

=head1 DESCRIPTION

A module for splitting bam files into chromosome bams.  Essentially calls
VertRes::Utils::split_bam_by_sequence on each input bam file.
Input bam files specified by a file of bam filenames.

=head1 AUTHOR

Martin Hunt: mh12@sanger.ac.uk

=cut

package VertRes::Pipelines::SplitBam;

use strict;
use warnings;
use VertRes::IO;
use VertRes::Utils::FileSystem;
use File::Basename;
use File::Spec;
use File::Copy;
use Cwd 'abs_path';
use LSF;

use base qw(VertRes::Pipeline);

our $actions = [{ name     => 'split',
                  action   => \&split,
                  requires => \&split_requires,
                  provides => \&split_provides } ];

our %options = (simultaneous_splits => 5,
                bsub_opts => '',);


=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Pipelines::SplitBam->new(lane => '/path/to/lane');
 Function: Create a new VertRes::Pipelines::SplitBam object.
 Returns : VertRes::Pipelines::SplitBam object
 Args    : lane_path => '/path/to/output_dir' (REQUIRED, set by
                         run-pipeline automatically)

           bams_fofn => 'bams.fofn' (REQUIRED, file of bam filenames you want split)

           simultaneous_splits => int (default 5; the number of splits to do at once -
                                     limited to avoid IO problems)

           split_bam_by_sequence_opts => hash (optional; hash of options to be passed
                                into VertRes::Utils::Sam::split_bam_by_sequence.
                                For possible options, see the usage of this subroutine.
                                output_dir will be ignored; the output directory
                                corresponding to each bam /path/to/foo.bam to be split is
                                /path/to/foo.split)

           other optional args as per VertRes::Pipeline

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(%options, actions => $actions, @args);

    $self->{lane_path} || $self->throw("lane_path (misnomer; actually split output) directory not supplied, can't continue");
    $self->{fsu} = VertRes::Utils::FileSystem->new;
    return $self;
}

=head2 split_requires

 Title   : split_requires
 Usage   : my $required_files = $obj->split_requires('/path/to/lane');
 Function: Find out what files the split action needs before
           it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub split_requires {
    my $self = shift;
    return [];
}

=head2 split_provides

 Title   : split_provides
 Usage   : my $provided_files = $obj->split_provides('/path/to/lane');
 Function: Find out what files the split action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub split_provides {
    my ($self, $lane_path) = @_;
    return ['split.done'];
}

=head2 split

 Title   : split
 Usage   : $obj->split('/path/to/lane', 'lock_filename');
 Function: Splits a bam file into bam files per sequence
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : Output directory of split bams, name of lock file

=cut

sub split {
    my ($self, $out_dir, $action_lock) = @_;
    my $perl_out = $self->{fsu}->catfile($out_dir, $self->{prefix}.'split.pl');

    # override the option output_dir to split_bam_by_sequence
    # specifying the ouptput directory.
    $self->{split_bam_by_sequence_opts}->{output_dir} = $out_dir;

    # if we're only pretending, then just print the list of files
    if ($self->{pretend}) {
        my $samtools = VertRes::Utils::Sam->new();
        $self->{split_bam_by_sequence_opts}->{pretend} = 1;
        my @expected = $samtools->split_bam_by_sequence($self->{bam}, %{$self->{split_bam_by_sequence_opts}});
        print "  Files which would be made from $self->{bam}\n";
        print join "\n\t", @expected, "\n";

        my $done_file = $self->{fsu}->catfile($out_dir, 'split.done');
        Utils::CMD("touch $done_file");
        return $self->{Yes};
    }

    # we override the pretend option to split_bam_by_sequence
    if ($self->{split_bam_by_sequence_opts}->{pretend}) {
        $self->{split_bam_by_sequence_opts}->{pretend} = 0;
        $self->warn("'pretend' option in split_bam_by_sequence_opts must be false, changing to false and continuing...\n");
    }

    # write and bsub perl script to do the split
    my $split_opts_dump = Data::Dumper->new([$self->{split_bam_by_sequence_opts}], ["opts"]);

    open my $fh, '>', $perl_out  or $self->throw("Couldn't write to $perl_out");

    print $fh qq[
use strict;
use warnings;
use VertRes::Utils::Sam;
my \$o = VertRes::Utils::Sam->new();
my ];
    print $fh $split_opts_dump->Dump;
    print $fh "my \@expected_bams = \$o->split_bam_by_sequence(qq[$self->{bam}], \%\$opts);\n";
    print $fh q[
open my $ofh, '>', '.split_expected' or die "Couldn't write to .split_epected";
foreach my $out_bam (@expected_bams) {
    print $ofh $out_bam, "\n";
}
my $expected = @expected_bams;
print $ofh "# expecting $expected\n";
close $ofh;
];

    close $fh;

    LSF::run($action_lock, $out_dir, $perl_out, $self, qq[perl $perl_out]);
    return $self->{No};
}


sub is_finished {
    my ($self, $out_dir, $action) = @_;

    my $expected_file = $self->{fsu}->catfile($out_dir, '.split_expected');
    my $done_file = $self->{fsu}->catfile($out_dir, 'split.done');

    if (-s $expected_file && ! -e $done_file) {
        my $done_bams = 0;
        my $expected_bams = 0;
        my $written_expected;
        open(my $efh, $expected_file) || $self->throw("Couldn't open $expected_file");
        my @bams;
        my @skipped;
        while (<$efh>) {
            chomp;
            /\S/ || next;
            if (/^# expecting (\d+)/) {
                $written_expected = $1;
                next;
            }
            $expected_bams++;
            if ($self->{fsu}->file_exists($_)) {
                $done_bams++;
                push(@bams, $_);
            }
            else {
                push(@skipped, $_);
            }
        }
        
        if ($written_expected == $expected_bams && $done_bams == $expected_bams) {
            move($expected_file, $done_file) || $self->throw("Failed to move $expected_file -> $done_file");
        }
    }

    return $self->SUPER::is_finished($out_dir, $action);
}

1;
