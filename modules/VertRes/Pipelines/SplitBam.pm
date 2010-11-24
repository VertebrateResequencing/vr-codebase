=head1 NAME

VertRes::Pipelines::SplitBam - pipeline to split BAMs by reference sequence (wrapper for VertRes::Utils::Sam::split_bam_by_sequence())

=head1 SYNOPSIS

# *** To do: take care of 454 bams made with ssaha2 and our cigar->sam code.
#
# Make a file of BAM filenames. Each BAM file will be split into separate BAM
# files in a manner specified in the config file.  Files can be split based
# on the sequence to which each read was mapped, and also make an unmapped
# reads file (plus anything else which
# VertRes::Utils::Sam::split_bam_by_sequence() can do)
#
# Make a config file specifying your BAM file of filenames, plus
# other optional arguments.
#
# An example config file is:
root    => '/abs/path/to/output/dir',  
module  => 'VertRes::Pipelines::SplitBam',
prefix  => '_',
bams_fofn => 'bams.fofn',
output_dir => 'split',
simultaneous_splits => 100,
data => {
    pretend => 0,
    split_bam_by_sequence_opts => {
        only => '^[0-9]'},
    }
}
#
# Options outside the data section include:
# 'root' can be anything, just there to keep run-pipeline happy
#
# output_dir => string (Optional.  If not given, the output BAMs will be
#                       in same directory as the original BAM to be split.
#                       If given, e.g. output_dir => 'split', then the output
#                       directory for /foo/in.bam will be /foo/split)
#
# Options which can go in the data {} section are:
#
# pretend => boolean (default false; if set to true, will print the
#                     BAM files which will be made, but not make them)
#
# split_bam_by_sequence_opts => hash of options to be passed into
#    VertRes::Utils::Sam::split_bam_by_sequence
#    For possible options, see the usage of this subroutine.
#    Exceptions:
#     1) output_dir will be ignored; the output directory
#        corresponding to each BAM to be split is
#        specified by 'output_dir' outside the data section.
#     2) pretend ignored; if you want this, use pretend option
#        in the data section, not here.
#
#
# make a pipeline file:
echo "__SPLITBAM__ splitbam.conf" > splitbam.pipeline

# run the pipeline:
run-pipeline -c splitbam.pipeline -v

# (and make sure it keeps running by adding that last to a regular cron job)

=head1 DESCRIPTION

A module for splitting BAM files into BAMs by reference sequence name,
and also unmapped reads.  e.g. split one BAM into 1 new BAM for
each reference sequence, plus a file of unmapped reads.
Essentially calls VertRes::Utils::split_bam_by_sequence on each
input BAM file. Input BAM files specified by a file of BAM filenames.

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
                  provides => \&split_provides },
                { name     => 'cleanup',
                  action   => \&cleanup,
                  requires => \&cleanup_requires,
                  provides => \&cleanup_provides }] ;

our %options = (simultaneous_splits => 100,
                bsub_opts => '',);


=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Pipelines::SplitBam->new(lane => '/path/to/lane');
 Function: Create a new VertRes::Pipelines::SplitBam object.
 Returns : VertRes::Pipelines::SplitBam object
 Args    : lane_path => '/path/to/output_dir' (REQUIRED, set by
                         run-pipeline automatically)

           bams_fofn => 'bams.fofn' (REQUIRED, file of BAM filenames you want split)

           output_dir => 'outdir' (Optional.  For a BAM file /path/in.bam to be
                                 split, will put the split BAMs in /path/outdir/.
                                 If this option is not used, output BAMs in same
                                 directory as input BAM, e.g. /path/ for /path/in.bam)


           simultaneous_splits => int (default 100; the number of splits to do at once -
                                     limited to avoid IO problems)

           split_bam_by_sequence_opts => hash (optional; hash of options to be passed
                                into VertRes::Utils::Sam::split_bam_by_sequence.
                                For possible options, see the usage of this subroutine.
                                output_dir will be ignored; the output directory
                                corresponding to each BAM is specified by outpt_dir
                                (see above))

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
 Function: Splits a BAM file into BAM files per sequence
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
        foreach my $file (@expected){
            my ($basename, $path) = fileparse($file);
            my $up_dir = dirname($path);
            my $outfile = File::Spec->catfile($up_dir, $basename);
            print "    $outfile\n";
        }

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

    print $fh qq[ use strict;
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

=head2 cleanup_requires

 Title   : cleanup_requires
 Usage   : my $required_files = $obj->cleanup_requires('/path/to/lane');
 Function: Find out what files the cleanup action needs before
           it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub cleanup_requires {
    my $self = shift;
    return [];
}

=head2 cleanup_provides

 Title   : cleanup_provides
 Usage   : my $provided_files = $obj->cleanup_provides('/path/to/lane');
 Function: Find out what files the cleanup action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub cleanup_provides {
    my ($self, $lane_path) = @_;
    return [];
}

=head2 cleanup

 Title   : cleanup
 Usage   : $obj->cleanup('/path/to/lane', 'lock_filename');
 Function: Unlink all the pipeline-related files (_*) in a lane.
           NB: do_cleanup => 1 must have been supplied to new();
 Returns : $VertRes::Pipeline::Yes
 Args    : lane path, name of lock file to use

=cut

sub cleanup {
    my ($self, $lane_path, $action_lock) = @_;
    return $self->{Yes} unless $self->{do_cleanup};
    
    my $prefix = $self->{prefix};

    foreach my $file (qw(log job_status)) {
        unlink($self->{fsu}->catfile($lane_path, $prefix.$file));
    }

    my $file_base = $self->{fsu}->catfile($lane_path, $prefix);

    foreach my $job_base (qw(split.pl)) {
        unlink "$file_base$job_base";
        foreach my $suffix ('o', 'e') {
            unlink("$file_base$job_base.$suffix");
        }
    }

   return $self->{Yes};
}

sub is_finished {
    my ($self, $out_dir, $action) = @_;
    my $action_name = $action->{name};

    if ($action_name eq 'cleanup') {
        return $self->{No};
    }

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
            foreach my $bam (@bams){
                my ($bam_basename, $path) = fileparse($bam);
                my $up_dir = dirname($path);
                move($bam, $up_dir) || $self->throw("Failed to move $bam up one directory");
            }
            move($expected_file, $done_file) || $self->throw("Failed to move $expected_file -> $done_file");
        }
    }

    return $self->SUPER::is_finished($out_dir, $action);
}

1;
