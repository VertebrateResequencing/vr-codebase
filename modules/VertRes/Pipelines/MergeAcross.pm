=head1 NAME

VertRes::Pipelines::MergeAcross - pipeline for merging groups of bam files

=head1 SYNOPSIS

# Make one file of filenames for each group of bam files you'd like merged into
# one file.  All the files named in one fofn will be merged into one bam file.
# e.g. if you want to merge all your bams into one new bam file, then just
# give one fofn containing all the bam filenames.
#
# Make a conf file with root pointing to where you'd like the merged bams, and
# that specifies the group => fofn mapping.
# Optional settings also go here.
#
# Example mergeAcross.conf:
root    => '/abs/path/to/output/dir',
module  => 'VertRes::Pipelines::MergeAcross',
prefix  => '_',
data => {
    group_merge => {
        group_1 => 'group_1.fofn',
        group_2 => 'group_2.fofn',
        group_3 => 'group_3.fofn',
    },
}
# The result of this would be three merged bam files, one for each group.
# /abs/path/to/output/dir/group_{1,2,3}.bam
#
# Merges can be done explicitly as above, or using one of these options
# instead of group_merge:
# regex_merge => 'bams.fofn',
# population_merge => 'bams.fofn', (not yet implemented)
#
# Other options which can go in the data {} section are:
# max_merges => int (default 5, the number of simultanous merges
#                    to be run.  Limited to avoid IO problems)
# run_index => bool (default false, run samtools index on the merged bams)
#
# make a pipeline file:
echo "/abs/path/to/output/dir mergeAcross.conf" > mergeAcross.pipeline

# run the pipeline:
run-pipeline -c mergeAcross.pipeline -v

# (and make sure it keeps running by adding that last to a regular cron job)

=head1 DESCRIPTION

A module for merging groups of bam files into one file.  Each group is
specified by a file of bam filenames or by chromosome name in the
bam files.

=head1 AUTHOR

Martin Hunt: mh12@sanger.ac.uk

=cut

package VertRes::Pipelines::MergeAcross;

use strict;
use warnings;
use VertRes::IO;
use VertRes::Utils::FileSystem;
use VertRes::Parser::sam;
use File::Basename;
use File::Spec;
use File::Copy;
use Cwd 'abs_path';
use LSF;
use Data::Dumper;

use base qw(VertRes::Pipeline);

our $actions = [{ name     => 'merge',
                  action   => \&merge,
                  requires => \&merge_requires,
                  provides => \&merge_provides } ];

our %options = (max_merges => 5,
                bsub_opts => '',
                run_index => 0);


=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Pipelines::MergeAcross->new(lane => '/path/to/lane');
 Function: Create a new VertRes::Pipelines::MergeAcross object.
 Returns : VertRes::Pipelines::MergeAcross object
 Args    : lane_path => '/path/to/dindel_group_dir' (REQUIRED, set by
                         run-pipeline automatically)

           REQUIRED: one of groups or regex_merge (but not both)

           groups => {group1 => 'group1.fofn', ...}  (specify
                      the bams which wil be grouped together to make merged
                      bam.  Result is one merged bam file per fofn)
                               
           regex_merge => 'bams.fofn' (file of bam file names, bams will be
                      merged by matches in basename of file,
                      using regex_merge_regex option to define the groups.) 
          
           regex_merge_regex => 'regular expresion' (default '^chrom(.*?)\.';
                                 files are grouped by those that match the regex.
                                 Match is looked for in basename of each the bams.
                                 Default works for bams made by splitBam pipeline)

           max_merges => int (default 5; the number of merges to do at once -
                                     limited to avoid IO problems)

           run_index => bool (default false; index the merged bam file) 

           other optional args as per VertRes::Pipeline

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(%options, actions => $actions, @args);

    $self->{lane_path} || $self->throw("lane_path (misnomer; actually dindel group output) directory not supplied, can't continue");
    $self->{io} = VertRes::IO->new;
    $self->{fsu} = VertRes::Utils::FileSystem->new;
    
    my %bam_groups; # will be filled with group_name => array of filenames

    # work out the bam grouping
    if ($self->{groups}){
        while (my ($group, $fofn) = each(%{$self->{groups}})) {
            my @files = $self->{io}->parse_fofn($fofn, "/");
            $bam_groups{$group} = \@files;
        }
    }
    elsif ($self->{regex_merge}){
        $self->{regex_merge_regex} or $self->{regex_merge_regex} = '^chrom(.*?)\.';
        my @files = $self->{io}->parse_fofn($self->{regex_merge}, "/");
        foreach my $file (@files) {
            my ($basename, $path) = fileparse($file);
            if ($basename =~ m/$self->{regex_merge_regex}/) {
              unless ($bam_groups{$1}) {
                  $bam_groups{$1} = ();
              }
              push @{$bam_groups{$1}}, $file;
            }
            else {
                $self->throw("No match in regular expression $self->{regex_merge_regex} to bam file $file");
            }
        }
    }
    else {
        $self->throw("One of groups, regex_merge must be supplied");
    }

    $self->{bams_by_group} = \%bam_groups;
    return $self;
}

=head2 merge_requires

 Title   : merge_requires
 Usage   : my $required_files = $obj->merge_requires('/path/to/lane');
 Function: Find out what files the merge action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub merge_requires {
    my $self = shift;
    return [];
}

=head2 merge_provides

 Title   : merge_provides
 Usage   : my $provided_files = $obj->merge_provides('/path/to/lane');
 Function: Find out what files the merge action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub merge_provides {
    my ($self, $lane_path) = @_;
    return ['merge.done'];
}

=head2 merge

 Title   : merge
 Usage   : $obj->merge('/path/to/lane', 'lock_filename');
 Function: Merge results and output final result file.
 Returns : Nothing, writes merge.done when all merges done
 Args    : lane path, name of lock file to use

=cut

sub merge {
    my ($self, $work_dir, $action_lock) = @_;
    my $jobs_running = 0;
    my $jobs_done = 0;

    # for each group of bam files, run the merge if not done already
    # (subject to max_jobs constraint)
    while (my ($group, $bams) = each(%{$self->{bams_by_group}})) {
        my $bam_out = File::Spec->catfile($work_dir, "$group.bam");
        my $tmp_bam_out = File::Spec->catfile($work_dir, "$self->{prefix}tmp.$group.bam");
        my $bai = "$bam_out.bai";
        my $tmp_bai = "$tmp_bam_out.bai";
        my $jids_file = File::Spec->catfile($work_dir, "$self->{prefix}$group.jids");
        my $perl_out = File::Spec->catfile($work_dir, "$self->{prefix}$group.pl");
        my $status = LSF::is_job_running($jids_file);

        if ($status & $LSF::Error) { 
            $self->warn("The command failed: $perl_out\n");
        }
        elsif ($status & $LSF::Running) {
            next;
        } 
        elsif ($status & $LSF::Done and $$self{fsu}->file_exists($bam_out)) {
            $jobs_done++;
            next;
        }
            
        if ($$self{max_merges} and $jobs_running >= $$self{max_merges}) {
            $self->debug("Max job limit reached, $$self{max_merges} jobs are running.");
            last;
        }

        $jobs_running++;

        # if here, then need to make and bsub a perl script to run the merging
        # command on the current group of bams
        my $d = Data::Dumper->new([$bams], ["bams"]);
        open my $fh, ">", $perl_out or $self->throw("$perl_out: $!");
        print $fh qq[use VertRes::Wrapper::picard;
use VertRes::Wrapper::samtools;
my \$o = VertRes::Wrapper::picard->new();
my ];
        print $fh $d->Dump;
        print $fh qq[
\$o->merge_and_check('$tmp_bam_out', \$bams);
\$o = VertRes::Wrapper::samtools->new();
];

        if ($self->{run_index}) {
            print $fh qq[
\$o->index('$tmp_bam_out', '$tmp_bai');
rename '$tmp_bai', '$bai';
];
        }
        print $fh qq[
rename '$tmp_bam_out', '$bam_out';
];
        close $fh;

        $self->archive_bsub_files($work_dir, "$self->{prefix}$group.pl");
        LSF::run($jids_file, $work_dir, $perl_out, $self, "perl -w $perl_out");
        print STDERR "    Submitted $perl_out\n"
    }

    # all jobs successfully completed?
    if ($jobs_done == scalar keys %{$self->{bams_by_group}}) {
        Utils::CMD("touch $work_dir/merge.done");
    }
}

1;
