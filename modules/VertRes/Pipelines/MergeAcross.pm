=head1 NAME

VertRes::Pipelines::MergeAcross - pipeline for merging groups of bam files

=head1 SYNOPSIS

# Make one file of filenames for each group of bam files you'd like merged.
# The BAMs in a given fofn will be merged either into one file, or 
# merged in other ways depending on the options chosen.
#
# Make a conf file with root pointing to where you'd like the merged bams and
# that specifies how to group the bam files.
# Optional settings also go here.
#
# Example mergeAcross.conf:
root    => '/abs/path/to/output/dir',
module  => 'VertRes::Pipelines::MergeAcross',
prefix  => '_',
data => {
    groups => {
        group_1 => 'group_1.fofn',
        group_2 => 'group_2.fofn',
    },
}
# The result of this would be two merged bam files, one for each group:
# /abs/path/to/output/dir/group_1.bam (merge of all files in group_1.fofn)
# /abs/path/to/output/dir/group_2.bam (merge of all files in group_2.fofn)
#
# In addition to merging above, files within a group can be subgrouped
# according to their basename using
common_prefix_merge => 1
# in the data section.  By default everything up to the first dot in the
# basename is used to determine the groups (can be changed with the 
# option common_prefix_regex).  For example, with
# groups => {...} as above, if group_1.fofn had these files:
#
# /path_a/chrom1.bam
# /path_a/chrom2.bam
# /path_b/chrom1.bam
# /path_b/chrom2.bam
#
# and group_2.fofn had these:
#
# /path_x/chrom1.bam
# /path_x/chrom2.bam
# /path_y/chrom1.bam
# /path_y/chrom2.bam
#
# then the resulting output files would be:
# /abs/path/to/output/dir/group_1/chrom1.bam  (merge of /path_{a,b}/chrom1.bam)
# /abs/path/to/output/dir/group_1/chrom2.bam  (merge of /path_{a,b}/chrom2.bam)
# /abs/path/to/output/dir/group_2/chrom1.bam  (merge of /path_{x,y}/chrom1.bam)
# /abs/path/to/output/dir/group_2/chrom2.bam  (merge of /path_{x,y}/chrom2.bam)
#
# Other options which can go in the data {} section are:
# max_merges => int (default 10, the number of simultanous merges
#                    to be run.  Limited to avoid IO problems)
# run_index => bool (default false, run samtools index on the merged bams)
#
# make a pipeline file:
echo "/abs/path/to/output/dir mergeAcross.conf" > mergeAcross.pipeline

# run the pipeline:
run-pipeline -c mergeAcross.pipeline -v

# (and make sure it keeps running by adding that last to a regular cron job)

=head1 DESCRIPTION

A module for merging groups of BAM files into one file per group.

=head1 AUTHOR

Martin Hunt: mh12@sanger.ac.uk

=cut

package VertRes::Pipelines::MergeAcross;

use strict;
use warnings;
use VertRes::IO;
use VertRes::Utils::FileSystem;
use File::Basename;
use File::Copy;
use Cwd 'abs_path';
use VertRes::LSF;
use Data::Dumper;
use Utils;

use base qw(VertRes::Pipeline);

our $actions = [{ name     => 'merge',
                  action   => \&merge,
                  requires => \&merge_requires,
                  provides => \&merge_provides },
                { name     => 'cleanup',
                  action   => \&cleanup,
                  requires => \&cleanup_requires,
                  provides => \&cleanup_provides } ];

our %options = (max_merges => 10,
                bsub_opts => '',
                run_index => 0,
                do_cleanup => 0);


=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Pipelines::MergeAcross->new(lane => '/path/to/lane');
 Function: Create a new VertRes::Pipelines::MergeAcross object.
 Returns : VertRes::Pipelines::MergeAcross object
 Args    : lane_path => '/path/to/dindel_group_dir' (REQUIRED, set by
                         run-pipeline automatically)

           groups => {group1 => 'group1.fofn', ...}  (REQUIRED; specify
                      the BAMs which wil be grouped together to make merged
                      BAMs.  Grouping depends on common_prefix_merge)
  
           common_prefix_merge => bool (default false; if true, then BAMs which
                                        are in the same group.fofn file
                                        and share a common prefix in their
                                        basename will be merged together.
                                        Definiton of prefix can be changed with
                                        common_prefix_regex option.
                                        See synopsis for an example.)

           common_prefix_regex => 'regular expression' (default '^(.*?)\.', which
                                        gets everything before the first . in the name.
                                        BAMs will be put in the same group when the match
                                        to this regular expression ($1) is the same.
                                        See synopsis for an example with the default.)

           max_merges => int (default 10; the number of merges to do at once -
                                     limited to avoid IO problems)

           run_index => bool (default false; index the merged bam file) 

           memory => int (default 500; Memory in MB to reserve for merge jobs.)
           queue => string (default 'normal'; bsub queue to use for merge jobs.  This
                            is forced to 'hugemem' if memory > 16000)

           do_cleanup => bool (cleanup pipeline related files upon completion; 
                               default false)

           other optional args as per VertRes::Pipeline

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(%options, actions => $actions, @args);

    $self->{lane_path} || $self->throw("lane_path (misnomer; actually mergeAcross group output) directory not supplied, can't continue");
    $self->{groups} || $self->throw("groups hash not supplied in data => {} hash, can't continue");
    $self->{queue} = 'normal' unless $self->{queue};
    $self->{memory} = 500 unless $self->{memory};
    $self->{common_prefix_regex} = '^(.*?)\.' unless $self->{common_prefix_regex};
    $self->{io} = VertRes::IO->new;
    $self->{fsu} = VertRes::Utils::FileSystem->new;
    
    my %bam_groups; # will be filled with group output directory => array of filenames

    # work out the bam grouping
    while (my ($group, $fofn) = each(%{$self->{groups}})) {
        my @files = $self->{io}->parse_fofn($fofn, "/");
        my $group_key = $group;

        foreach my $bam (@files) {
            if ($self->{common_prefix_merge}) {
                my ($bam_base, $bam_path) = fileparse($bam);

                if ($bam_base =~ m/$self->{common_prefix_regex}/){
                    $group_key = $self->{fsu}->catfile($group, $1);
                }
                else {
                    $self->throw("No match to regular expression '$self->{common_prefix_regex}' in bam file $bam");
                }
            } 

            $bam_groups{$group_key} = () unless defined $bam_groups{$group_key};
            push @{$bam_groups{$group_key}}, $bam;
        }
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
    my $memory = $self->{memory}; 
    my $java_mem = int($memory * 0.9);
    my $queue = $memory >= 16000 ? "hugemem" : $self->{queue};
    my $orig_bsub_opts = $self->{bsub_opts};
    $self->{bsub_opts} = "-q $queue -M${memory} -R 'select[mem>$memory] rusage[mem=$memory]'";

    # for each group of bam files, run the merge if not done already
    # (subject to max_jobs constraint)
    while (my ($group, $bams) = each(%{$self->{bams_by_group}})) {
        my $bam_out = $self->{fsu}->catfile($work_dir, "$group.bam");
        my ($bam_out_base, $bam_out_dir) = fileparse($bam_out);
        my $prefix = $bam_out_base;
        $prefix =~ s/\.bam$//;
        my $job_name = $self->{prefix} . $prefix. '.pl';
        $prefix = $self->{fsu}->catfile($bam_out_dir, "$self->{prefix}$prefix");
        my $tmp_bam_out = "$prefix.tmp.bam";
        my $bai = "$bam_out.bai";
        my $tmp_bai = "$tmp_bam_out.bai";
        my $jids_file = "$prefix.jids";
        my $perl_out = "$prefix.pl";
        my $status =  VertRes::LSF::is_job_running($jids_file);
        
        # Check that all BAMs are present before proceeding
        if ($status & $VertRes::LSF::No) {
            my $ready = 1;
            foreach my $bam (@{$bams}) {
                unless ($$self{fsu}->file_exists($bam)) {
                    $ready = 0;
                    $self->debug(qq[All BAMs for group $group are not present... skipping]);
                    last;
                }
            }
            next unless $ready;
        }

        unless (-d $bam_out_dir) {
            File::Path::make_path($bam_out_dir) || $self->throw("Error making directory '$bam_out_dir'");
        }

        if ($status & $VertRes::LSF::Done and $$self{fsu}->file_exists($bam_out)) {
            ++$jobs_done;
            next;
        }
        elsif ($status & $VertRes::LSF::Running) {
            $jobs_running++;
            next;
        } 
        if ($status & $VertRes::LSF::Error) { 
            $self->warn("The command failed: $perl_out\n");
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
my \$o = VertRes::Wrapper::picard->new(java_memory => $java_mem);
my ];
        print $fh $d->Dump;
        print $fh qq[
\$o->merge_and_check('$tmp_bam_out', \$bams);
];

        if ($self->{run_index}) {
            print $fh qq[
\$o = VertRes::Wrapper::samtools->new();
\$o->index('$tmp_bam_out', '$tmp_bai');
rename '$tmp_bai', '$bai';
];
        }
        print $fh qq[
rename '$tmp_bam_out', '$bam_out';
];
        close $fh;

        $self->archive_bsub_files($work_dir, "$self->{prefix}$group.pl");
        VertRes::LSF::run($jids_file, $bam_out_dir, $job_name, $self, "perl -w $perl_out");
        print STDERR "    Submitted $perl_out\n"
    }

    $self->{bsub_opts} = $orig_bsub_opts;
    
    return $self->{No};
}

=head2 cleanup_requires

 Title   : cleanup_requires
 Usage   : my $required_files = $obj->cleanup_requires('/path/to/lane');
 Function: Find out what files the cleanup action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub cleanup_requires {
    my $self = shift;
    return ['merge.done'];
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
    my ($self, $work_dir, $action_lock) = @_;
    return $self->{Yes} unless $self->{do_cleanup};
    
    my $prefix = $self->{prefix};
    my $fsu = $self->{fsu};
    
    foreach my $file (qw(log job_status)) {
        $self->clean_file($fsu->catfile($work_dir, $prefix.$file));
    }

    
    foreach my $group (keys %{$self->{bams_by_group}}) {
        my $grp_base = $fsu->catfile($work_dir, $group);
        my ($base, $dir) = fileparse($grp_base);
        $base = $fsu->catfile($dir, $prefix.$base);
        
        foreach my $suffix (qw(.pl.o .pl.e .jids .pl)) {
            $self->clean_file($base.$suffix);
        }
    }
    return $self->{Yes};
}

sub clean_file
{
    my ($self, $file) = @_;
    if ( -e $file ) { unlink($file); }
    $self->{fsu}->file_exists($file, wipe_out => 1);
}



sub is_finished {
    my ($self, $work_dir, $action) = @_;
    
    # all merges successfully completed?
    if ($action->{name} eq 'merge') {
        my $jobs_done = 0;
        foreach my $group (keys %{$self->{bams_by_group}}) {
            my $bam = $self->{fsu}->catfile($work_dir, "$group.bam");
            $jobs_done++ if $self->{fsu}->file_exists($bam);
        }
        if ( $jobs_done == scalar keys %{$self->{bams_by_group}} ) {
            Utils::CMD("touch $work_dir/merge.done");
        }
    }
    elsif ($action->{name} eq 'cleanup') {
        return $self->{No};
    }
    
    return $self->SUPER::is_finished($work_dir, $action);
}


1;
