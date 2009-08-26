=head1 NAME

VertRes::Pipelines::Release - pipeline for creating a release from mapped lanes

=head1 SYNOPSIS

# make a file of absolute paths to the lane directories you want in the release,
# eg:
find $G1K/mapping -type d -name "*RR*" | grep -v "SOLID" > lanes.fofn

# make a config file with a single directory path pointing to where you want
# the release built, and setting the lanes option
echo '/path/to/release verbose=1;do_cleanup=1;lanes=lanes.fofn' > pipeline.config

# other options you could add to the above config line include:
# do_chr_splits=1
# do_recalibration=1

# run the pipeline:
mkdir /path/to/release
run-pipeline -c pipeline.config -v -m VertRes::Pipelines::Release

# (and make sure it keeps running by adding that last to a regular cron job)

=head1 DESCRIPTION

A module for creating a release (merging, splitting, recalibration etc.) from
lanes mapped by VertRes::Pipelines::Mapping.


0) List of lanes -> create a hierarchy & symlink lane bams
1) Lanes may contain paired and/or single ended bams
2) lib-level merge, creating seperate [ps]e_raw.bam files and doing
   rmdup/rmdupse as appropriate to create [ps]e_rmdup.bam
3) platform-level merge, combining both [ps]e_rmdup.bam into raw.bam.
   Optional broad recal on all these (with option to force repeat recals).
   Option to split by chr afterwards.
   Option to filter out lanes at this stage?
   Create release.bam symlinks. bai and md5 and bas file creation.
     *** or better, hardlinks to DCC-style filenames for the bams, bais and bas
         files, then make a single md5 file for everything, making for an easy
         aspera upload
4) sample-level merge, making raw.bam (optionally from platform recals).
   Option to split by chr. for sample-level have 'latest.bam' symlinks.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Pipelines::Release;

use strict;
use warnings;
use VertRes::Utils::Hierarchy;
use HierarchyUtilities;
use VertRes::IO;
use File::Spec;
use File::Copy;
use Cwd 'abs_path';
use LSF;

use base qw(VertRes::Pipeline);

our $actions = [{ name     => 'create_release_hierarchy',
                  action   => \&create_release_hierarchy,
                  requires => \&create_release_hierarchy_requires, 
                  provides => \&create_release_hierarchy_provides },
                { name     => 'library',
                  action   => \&library,
                  requires => \&library_requires, 
                  provides => \&library_provides },
                { name     => 'lib_rmdup',
                  action   => \&lib_rmdup,
                  requires => \&lib_rmdup_requires, 
                  provides => \&lib_rmdup_provides },
                { name     => 'platform',
                  action   => \&platform,
                  requires => \&platform_requires, 
                  provides => \&platform_provides },
                { name     => 'sample',
                  action   => \&sample,
                  requires => \&sample_requires, 
                  provides => \&sample_provides },
                { name     => 'cleanup',
                  action   => \&cleanup,
                  requires => \&cleanup_requires, 
                  provides => \&cleanup_provides }];

our %options = (sequence_index => '/nfs/sf8/G1K/misc/sequence.index',
                do_cleanup => 0,
                do_chr_splits => 0,
                do_recalibration => 0,
                bsub_opts => '');

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Pipelines::Release->new(lane => '/path/to/lane');
 Function: Create a new VertRes::Pipelines::Release object.
 Returns : VertRes::Pipelines::Release object
 Args    : lane => '/path/to/desired/release_dir' (REQUIRED)
           lanes => fofn (REQUIRED; a file listing all the mapped lanes you want
                          included in the release)
           sequence_index => '/path/to/sequence.index' (there is a G1K default)
           do_cleanup => boolean (default false: don't do the cleanup action)
           do_chr_splits => boolean (default false: don't split platform-level
                                     or individual-level bams by chr)
           do_recalibration => boolean (default false: don't recalibrate any
                                        bams)
           other optional args as per VertRes::Pipeline

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(%options, actions => $actions, @args);
    
    $self->{lane} || $self->throw("lane directory not supplied, can't continue");
    $self->{lanes} || $self->throw("lanes fofn not supplied, can't continue");
    
    $self->{io} = VertRes::IO->new;
    
    return $self;
}

=head2 create_release_hierarchy_requires

 Title   : create_release_hierarchy_requires
 Usage   : my $required_files = $obj->create_release_hierarchy_requires('/path/to/lane');
 Function: Find out what files the create_release_hierarchy action needs before
           it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub create_release_hierarchy_requires {
    my $self = shift;
    return [];
}

=head2 create_release_hierarchy_provides

 Title   : create_release_hierarchy_provides
 Usage   : my $provided_files = $obj->create_release_hierarchy_provides('/path/to/lane');
 Function: Find out what files the create_release_hierarchy action generates on
           success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub create_release_hierarchy_provides {
    my ($self, $lane_path) = @_;
    return ['.release_hierarchy_made'];
}

=head2 create_release_hierarchy

 Title   : create_release_hierarchy
 Usage   : $obj->create_release_hierarchy('/path/to/lane', 'lock_filename');
 Function: Create a release hierarchy with mapped lane-level bams symlinked in.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub create_release_hierarchy {
    my ($self, $lane_path, $action_lock) = @_;
    
    my @mapped_lanes = $self->{io}->parse_fod($self->{lanes});
    
    my $hu = VertRes::Utils::Hierarchy->new(verbose => $self->verbose);
    $hu->check_lanes_vs_sequence_index(\@mapped_lanes, $self->{sequence_index}) || $self->throw("There was inconsistency between the mapped lanes and the sequence index, can't continue");
    
    my @bams_in_rel = $hu->create_release_hierarchy(\@mapped_lanes, $lane_path);
    @bams_in_rel || $self->throw("Failed to create the release hierarchy");
    
    my $done_file = $self->{io}->catfile($lane_path, '.release_hierarchy_made');
    open(my $dfh, '>', $done_file) || $self->throw("Could not write to $done_file");
    foreach my $bam_in_rel (@bams_in_rel) {
        print $dfh $bam_in_rel, "\n";
    }
    my $expected = @bams_in_rel;
    print $dfh "# expecting $expected\n";
    close($dfh);
    
    # just check that the write was successful
    sleep(5);
    my @got = $self->{io}->parse_fofn($done_file);
    unless (@got == @bams_in_rel) {
        unlink($done_file);
        $self->throw("Write to $done_file gave inconsistent data");
    }
    
    return $self->{YES};
}

=head2 library_requires

 Title   : library_requires
 Usage   : my $required_files = $obj->library_requires('/path/to/lane');
 Function: Find out what files the library action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub library_requires {
    my $self = shift;
    return ['.release_hierarchy_made'];
}

=head2 library_provides

 Title   : library_provides
 Usage   : my $provided_files = $obj->library_provides('/path/to/lane');
 Function: Find out what files the library action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub library_provides {
    my $self = shift;
    return ['.library_merge_done'];
}

=head2 library

 Title   : library
 Usage   : $obj->library('/path/to/lane', 'lock_filename');
 Function: Merge lane-level bams to the library level.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub library {
    my ($self, $lane_path, $action_lock) = @_;
    
    $self->merge_up_one_level($lane_path,
                              $action_lock,
                              '.release_hierarchy_made',
                              undef,
                              'library_merge',
                              '.library_merge_expected');
    
    # we've only submitted to LSF, so it won't have finished; we always return
    # that we didn't complete
    return $self->{No};
}


=head2 lib_rmdup_requires

 Title   : lib_rmdup_requires
 Usage   : my $required_files = $obj->lib_rmdup_requires('/path/to/lane');
 Function: Find out what files the lib_rmdup action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub lib_rmdup_requires {
    my $self = shift;
    return ['.library_merge_done'];
}

=head2 lib_rmdup_provides

 Title   : lib_rmdup_provides
 Usage   : my $provided_files = $obj->lib_rmdup_provides('/path/to/lane');
 Function: Find out what files the lib_rmdup action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub lib_rmdup_provides {
    my $self = shift;
    return ['.lib_rmdup_done'];
}

=head2 lib_rmdup

 Title   : lib_rmdup
 Usage   : $obj->lib_rmdup('/path/to/lane', 'lock_filename');
 Function: Remove duplicates at the library level.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub lib_rmdup {
    my ($self, $lane_path, $action_lock) = @_;
    
    my $fofn = $self->{io}->catfile($lane_path, '.library_merge_done');
    my @files = $self->{io}->parse_fofn($fofn, $lane_path);
    my $verbose = $self->verbose();
    
    my @rmdup_bams;
    foreach my $merge_bam (@files) {
        $merge_bam = $self->{io}->catfile($lane_path, $merge_bam);
        my ($library, $basename) = (File::Spec->splitdir($merge_bam))[-2..-1];
        
        my $rmdup_bam = $merge_bam;
        $rmdup_bam =~ s/\.bam$/.rmdup.bam/;
        push(@rmdup_bams, $rmdup_bam);
        next if -s $rmdup_bam;
        
        my $single_ended = $basename =~ /^se_/ ? 1 : 0;
        
        LSF::run($action_lock, $lane_path, $self->{prefix}.'lib_rmdup', $self,
                 qq{perl -MVertRes::Utils::Sam -Mstrict -e "VertRes::Utils::Sam->new(verbose => $verbose)->rmdup(qq[$merge_bam], qq[$rmdup_bam], single_ended => $single_ended) || die qq[rmdup failed for $merge_bam\n];"});
    }
    
    my $out_fofn = $self->{io}->catfile($lane_path, '.lib_rmdup_expected');
    open(my $ofh, '>', $out_fofn) || $self->throw("Couldn't write to $out_fofn");
    foreach my $out_bam (@rmdup_bams) {
        print $ofh $out_bam, "\n";
    }
    my $expected = @rmdup_bams;
    print $ofh "# expecting $expected\n";
    close($ofh);
    
    return $self->{No};
}

=head2 platform_requires

 Title   : platform_requires
 Usage   : my $required_files = $obj->platform_requires('/path/to/lane');
 Function: Find out what files the platform action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub platform_requires {
    my $self = shift;
    return ['.lib_rmdup_done'];
}

=head2 platform_provides

 Title   : platform_provides
 Usage   : my $provided_files = $obj->platform_provides('/path/to/lane');
 Function: Find out what files the platform action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub platform_provides {
    my $self = shift;
    return ['.platform_merge_done'];
}

=head2 platform

 Title   : platform
 Usage   : $obj->platform('/path/to/lane', 'lock_filename');
 Function: Merge library-level bams to the platform level.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub platform {
    my ($self, $lane_path, $action_lock) = @_;
    
    $self->merge_up_one_level($lane_path,
                              $action_lock,
                              '.lib_rmdup_done',
                              'raw.bam',
                              'platform_merge',
                              '.platform_merge_expected');
    
    return $self->{No};
}

=head2 sample_requires

 Title   : sample_requires
 Usage   : my $required_files = $obj->sample_requires('/path/to/lane');
 Function: Find out what files the sample action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub sample_requires {
    my $self = shift;
    return ['.platform_merge_done'];
}

=head2 sample_provides

 Title   : sample_provides
 Usage   : my $provided_files = $obj->sample_provides('/path/to/lane');
 Function: Find out what files the sample action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub sample_provides {
    my $self = shift;
    return ['.sample_merge_done'];
}

=head2 sample

 Title   : sample
 Usage   : $obj->sample('/path/to/lane', 'lock_filename');
 Function: Merge platform-level bams to the sample level.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub sample {
    my ($self, $lane_path, $action_lock) = @_;
    
    $self->merge_up_one_level($lane_path,
                              $action_lock,
                              '.platform_merge_done',
                              'raw.bam',
                              'sample_merge',
                              '.sample_merge_expected');
    
    return $self->{No};
}

=head2 merge_up_one_level

 Title   : merge_up_one_level
 Usage   : $obj->merge_up_one_level('working_directory',
                                    'lock_file_name',
                                    'fofn',
                                    'output_basename',
                                    'job_name',
                                    'out_fofn');
 Function: Given a fofn of absolute paths to bams in multiple different
           directories but at the same depth within a hierarchy, groups them
           by their containing directory's parent and basename prefix, and
           merges the groups into bam files placed in that parent directory.

           Eg. given a fofn listing:
           /.../mapping/proj1/ind1/SLX/lib1/lane1/se_map.bam
           /.../mapping/proj1/ind1/SLX/lib1/lane2/se_map.bam
           /.../mapping/proj1/ind1/SLX/lib2/lane3/pe_map.bam
           /.../mapping/proj1/ind1/SLX/lib2/lane4/se_map.filtered.bam
           /.../mapping/proj1/ind1/SLX/lib2/lane5/pe_map.bam
           /.../mapping/proj1/ind1/SLX/lib2/lane6/se_map.bam
           with output_basename undef it will create:
           /.../mapping/proj1/ind1/SLX/lib1/se_map.bam (lane1+2 merge)
           /.../mapping/proj1/ind1/SLX/lib2/pe_map.bam (lane3+5 merge)
           /.../mapping/proj1/ind1/SLX/lib2/se_map.bam (lane4+6 merge)
           and with the output_basename set to 'merged.bam' it will create:
           /.../mapping/proj1/ind1/SLX/lib1/merged.bam (lane1+2 merge)
           /.../mapping/proj1/ind1/SLX/lib2/merged.bam (lane3,4,5,6 merge)

 Returns : n/a (creates out_fofn which will be a fofn of absolute paths to
           the merged bams it was supposed to create)
 Args    : a file basename to name output bams after

=cut

sub merge_up_one_level {
    my ($self, $lane_path, $action_lock, $fofn, $output_basename, $job_name, $out_fofn) = @_;
    $fofn = $self->{io}->catfile($lane_path, $fofn);
    $out_fofn = $self->{io}->catfile($lane_path, $out_fofn);
    
    unlink($out_fofn);
    
    my $group_by_basename = ! defined $output_basename;
    my %grouped_bams = $self->_fofn_to_bam_groups($lane_path, $fofn, $group_by_basename);
    
    my $verbose = $self->verbose;
    
    my @out_bams;
    while (my ($group, $bams) = each %grouped_bams) {
        my @bams = @{$bams};
        @bams || next;
        
        my ($out_dir, $out_prefix) = split(':!:', $group);
        my $out_bam;
        if ($out_prefix) {
            $out_bam = $out_prefix.'.bam';
        }
        else {
            $out_bam = $output_basename;
        }
        $out_bam || $self->throw("Couldn't work out what to call the merged bam given ($lane_path, $output_basename, $job_name, $group)!");
        
        $out_bam = $self->{io}->catfile($lane_path, $out_dir, $out_bam);
        push(@out_bams, $out_bam);
        
        next if -s $out_bam;
        
        LSF::run($action_lock, $lane_path, $self->{prefix}.$job_name, $self,
                 qq{perl -MVertRes::Utils::Sam -Mstrict -e "VertRes::Utils::Sam->new(verbose => $verbose)->merge(qq[$out_bam], qw(@bams)) || die qq[merge failed for (@bams) -> $out_bam\n];"});
    }
    
    open(my $ofh, '>', $out_fofn) || $self->throw("Couldn't write to $out_fofn");
    foreach my $out_bam (@out_bams) {
        print $ofh $out_bam, "\n";
    }
    my $expected = @out_bams;
    print $ofh "# expecting $expected\n";
    close($ofh);
    
    return;
}

sub _fofn_to_bam_groups {
    my ($self, $lane_path, $fofn, $group_by_basename) = @_;
    
    # /.../mapping/proj1/ind1/SLX/lib1/lane1/pe_raw.sorted.bam
    
    my @files = $self->{io}->parse_fofn($fofn, $lane_path);
    
    my %bams_groups;
    foreach my $path (@files) {
        my @parts = File::Spec->splitdir($path);
        my $basename = pop @parts;
        pop @parts;
        my $parent = File::Spec->catfile(@parts);
        
        my $group;
        if ($group_by_basename) {
            $basename =~ /^([^.]+).*\.bam$/ || $self->throw("Got a bad bamfile name in the fofn: $_");
            $group = "$parent:!:$1";
        }
        else {
            $group = $parent;
        }
        
        push(@{$bams_groups{$group}}, $self->{io}->catfile($lane_path, $path));
    }
    
    return %bams_groups;
}

=head2 cleanup_requires

 Title   : cleanup_requires
 Usage   : my $required_files = $obj->cleanup_requires('/path/to/lane');
 Function: Find out what files the cleanup action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub cleanup_requires {
    my ($self, $lane_path) = @_;
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
 Function: Unlink all the pipeline-related files (_*) in a lane, as well
           as the .sai files and split directory.
           NB: do_cleanup => 1 must have been supplied to new();
 Returns : $VertRes::Pipeline::Yes
 Args    : lane path, name of lock file to use

=cut

sub cleanup {
    my ($self, $lane_path, $action_lock) = @_;
    return $self->{Yes} unless $self->{do_cleanup};
    
    my $prefix = $self->{prefix};
    
    foreach my $file (qw(log)) {
        unlink($self->{io}->catfile($lane_path, $prefix.$file));
    }
    
    my $file_base = $self->{io}->catfile($lane_path, $prefix);
    foreach my $job_base (qw(library_merge lib_rmdup platform_merge sample_merge)) {
        system("rm $file_base$job_base.o $file_base$job_base.e");
    }
    
    return $self->{Yes};
}

sub is_finished {
    my ($self, $lane_path, $action) = @_;
    
    my $action_name = $action->{name};
    
    if ($action_name eq 'cleanup') {
        return $self->{No};
    }
    elsif ($action_name eq 'lib_rmdup') {
        my $expected_file = $self->{io}->catfile($lane_path, ".lib_rmdup_expected");
        my $done_file = $self->{io}->catfile($lane_path, ".lib_rmdup_done");
        $self->_merge_check($expected_file, $done_file);
    }
    elsif ($action_name eq 'library' || $action_name eq 'platform' || $action_name eq 'sample') {
        my $expected_file = $self->{io}->catfile($lane_path, ".${action_name}_merge_expected");
        my $done_file = $self->{io}->catfile($lane_path, ".${action_name}_merge_done");
        $self->_merge_check($expected_file, $done_file);
    }
    
    return $self->SUPER::is_finished($lane_path, $action);
}

sub _merge_check {
    my ($self, $expected_file, $done_file) = @_;
    
    if (-s $expected_file && ! -e $done_file) {
        my $done_bams = 0;
        my $expected_bams = 0;
        my $written_expected;
        open(my $efh, $expected_file) || $self->throw("Couldn't open $expected_file");
        my @bams;
        while (<$efh>) {
            chomp;
            /\S/ || next;
            if (/^# expecting (\d+)/) {
                $written_expected = $1;
                next;
            }
            $expected_bams++;
            $done_bams++ if -s $_;
        }
        
        if ($written_expected == $expected_bams && $done_bams == $expected_bams) {
            move($expected_file, $done_file) || $self->throw("Failed to move $expected_file -> $done_file");
        }
    }
}

1;

