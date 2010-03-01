=head1 NAME

VertRes::Pipelines::Release - pipeline for creating a release from mapped lanes

=head1 SYNOPSIS

# make a file of absolute paths to the lane directories you want in the release,
# eg:
release_lanes.pl --root /abs/META --db g1k_meta --all_samples
                 --ignore_platform SOLID > lanes.fofn

# make a conf file with root pointing to where you want the release built, and
# setting the lanes, mapper and assembly_name options. Optional settings also go
# here. Your release.conf file may look like:
root    => '/path/to/rel_dir',
module  => 'VertRes::Pipelines::Release',
prefix  => '_',
data => {
    lanes => 'lanes.fofn',
    slx_mapper => 'bwa',
    '454_mapper' => 'ssaha',
    assembly_name => 'NCBI37',
    release_date => '20100208',
    simultaneous_merges => 200,
    
    dcc_hardlinks => 1,
    do_cleanup => 1,
    
    db => {
        database => 'g1k_meta',
        host     => 'mcs4a',
        port     => 3306,
        user     => 'vreseq_ro'
    }
}

# note that slx_mapper, 454_mapper and assembly_name options must be the same
# as those you used in the VertRes::Pipelines::Mapping pipeline config file
# used to do the mapping you want to now make a release from

# other options you could add to the release.conf include:
# do_chr_splits => 1
# do_sample_merge => 1
# skip_fails => 1
# previous_release_root => '/path/to/rel'

# the previous_release_root will result in the pipline checking if a particular
# bam needs to be (re)made depending on if lanes have been removed/added/altered
# since that previous release. If a bam doesn't need to be made, the new
# release directory will contain a symlink to the bam in the old directory

# make another config file that simply mentions the release directory and the
# previous config file:
echo "/path/to/rel_dir release.conf" > rel.pipeline

# make your release directory if it doesn't already exist:
mkdir /path/to/rel_dir

# run the pipeline:
run-pipeline -c release.pipeline -v

# (and make sure it keeps running by adding that last to a regular cron job)

=head1 DESCRIPTION

A module for creating a release (merging, splitting etc.) from
lanes mapped by VertRes::Pipelines::Mapping.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Pipelines::Release;

use strict;
use warnings;
use VertRes::Utils::Hierarchy;
use VertRes::IO;
use VertRes::Utils::Sam;
use VertRes::Utils::FileSystem;
use VRTrack::VRTrack;
use VRTrack::Lane;
use File::Basename;
use File::Spec;
use File::Copy;
use Cwd 'abs_path';
use Date::Parse;
use Time::Format;
use LSF;

use base qw(VertRes::Pipeline);

our $actions = [{ name     => 'create_release_hierarchy',
                  action   => \&create_release_hierarchy,
                  requires => \&create_release_hierarchy_requires, 
                  provides => \&create_release_hierarchy_provides },
                { name     => 'library_merge',
                  action   => \&library_merge,
                  requires => \&library_merge_requires, 
                  provides => \&library_merge_provides },
                { name     => 'lib_markdup',
                  action   => \&lib_markdup,
                  requires => \&lib_markdup_requires, 
                  provides => \&lib_markdup_provides },
                { name     => 'platform_merge',
                  action   => \&platform_merge,
                  requires => \&platform_merge_requires, 
                  provides => \&platform_merge_provides },
                { name     => 'create_release_files',
                  action   => \&create_release_files,
                  requires => \&create_release_files_requires, 
                  provides => \&create_release_files_provides},
                { name     => 'sample_merge',
                  action   => \&sample_merge,
                  requires => \&sample_merge_requires, 
                  provides => \&sample_merge_provides },
                { name     => 'cleanup',
                  action   => \&cleanup,
                  requires => \&cleanup_requires, 
                  provides => \&cleanup_provides }];

our %options = (do_cleanup => 0,
                do_chr_splits => 0,
                do_sample_merge => 0,
                simultaneous_merges => 200,
                dcc_mode => 0,
                bsub_opts => '',
                dont_wait => 1,
                previous_release_root => '');

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Pipelines::Release->new(lane => '/path/to/lane');
 Function: Create a new VertRes::Pipelines::Release object.
 Returns : VertRes::Pipelines::Release object
 Args    : lane => '/path/to/desired/release_dir' (REQUIRED)
           lanes => fofn (REQUIRED; a file listing all the mapped lanes you want
                          included in the release)
           db => {} (REQUIRED, meta database connection details)
           slx_mapper => string (REQUIRED, the mapper used for SLX lanes)
           '454_mapper' => string (REQUIRED, the mapper used for 454 lanes)
           assembly_name => string (REQUIRED, the name of the assembly)

           do_cleanup => boolean (default false: don't do the cleanup action)
           do_chr_splits => boolean (default false: don't split platform-level
                                     bams by chr)
           dcc_mode => boolean (default false; when true, implies do_chr_splits
                                and renames the per-chr bams to the DCC naming
                                convention)
           simultaneous_merges => int (default 200; the number of merge jobs to
                                       do at once - limited to avoid IO
                                       problems)
           do_sample_merge => boolean (default false: don't create sample-level
                                       bams)
           skip_fails => boolean (default false; when true, if some jobs fail,
                                  continue to the next action anyway)
           previous_release_root => string (if a part of the hierarchy didn't
                                            change since the last release,
                                            bams are symlinked instead of being
                                            remade)
           release_date => 'YYYYMMDD' (defaults to todays date; used for getting
                                       the DCC filename correct in .bas files
                                       and getting the release filenames correct
                                       - should be the same date as of the DCC
                                       sequence.index used. Not important if
                                       not doing a DCC release (dcc_mode is
                                       off); required if dcc_mode is on
           other optional args as per VertRes::Pipeline

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(%options, actions => $actions, @args);
    
    $self->{lane} || $self->throw("lane directory not supplied, can't continue");
    $self->{lanes} || $self->throw("lanes fofn not supplied, can't continue");
    $self->{db} || $self->throw("db connection details not supplied, can't continue");
    $self->{slx_mapper} || $self->throw("slx_mapper not supplied, can't continue");
    $self->{'454_mapper'} || $self->throw("454_mapper not supplied, can't continue");
    $self->{assembly_name} || $self->throw("assembly_name not supplied, can't continue");
    
    if ($self->{dcc_mode}) {
        $self->{do_chr_splits} = 1;
    }
    
    $self->{io} = VertRes::IO->new;
    $self->{fsu} = VertRes::Utils::FileSystem->new;
    
    if ($self->{dcc_mode}) {
        $self->throw("release_date must be supplied in dcc_mode") unless defined $self->{release_date};
    }
    unless (defined $self->{release_date}) {
        $self->{release_date} = "$time{'yyyymmdd'}"; 
    }
    
    $self->{vrtrack} = VRTrack::VRTrack->new($self->{db});
    
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
    # we don't supply true as the third arg to check we have all lanes, since
    # for releases we expect new lanes to be added after we start building, but
    # a release is a freeze and this is OK. release_lanes.pl, which makes our
    # list of lanes to include in the release, does this check.
    $hu->check_lanes_vs_database(\@mapped_lanes, $self->{vrtrack}) || $self->throw("There was inconsistency between the mapped lanes and the database, can't continue");
    
    my @bams_in_rel = $hu->create_release_hierarchy(\@mapped_lanes, $lane_path, %{$self});
    @bams_in_rel || $self->throw("Failed to create the release hierarchy");
    
    my $done_file = $self->{fsu}->catfile($lane_path, '.release_hierarchy_made');
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
    
    return $self->{Yes};
}

=head2 library_merge_requires

 Title   : library_merge_requires
 Usage   : my $required_files = $obj->library_merge_requires('/path/to/lane');
 Function: Find out what files the library_merge action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub library_merge_requires {
    my $self = shift;
    return ['.release_hierarchy_made'];
}

=head2 library_merge_provides

 Title   : library_merge_provides
 Usage   : my $provided_files = $obj->library_merge_provides('/path/to/lane');
 Function: Find out what files the library_merge action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub library_merge_provides {
    my $self = shift;
    return ['.library_merge_done'];
}

=head2 library_merge

 Title   : library_merge
 Usage   : $obj->library_merge('/path/to/lane', 'lock_filename');
 Function: Merge lane-level bams to the library level.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub library_merge {
    my ($self, $lane_path, $action_lock) = @_;
    
    $self->merge_up_one_level($lane_path,
                              $action_lock,
                              '.release_hierarchy_made',
                              undef,
                              'library_merge',
                              '.library_merge_expected',
                              'long');
    
    # we've only submitted to LSF, so it won't have finished; we always return
    # that we didn't complete
    return $self->{No};
}


=head2 lib_markdup_requires

 Title   : lib_markdup_requires
 Usage   : my $required_files = $obj->lib_markdup_requires('/path/to/lane');
 Function: Find out what files the lib_markdup action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub lib_markdup_requires {
    my $self = shift;
    return ['.library_merge_done'];
}

=head2 lib_markdup_provides

 Title   : lib_markdup_provides
 Usage   : my $provided_files = $obj->lib_markdup_provides('/path/to/lane');
 Function: Find out what files the lib_markdup action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub lib_markdup_provides {
    my $self = shift;
    return ['.lib_markdup_done'];
}

=head2 lib_markdup

 Title   : lib_markdup
 Usage   : $obj->lib_markdup('/path/to/lane', 'lock_filename');
 Function: Marks duplicates at the library level.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub lib_markdup {
    my ($self, $lane_path, $action_lock) = @_;
    
    my $fofn = $self->{fsu}->catfile($lane_path, '.library_merge_done');
    my @files = $self->{io}->parse_fofn($fofn, $lane_path);
    my $verbose = $self->verbose();
    
    my $orig_bsub_opts = $self->{bsub_opts};
    $self->{bsub_opts} = '-q basement -M6100000 -R \'select[mem>6100] rusage[mem=6100]\'';
    
    my @markdup_bams;
    foreach my $merge_bam (@files) {
        $merge_bam = $self->{fsu}->catfile($lane_path, $merge_bam);
        my($basename, $path) = fileparse($merge_bam);
        
        my $markdup_bam = $merge_bam;
        $markdup_bam =~ s/\.bam$/.markdup.bam/;
        push(@markdup_bams, $markdup_bam);
        next if -s $markdup_bam;
        
        my $single_ended = $basename =~ /^se_/ ? 1 : 0;
        
        my $job_name = $self->{prefix}.'lib_markdup';
        $self->archive_bsub_files($path, $job_name);
        $job_name = $self->{fsu}->catfile($path, $job_name);
        
        LSF::run($action_lock, $lane_path, $job_name, $self,
                 qq{perl -MVertRes::Utils::Sam -Mstrict -e "VertRes::Utils::Sam->new(verbose => $verbose)->markdup(qq[$merge_bam], qq[$markdup_bam]) || die qq[markdup failed for $merge_bam\n];"});
    }
    
    my $out_fofn = $self->{fsu}->catfile($lane_path, '.lib_markdup_expected');
    open(my $ofh, '>', $out_fofn) || $self->throw("Couldn't write to $out_fofn");
    foreach my $out_bam (@markdup_bams) {
        print $ofh $out_bam, "\n";
    }
    my $expected = @markdup_bams;
    print $ofh "# expecting $expected\n";
    close($ofh);
    
    $self->{bsub_opts} = $orig_bsub_opts;
    
    return $self->{No};
}

=head2 platform_merge_requires

 Title   : platform_merge_requires
 Usage   : my $required_files = $obj->platform_merge_requires('/path/to/lane');
 Function: Find out what files the platform_merge action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub platform_merge_requires {
    my $self = shift;
    return ['.lib_markdup_done'];
}

=head2 platform_merge_provides

 Title   : platform_merge_provides
 Usage   : my $provided_files = $obj->platform_merge_provides('/path/to/lane');
 Function: Find out what files the platform_merge action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub platform_merge_provides {
    my $self = shift;
    return ['.platform_merge_done'];
}

=head2 platform_merge

 Title   : platform_merge
 Usage   : $obj->platform_merge('/path/to/lane', 'lock_filename');
 Function: Merge library-level bams to the platform level.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub platform_merge {
    my ($self, $lane_path, $action_lock) = @_;
    
    $self->merge_up_one_level($lane_path,
                              $action_lock,
                              '.lib_markdup_done',
                              'raw.bam',
                              'platform_merge',
                              '.platform_merge_expected',
                              'long');
    
    return $self->{No};
}

=head2 sample_merge_requires

 Title   : sample_merge_requires
 Usage   : my $required_files = $obj->sample_merge_requires('/path/to/lane');
 Function: Find out what files the sample_merge action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub sample_merge_requires {
    my $self = shift;
    return ['.create_release_files_done'];
}

=head2 sample_merge_provides

 Title   : sample_merge_provides
 Usage   : my $provided_files = $obj->sample_merge_provides('/path/to/lane');
 Function: Find out what files the sample_merge action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub sample_merge_provides {
    my $self = shift;
    return $self->{do_sample_merge} ? ['.sample_merge_done'] : [];
}

=head2 sample_merge

 Title   : sample_merge
 Usage   : $obj->sample_merge('/path/to/lane', 'lock_filename');
 Function: Merge platform-level bams to the sample level.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub sample_merge {
    my ($self, $lane_path, $action_lock) = @_;
    return $self->{Yes} unless $self->{do_sample_merge};
    
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
                                    'out_fofn',
                                    'queue_name');
 Function: Given a fofn of absolute paths to bams in multiple different
           directories but at the same depth within a hierarchy, groups them
           by their containing directory's parent and basename prefix, and
           merges the groups into bam files placed in that parent directory.

           Eg. given a fofn listing:
           /.../mapping/proj1/ind1/SLX/lib1/lane1/se.map.bam
           /.../mapping/proj1/ind1/SLX/lib1/lane2/se.map.bam
           /.../mapping/proj1/ind1/SLX/lib2/lane3/pe.map.bam
           /.../mapping/proj1/ind1/SLX/lib2/lane4/se.map.filtered.bam
           /.../mapping/proj1/ind1/SLX/lib2/lane5/pe.map.bam
           /.../mapping/proj1/ind1/SLX/lib2/lane6/se.map.bam
           with output_basename undef it will create:
           /.../mapping/proj1/ind1/SLX/lib1/se.map.bam (lane1+2 merge)
           /.../mapping/proj1/ind1/SLX/lib2/pe.map.bam (lane3+5 merge)
           /.../mapping/proj1/ind1/SLX/lib2/se.map.bam (lane4+6 merge)
           and with the output_basename set to 'merged.bam' it will create:
           /.../mapping/proj1/ind1/SLX/lib1/merged.bam (lane1+2 merge)
           /.../mapping/proj1/ind1/SLX/lib2/merged.bam (lane3,4,5,6 merge)

 Returns : n/a (creates out_fofn which will be a fofn of absolute paths to
           the merged bams it was supposed to create)
 Args    : see usage. last option of queue name is optional, and defaults to
           'normal'

=cut

sub merge_up_one_level {
    my ($self, $lane_path, $action_lock, $fofn, $output_basename, $job_name, $out_fofn, $queue) = @_;
    $fofn = $self->{fsu}->catfile($lane_path, $fofn);
    $out_fofn = $self->{fsu}->catfile($lane_path, $out_fofn);
    
    unlink($out_fofn);
    
    my $orig_bsub_opts = $self->{bsub_opts};
    $queue ||= 'normal';
    $self->{bsub_opts} = '-q '.$queue.' -M5000000 -R \'select[mem>5000] rusage[mem=5000]\'';
    
    my $group_by_basename = ! defined $output_basename;
    my %grouped_bams = $self->_fofn_to_bam_groups($lane_path, $fofn, $group_by_basename);
    
    my $verbose = $self->verbose;
    
    my @out_bams;
    my %lane_paths;
    my $jobs = 0;
    while (my ($group, $bams) = each %grouped_bams) {
        my @bams = @{$bams};
        @bams || next;
        
        my ($out_dir, $out_prefix) = split(':!:', $group);
        my $out_bam_name;
        if ($out_prefix) {
            $out_bam_name = $out_prefix.'.bam';
        }
        else {
            $out_bam_name = $output_basename;
        }
        $out_bam_name || $self->throw("Couldn't work out what to call the merged bam given ($lane_path, $output_basename, $job_name, $group)!");
        
        my $current_path = $self->{fsu}->catfile($lane_path, $out_dir);
        my $out_bam = $self->{fsu}->catfile($current_path, $out_bam_name);
        push(@out_bams, $out_bam);
        
        # no need to merge if the hierarchy beneath hasn't changed since the
        # previous release
        if ($self->{previous_release_root}) {
            my $previous_path = $self->{fsu}->catfile($self->{previous_release_root}, $out_dir);
            my $previous_bam =  $self->{fsu}->catfile($previous_path, $out_bam_name);
            
            if ($self->{fsu}->directory_structure_same($previous_path, $current_path, leaf_mtimes => \%lane_paths)) {
                # check that none of the lanes in this part of the hierarchy
                # have been deleted and recreated since the previous release
                my $lanes_changed = 0;
                foreach my $lane_hname (keys %{$lane_paths{$previous_path}}) {
                    my $lane = VRTrack::Lane->new_by_hierarchy_name($self->{vrtrack}, $lane_hname) || $self->throw("Lane $lane_hname was in the release, but not in the db!");
                    my $changed = str2time($lane->changed);
                    
                    if ($changed > $lane_paths{$previous_path}->{$lane_hname}) {
                        $lanes_changed = 1;
                    }
                }
                
                unless ($lanes_changed) {
                    symlink($previous_bam, $out_bam);
                    
                    # symlink the other files we might have made from that bam
                    # in other actions
                    my $prev_markdup_bam = $previous_bam;
                    $prev_markdup_bam =~ s/\.bam$/.markdup.bam/;
                    if (-s $prev_markdup_bam) {
                        my $cur_markdup_bam = $out_bam;
                        $cur_markdup_bam =~ s/\.bam$/.markdup.bam/;
                        symlink($prev_markdup_bam, $cur_markdup_bam);
                    }
                    
                    foreach my $suffix2 ('', '.md5') {
                        foreach my $suffix ('', '.bai', '.bas') {
                            my $prev = $previous_bam.$suffix.$suffix2;
                            if (-s $prev) {
                                my $cur = $out_bam.$suffix.$suffix2;
                                symlink($prev, $cur);
                            }
                        }
                    }
                    
                    # we don't merge split bams; we merge the parent, and split
                    # the merge if we want to
                }
            }
        }
        
        next if -s $out_bam;
        
        # don't do more than desired merges at once, or we'll kill IO and jobs
        # will fail
        next if $jobs >= $self->{simultaneous_merges};
	
        my (undef, $path) = fileparse($out_bam);
        my $this_job_name = $self->{prefix}.$job_name;
        my $pathed_job_name = $self->{fsu}->catfile($path, $this_job_name);
        my $lock_file = $pathed_job_name.'.jids';
        
        # keep simultaneous_merges jobs running all the time
        my $is_running = LSF::is_job_running($lock_file);
        if ($is_running & $LSF::Error) {
            warn "$job_name failed!\n";
            next;
        }
        elsif ($is_running & $LSF::Running) {
            $jobs++;
            next;
        }
        elsif ($is_running & $LSF::Done) {
            next;
        }
        else {
            $jobs++;
            $self->archive_bsub_files($path, $this_job_name);
            
            LSF::run($action_lock, $lane_path, $pathed_job_name, $self,
                 qq{perl -MVertRes::Utils::Sam -Mstrict -e "VertRes::Utils::Sam->new(verbose => $verbose)->merge(qq[$out_bam], qw(@bams)) || die qq[merge failed for (@bams) -> $out_bam\n];"});
        }
    }
    
    open(my $ofh, '>', $out_fofn) || $self->throw("Couldn't write to $out_fofn");
    foreach my $out_bam (@out_bams) {
        print $ofh $out_bam, "\n";
    }
    my $expected = @out_bams;
    print $ofh "# expecting $expected\n";
    close($ofh);
    
    $self->{bsub_opts} = $orig_bsub_opts;
    
    return;
}

sub _fofn_to_bam_groups {
    my ($self, $lane_path, $fofn, $group_by_basename) = @_;
    
    # /.../mapping/proj1/ind1/SLX/lib1/lane1/pe.raw.sorted.bam
    
    my @files = $self->{io}->parse_fofn($fofn, $lane_path);
    
    my %bams_groups;
    foreach my $path (@files) {
        my @parts = File::Spec->splitdir($path);
        my $basename = pop @parts;
        pop @parts;
        my $parent = $self->{fsu}->catfile(@parts);
        
        my $group;
        if ($group_by_basename) {
            $basename =~ /^([^.]+).*\.bam$/ || $self->throw("Got a bad bamfile name in the fofn: $_");
            $group = "$parent:!:$1";
        }
        else {
            $group = $parent;
        }
        
        push(@{$bams_groups{$group}}, $self->{fsu}->catfile($lane_path, $path));
    }
    
    return %bams_groups;
}

=head2 create_release_files_requires

 Title   : create_release_files_requires
 Usage   : my $required_files = $obj->create_release_files_requires('/path/to/lane');
 Function: Find out what files the create_release_files action needs before
           it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub create_release_files_requires {
    my $self = shift;
    return ['.platform_merge_done'];
}

=head2 create_release_files_provides

 Title   : create_release_files_provides
 Usage   : my $provided_files = $obj->create_release_files_provides('/path/to/lane');
 Function: Find out what files the create_release_files action generates on
           success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub create_release_files_provides {
    my ($self, $lane_path) = @_;
    my @files = ('.create_release_files_done');
    if ($self->{dcc_mode}) {
        push(@files, 'release_files.fofn', 'release_files.md5');
    }
    return \@files;
}

=head2 create_release_files

 Title   : create_release_files
 Usage   : $obj->create_release_files('/path/to/lane', 'lock_filename');
 Function: Creates the set of platform-level files that form the release. That
           is, .bai and .bas files are made for each platform-level bam. .md5
           files are made for all 3. Also makes softlinks to the release files
           prefixed with 'release' (release.bam, release.bam.bai,
           release.bam.bas), so that it is consistent to find them.

           If the dcc_mode option is supplied, will also do per-chr splits and
           rename these to DCC-style filenames. Otherwise, if do_chr_splits was
           supplied, will just do the splits and leave them with Sanger-style
           filenames. In either case, .bai and .bas files are made for each
           per-chr bam.

           At the end, if dcc_mode option was supplied, will create two
           files in the root of the release directory: release_files.fofn
           (listing all the DCC-filenamed release files minus the md5s) and
           release_files.md5 (containing all the md5 checksums for the release
           files).
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub create_release_files {
    my ($self, $lane_path, $action_lock) = @_;
    
    my $fofn = $self->{fsu}->catfile($lane_path, '.platform_merge_done');
    my @in_bams = $self->{io}->parse_fofn($fofn, $lane_path);
    
    my $out_fofn = $self->{fsu}->catfile($lane_path, '.create_release_files_expected');
    unlink($out_fofn);
    if ($self->{dcc_mode}) {
        unlink($self->{fsu}->catfile($lane_path, 'release_files.fofn'));
        unlink($self->{fsu}->catfile($lane_path, 'release_files.md5'));
    }
    my $vuh = VertRes::Utils::Hierarchy->new();
    
    my $verbose = $self->verbose;
    
    my $orig_bsub_opts = $self->{bsub_opts};
    $self->{bsub_opts} = '-q long';
    
    my @release_files;
    foreach my $bam (@in_bams) {
        $bam = $self->{fsu}->catfile($lane_path, $bam);
        my ($basename, $path) = fileparse($bam);
        my $release_name = $self->{fsu}->catfile($path, 'release.bam');
        
        # md5 of bam & links
        my $bam_md5 = $bam.'.md5';
        unless (-s $bam_md5) {
            LSF::run($action_lock, $lane_path, $self->{prefix}.'create_release_files', $self,
                     qq{md5sum $bam > $bam_md5; ln -s $basename $release_name});
        }
        elsif (! -e $release_name) {
            symlink($basename, $release_name);
        }
        push(@release_files, $bam, $bam_md5);
        
        # bai & its md5 & links
        my $bai = $bam.'.bai';
        my $bai_md5 = $bai.'.md5';
        unless (-s $bai && -s $bai_md5) {
            LSF::run($action_lock, $lane_path, $self->{prefix}.'create_release_files', $self,
                     qq{perl -MVertRes::Wrapper::samtools -Mstrict -e "VertRes::Wrapper::samtools->new(verbose => $verbose)->index(qq[$bam], qq[$bai.tmp]); die qq[index failed for $bam\n] unless -s qq[$bai.tmp]; system(qq[mv $bai.tmp $bai; md5sum $bai > $bai_md5; ln -s $basename.bai $release_name.bai]);"});
        }
        elsif (! -e "$release_name.bai") {
            symlink("$basename.bai", "$release_name.bai");
        }
        push(@release_files, $bai, $bai_md5);
        
        # bas & its md5 & links
        my $bas = $bam.'.bas';
        my $bas_md5 = $bas.'.md5';
        unless (-s $bas && -s $bas_md5) {
            LSF::run($action_lock, $lane_path, $self->{prefix}.'create_release_files', $self,
                     qq{perl -MVertRes::Utils::Sam -Mstrict -e "VertRes::Utils::Sam->new(verbose => $verbose)->bas(qq[$bam], qq[$self->{release_date}], qq[$bas]); die qq[bas failed for $bam\n] unless -s qq[$bas]; system(qq[md5sum $bas > $bas_md5; ln -s $basename.bas $release_name.bas]);"}); # , qq[$self->{sequence_index}] bas() needs a database option?
        }
        elsif (! -e "$release_name.bas") {
            symlink("$basename.bas", "$release_name.bas");
        }
        push(@release_files, $bas, $bas_md5);
        
        # chr splits
        if ($self->{do_chr_splits}) {
            my $su = VertRes::Utils::Sam->new();
            my @expected_split_bams = $su->split_bam_by_sequence($bam, pretend => 1);
            my @orig_split_bams = @expected_split_bams;
            
            my $bams = 0;
            foreach my $ebam (@expected_split_bams) {
                push(@release_files, $ebam) unless $self->{dcc_mode};
                $bams += -s $ebam ? 1 : 0;
                
                unless ($self->{dcc_mode}) {
                    foreach my $suffix ('.md5', '.bai', '.bai.md5', '.bas', '.bas.md5') {
                        push(@release_files, $ebam.$suffix);
                    }
                }
            }
            
            my $dcc_bams = 0;
            if ($self->{dcc_mode}) {
                
                my @expected_dcc_bams;
                foreach my $ebam (@expected_split_bams) {
                    my $basename = basename($ebam);
                    my ($chr) = $basename =~ /^([^\.]+)/;
                    my $dcc_filename = $vuh->dcc_filename($bam, $self->{release_date}, $chr).'.bam';
                    
                    my $dccbam = $ebam;
                    $dccbam =~ s/$basename$/$dcc_filename/;
                    push(@release_files, $dccbam);
                    push(@expected_dcc_bams, $dccbam);
                    $dcc_bams += -s $dccbam ? 1 : 0;
                    
                    foreach my $suffix ('.md5', '.bai', '.bai.md5', '.bas', '.bas.md5') {
                        push(@release_files, $dccbam.$suffix);
                    }
                }
                @expected_split_bams = @expected_dcc_bams;
            }
            else {
                $dcc_bams = $bams;
            }
            
            unless ($dcc_bams == @expected_split_bams) {
                if ($bams == @expected_split_bams && $self->{dcc_mode}) {
                    # rename to dcc
                    foreach my $i (0..$#orig_split_bams) {
                        move($orig_split_bams[$i], $expected_split_bams[$i]);
                    }
                }
                elsif (-s $bai) {
                    LSF::run($action_lock, $lane_path, $self->{prefix}.'create_release_files', $self,
                             qq{perl -MVertRes::Utils::Sam -Mstrict -e "VertRes::Utils::Sam->new(verbose => $verbose)->split_bam_by_sequence(qq[$bam]);"});
                }
            }
            else {
                $self->{bsub_opts} = '-q normal';
                foreach my $ebam (@expected_split_bams) {
                    # md5
                    my $emd5 = $ebam.'.md5';
                    unless (-s $emd5) {
                        $self->{bsub_opts} = '-q small';
                        LSF::run($action_lock, $lane_path, $self->{prefix}.'create_release_files', $self, qq{md5sum $ebam > $emd5});
                        $self->{bsub_opts} = '-q normal';
                    }
                    
                    # bai & its md5
                    my $ebai = $ebam.'.bai';
                    my $ebai_md5 = $ebai.'.md5';
                    unless (-s $ebai && -s $ebai_md5) {
                        LSF::run($action_lock, $lane_path, $self->{prefix}.'create_release_files', $self,
                                 qq{perl -MVertRes::Wrapper::samtools -Mstrict -e "VertRes::Wrapper::samtools->new(verbose => $verbose)->index(qq[$ebam], qq[$ebai.tmp]); die qq[index failed for $ebam\n] unless -s qq[$ebai.tmp]; system(qq[mv $ebai.tmp $ebai; md5sum $ebai > $ebai_md5]);"});
                    }
                    
                    # bas & its md5
                    my $ebas = $ebam.'.bas';
                    my $ebas_md5 = $ebas.'.md5';
                    unless (-s $ebas && -s $ebas_md5) {
                        LSF::run($action_lock, $lane_path, $self->{prefix}.'create_release_files', $self,
                                 qq{perl -MVertRes::Utils::Sam -Mstrict -e "VertRes::Utils::Sam->new(verbose => $verbose)->bas(qq[$ebam], qq[$self->{release_date}], qq[$ebas]); die qq[bas failed for $ebam\n] unless -s qq[$ebas]; system(qq[md5sum $ebas > $ebas_md5]);"});
                    }
                }
            }
        }
    }
    
    open(my $ofh, '>', $out_fofn) || $self->throw("Couldn't write to $out_fofn");
    foreach my $file (@release_files) {
        print $ofh $file, "\n";
    }
    my $expected = @release_files;
    print $ofh "# expecting $expected\n";
    close($ofh);
    
    $self->{bsub_opts} = $orig_bsub_opts;
    
    return $self->{NO};
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
    
    foreach my $file (qw(log job_status)) {
        unlink($self->{fsu}->catfile($lane_path, $prefix.$file));
    }
    
    my $file_base = $self->{fsu}->catfile($lane_path, $prefix);
    foreach my $job_base (qw(library_merge lib_markdup platform_merge create_release_files sample_merge)) {
        foreach my $suffix ('o', 'e') {
            unlink("$file_base$job_base.$suffix");
        }
    }
    
    return $self->{Yes};
}

sub is_finished {
    my ($self, $lane_path, $action) = @_;
    
    my $action_name = $action->{name};
    
    if ($action_name eq 'cleanup') {
        return $self->{No};
    }
    elsif ($action_name eq 'lib_markdup' ||
            $action_name eq 'library_merge' ||
            $action_name eq 'platform_merge' ||
            $action_name eq 'sample_merge' ||
            $action_name eq 'create_release_files') {
        my $expected_file = $self->{fsu}->catfile($lane_path, ".${action_name}_expected");
        my $not_yet_done = -e $expected_file;
        my $done_file = $self->{fsu}->catfile($lane_path, ".${action_name}_done");
        $self->_merge_check($expected_file, $done_file, $lane_path);
        
        if ($action_name eq 'lib_markdup' && $not_yet_done && -s $done_file) {
            # we've just completed marking dups; delete the original unmarked
            # bams to save space
            my $fofn = $self->{fsu}->catfile($lane_path, '.library_merge_done');
            my @files = $self->{io}->parse_fofn($fofn, $lane_path);
            foreach my $file (@files) {
                my $marked = $file;
                $marked =~ s/\.bam$/.markdup.bam/;
                if (-s $file && -s $marked) {
                    unlink($file);
                }
            }
        }
    }
    
    if ($action_name eq 'create_release_files' &&
            $self->{dcc_mode} &&
            -s $self->{fsu}->catfile($lane_path, ".create_release_files_done") &&
            ! -s $self->{fsu}->catfile($lane_path, "release_files.fofn")) {
        my $vuh = VertRes::Utils::Hierarchy->new();
        
        # make a release_files.fofn with all the bams and other files sans the
        # md5s, and cat all the md5s into 1 file
        my $rfofn = $self->{fsu}->catfile($lane_path, "release_files.fofn.tmp");
        open(my $rfofnfh, '>', $rfofn) || $self->throw("Couldn't write to $rfofn");
        my $rmd5 = $self->{fsu}->catfile($lane_path, "release_files.md5.tmp");
        open(my $rmd5fh, '>', $rmd5) || $self->throw("Couldn't write to $rmd5");
        
        my $rdone = $self->{fsu}->catfile($lane_path, ".create_release_files_done");
        open(my $rdfh, $rdone) || $self->throw("Couldn't open $rdone");
        local $| = 1;
        my $spin = 0;
        while (<$rdfh>) {
            chomp;
            /\S/ || next;
            /^#/ && next;
            my $file = $_;
            my ($base, $path) = fileparse($file);
            
            # only interested in splits
            $base =~ /\.(?:chrom(?:\d+|[XY]|MT)|nonchrom|unmapped)\./ || next;
            
            $spin++;
            if ($spin == 1) {
                print "\b-";
            }
            elsif ($spin == 2) {
                print "\b\\";
            }
            elsif ($spin == 3) {
                print "\b|";
            }
            elsif ($spin == 4) {
                print "\b/";
            }
            elsif ($spin == 5) {
                print "\b-";
                $spin = 0;
            }
            
            if ($file =~ /md5$/) {
                open(my $md5fh, $file) || $self->throw("Couldn't open $_");
                my $line = <$md5fh>;
                my ($md5, $path) = split(' ', $line);
                $base =~ s/\.md5$//;
                print $rmd5fh $md5, "\t", $base, "\n";
                close($md5fh);
            }
            else {
                print $rfofnfh $file, "\n";
            }
        }
        close($rdfh);
        print $rfofnfh $self->{fsu}->catfile($lane_path, "release_files.md5"), "\n";
        close($rmd5fh);
        close($rfofnfh);
        move($rfofn, $self->{fsu}->catfile($lane_path, "release_files.fofn"));
        move($rmd5, $self->{fsu}->catfile($lane_path, "release_files.md5"));
    }
    
    return $self->SUPER::is_finished($lane_path, $action);
}

sub _merge_check {
    my ($self, $expected_file, $done_file, $lane_path) = @_;
    
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
            if (-s $_) {
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
        elsif ($written_expected == $expected_bams && $self->{skip_fails}) {
            my $new_expected = @bams;
            my $diff = $expected_bams - $done_bams;
            $self->warn("$diff files are being skipped! See skipped_files.fofn in your release directory.\n");
            
            open(my $dfh, '>', $done_file) || $self->throw("Couldn't write to $done_file");
            foreach my $bam (@bams) {
                print $dfh $bam, "\n";
            }
            print $dfh "# expecting $new_expected\n";
            close($dfh);
            
            my $skip_file = $self->{fsu}->catfile($lane_path, 'skipped_files.fofn');
            open(my $sfh, '>>', $skip_file) || $self->throw("Couldn't write to $skip_file");
            foreach my $bam (@skipped) {
                print $sfh $bam, "\n";
            }
            close($sfh);
        }
    }
}

# we override running_status to allow *_merge() to be called even while some jobs
# are still running, because it is doing 200 jobs at a time and we don't want
# to wait for all 200 to finish before the next 200 are scheduled.
sub running_status {
    my ($self, $jids_file) = @_;
    if ($jids_file =~ /_merge/) {
        return $LSF::No;
    }
    
    return $self->SUPER::running_status($jids_file);
}

1;

