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
    
    do_recalibration => 1,
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

# make another file that simply also contains the release directory:
echo "/path/to/rel_dir" > rel.fod

# make another config file that mentions the previous two files:
echo "<rel.fod release.conf" > release.pipeline

# make your release directory if it doesn't already exist:
mkdir /path/to/rel_dir

# run the pipeline:
run-pipeline -c release.pipeline -v

# (and make sure it keeps running by adding that last to a regular cron job)

=head1 DESCRIPTION

A module for creating a release (merging, splitting, recalibration etc.) from
lanes mapped by VertRes::Pipelines::Mapping.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Pipelines::Release;

use strict;
use warnings;
use VertRes::Utils::Hierarchy;
use VertRes::IO;
use VRTrack::VRTrack;
use VRTrack::Lane;
use File::Basename;
use File::Spec;
use File::Copy;
use Cwd 'abs_path';
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
                { name     => 'lib_rmdup',
                  action   => \&lib_rmdup,
                  requires => \&lib_rmdup_requires, 
                  provides => \&lib_rmdup_provides },
                { name     => 'platform_merge',
                  action   => \&platform_merge,
                  requires => \&platform_merge_requires, 
                  provides => \&platform_merge_provides },
                { name     => 'recalibrate',
                  action   => \&recalibrate,
                  requires => \&recalibrate_requires, 
                  provides => \&recalibrate_provides },
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
                do_recalibration => 0,
                do_sample_merge => 0,
                bsub_opts => '',
                dont_wait => 1);

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
                                     or individual-level bams by chr)
           do_recalibration => boolean (default false: don't recalibrate any
                                        bams; when true, recalibrates the
                                        platform-level bams, which will be used
                                        to make the individual-level bams)
           dcc_hardlinks => boolean (default false: don't create hardlinks to
                                     release files with DCC-style filenames)
           do_sample_merge => boolean (default false: don't create sample-level
                                       bams)
           skip_fails => boolean (default false; when true, if some jobs fail,
                                  continue to the next action anyway)
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
    
    $self->{io} = VertRes::IO->new;
    
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

=head2 platform_merge_requires

 Title   : platform_merge_requires
 Usage   : my $required_files = $obj->platform_merge_requires('/path/to/lane');
 Function: Find out what files the platform_merge action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub platform_merge_requires {
    my $self = shift;
    return ['.lib_rmdup_done'];
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
                              '.lib_rmdup_done',
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
    return ['.release_files_done'];
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
    
    #*** probably going to have an issue with relative paths in .recalibrate_done ?
    
    $self->merge_up_one_level($lane_path,
                              $action_lock,
                              $self->{do_recalibration} ? '.recalibrate_done' : '.platform_merge_done',
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
    $fofn = $self->{io}->catfile($lane_path, $fofn);
    $out_fofn = $self->{io}->catfile($lane_path, $out_fofn);
    
    unlink($out_fofn);
    
    my $orig_bsub_opts = $self->{bsub_opts};
    $queue ||= 'normal';
    $self->{bsub_opts} = '-q '.$queue.' -M5000000 -R \'select[mem>5000] rusage[mem=5000]\'';
    
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

=head2 recalibrate_requires

 Title   : recalibrate_requires
 Usage   : my $required_files = $obj->recalibrate_requires('/path/to/lane');
 Function: Find out what files the recalibrate action needs before
           it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub recalibrate_requires {
    my $self = shift;
    return ['.platform_merge_done'];
}

=head2 recalibrate_provides

 Title   : recalibrate_provides
 Usage   : my $provided_files = $obj->recalibrate_provides('/path/to/lane');
 Function: Find out what files the recalibrate action generates on
           success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub recalibrate_provides {
    my ($self, $lane_path) = @_;
    return $self->{do_recalibration} ? ['.recalibrate_done'] : ['.platform_merge_done'];
}

=head2 recalibrate

 Title   : recalibrate
 Usage   : $obj->recalibrate('/path/to/lane', 'lock_filename');
 Function: Recalibrate the quality values in platform-level bams.
           Doesn't run unless do_recalibration has been supplied as a config
           option.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub recalibrate {
    my ($self, $lane_path, $action_lock) = @_;
    return $self->{Yes} unless $self->{do_recalibration};
    
    my $fofn = $self->{io}->catfile($lane_path, '.platform_merge_done');
    my @in_bams = $self->{io}->parse_fofn($fofn, $lane_path);
    
    my $out_fofn = $self->{io}->catfile($lane_path, '.recalibrate_expected');
    unlink($out_fofn);
    
    my $verbose = $self->verbose;
    
    my $orig_bsub_opts = $self->{bsub_opts};
    $self->{bsub_opts} = '-q normal -M6000000 -R \'select[mem>6000] rusage[mem=6000]\'';
    
    my @out_bams;
    foreach my $bam (@in_bams) {
        $bam = $self->{io}->catfile($lane_path, $bam);
        my $out_bam = $bam;
        $out_bam =~ s/\.bam$/.recal.bam/;
        push(@out_bams, $out_bam);
        
        next if -s $out_bam;
        
        LSF::run($action_lock, $lane_path, $self->{prefix}.'platform_recalibration', $self,
                 qq{perl -MVertRes::Wrapper::GATK -Mstrict -e "VertRes::Wrapper::GATK->new(verbose => $verbose)->recalibrate(qq[$bam], qq[$out_bam]); die qq[recalibration failed for $bam\n] unless -s qq[$out_bam];"});
    }
    
    open(my $ofh, '>', $out_fofn) || $self->throw("Couldn't write to $out_fofn");
    foreach my $out_bam (@out_bams) {
        print $ofh $out_bam, "\n";
    }
    my $expected = @out_bams;
    print $ofh "# expecting $expected\n";
    close($ofh);
    
    $self->{bsub_opts} = $orig_bsub_opts;
    return $self->{NO};
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
    return $self->{do_recalibration} ? ['.recalibrate_done'] : ['.platform_merge_done'];
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
    if ($self->{dcc_hardlinks}) {
        push(@files, 'release_files.fofn', 'release_files.md5');
    }
    return \@files;
}

=head2 create_release_files

 Title   : create_release_files
 Usage   : $obj->create_release_files('/path/to/lane', 'lock_filename');
 Function: Creates the set of platform-level files that form the release and
           would be uploaded to the DCC. That is, .bai and .bas files are made
           for each platform-level bam. .md5 files are made for all 3. If the
           dcc_hardlinks option is supplied, will also make hardlinks to these
           release files to give them DCC-style filenames. Always also makes
           softlinks to the release files prefixed with 'release' (release.bam,
           release.bam.bai, release.bam.bas), so that it is consistent to find
           them, regardless of if recalibration was done or not.

           At the end, if dcc_hardlinks option was supplied, will create two
           files in the root of the release directory: release_files.fofn
           (listing all the DCC-filenamed release files minus the md5s) and
           release_files.md5 (containing all the md5 checksums for the release
           files).
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub create_release_files {
    my ($self, $lane_path, $action_lock) = @_;
    
    my $fofn = $self->{io}->catfile($lane_path, $self->{do_recalibration} ? '.recalibrate_done' : '.platform_merge_done');
    my @in_bams = $self->{io}->parse_fofn($fofn, $lane_path);
    
    my $out_fofn = $self->{io}->catfile($lane_path, '.create_release_files_expected');
    unlink($out_fofn);
    if ($self->{dcc_hardlinks}) {
        unlink($self->{io}->catfile($lane_path, 'release_files.fofn'));
        unlink($self->{io}->catfile($lane_path, 'release_files.md5'));
    }
    
    my $verbose = $self->verbose;
    
    my @release_files;
    foreach my $bam (@in_bams) {
        $bam = $self->{io}->catfile($lane_path, $bam);
        my ($basename, $path) = fileparse($bam);
        my $release_name = File::Spec->catfile($path, 'release.bam');
        
        # md5 of bam & links
        my $bam_md5 = $bam.'.md5';
        unless (-s $bam_md5) {
            LSF::run($action_lock, $lane_path, $self->{prefix}.'create_release_files', $self,
                     qq{md5sum $bam > $bam_md5; ln -s $basename $release_name});
        }
        push(@release_files, $bam, $bam_md5);
        
        # bai & its md5 & links
        my $bai = $bam.'.bai';
        my $bai_md5 = $bai.'.md5';
        unless (-s $bai && -s $bai_md5) {
            LSF::run($action_lock, $lane_path, $self->{prefix}.'create_release_files', $self,
                     qq{perl -MVertRes::Wrapper::samtools -Mstrict -e "VertRes::Wrapper::samtools->new(verbose => $verbose)->index(qq[$bam], qq[$bai]); die qq[index failed for $bam\n] unless -s qq[$bai]; system(qq[md5sum $bai > $bai_md5; ln -s $basename.bai $release_name.bai]);"});
        }
        push(@release_files, $bai, $bai_md5);
        
        # bas & its md5 & links
        my $bas = $bam.'.bas';
        my $bas_md5 = $bas.'.md5';
        unless (-s $bas && -s $bas_md5) {
            LSF::run($action_lock, $lane_path, $self->{prefix}.'create_release_files', $self,
                     qq{perl -MVertRes::Utils::Sam -Mstrict -e "VertRes::Utils::Sam->new(verbose => $verbose)->bas(qq[$bam], qq[$bas]); die qq[bas failed for $bam\n] unless -s qq[$bas]; system(qq[md5sum $bas > $bas_md5; ln -s $basename.bas $release_name.bas]);"}); # , qq[$self->{sequence_index}] bas() needs a database option?
        }
        push(@release_files, $bas, $bas_md5);
    }
    
    open(my $ofh, '>', $out_fofn) || $self->throw("Couldn't write to $out_fofn");
    foreach my $file (@release_files) {
        print $ofh $file, "\n";
    }
    my $expected = @release_files;
    print $ofh "# expecting $expected\n";
    close($ofh);
    
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
        unlink($self->{io}->catfile($lane_path, $prefix.$file));
    }
    
    my $file_base = $self->{io}->catfile($lane_path, $prefix);
    foreach my $job_base (qw(library_merge lib_rmdup platform_merge platform_recalibration create_release_files sample_merge)) {
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
    elsif ($action_name eq 'lib_rmdup' ||
            $action_name eq 'library_merge' ||
            $action_name eq 'platform_merge' ||
            $action_name eq 'sample_merge' ||
            $action_name eq 'recalibrate' ||
            $action_name eq 'create_release_files') {
        my $expected_file = $self->{io}->catfile($lane_path, ".${action_name}_expected");
        my $done_file = $self->{io}->catfile($lane_path, ".${action_name}_done");
        $self->_merge_check($expected_file, $done_file, $lane_path);
    }
    
    if ($action_name eq 'create_release_files' &&
            $self->{dcc_hardlinks} &&
            -s $self->{io}->catfile($lane_path, ".create_release_files_done") &&
            ! -s $self->{io}->catfile($lane_path, "release_files.fofn")) {
        my $vuh = VertRes::Utils::Hierarchy->new();
        
        # make a release_files.fofn with all the bams and other files sans the
        # md5s, and cat all the md5s into 1 file
        my $rfofn = $self->{io}->catfile($lane_path, "release_files.fofn.tmp");
        open(my $rfofnfh, '>', $rfofn) || $self->throw("Couldn't write to $rfofn");
        my $rmd5 = $self->{io}->catfile($lane_path, "release_files.md5.tmp");
        open(my $rmd5fh, '>', $rmd5) || $self->throw("Couldn't write to $rmd5");
        
        my $rdone = $self->{io}->catfile($lane_path, ".create_release_files_done");
        open(my $rdfh, $rdone) || $self->throw("Couldn't open $rdone");
        my $dcc_filename;
        my $orig_name;
        local $| = 1;
        print "making links ";
        my $doti = 0;
        while (<$rdfh>) {
            chomp;
            /\S/ || next;
            /^#/ && next;
            my $file = $_;
            my ($base, $path) = fileparse($file);
            
            $doti++;
            if ($doti == 1) {
                print "\b\b\b.  ";
            }
            elsif ($doti == 2) {
                print "\b\b.. ";
            }
            if ($doti == 3) {
                print "\b...";
                $doti = 0;
            }
            
            if (/bam$/) {
                $dcc_filename = $vuh->dcc_filename($file).'.bam';
                $orig_name = $base;
            }
            
            $base =~ s/$orig_name/$dcc_filename/;
            
            if ($file =~ /md5$/) {
                open(my $md5fh, $file) || $self->throw("Couldn't open $_");
                my $line = <$md5fh>;
                my ($md5, $path) = split(' ', $line);
                $base =~ s/\.md5$//;
                print $rmd5fh $md5, "\t", $base, "\n";
                close($md5fh);
            }
            else {
                # make hardlink to DCC-style filename
                my $dcc_path = File::Spec->catfile($path, $base);
                unlink($dcc_path);
                system("ln $file $dcc_path");
                print $rfofnfh $dcc_path, "\n";
            }
        }
        close($rdfh);
        print $rfofnfh $self->{io}->catfile($lane_path, "release_files.md5"), "\n";
        close($rmd5fh);
        close($rfofnfh);
        move($rfofn, $self->{io}->catfile($lane_path, "release_files.fofn"));
        move($rmd5, $self->{io}->catfile($lane_path, "release_files.md5"));
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
            
            my $skip_file = $self->{io}->catfile($lane_path, 'skipped_files.fofn');
            open(my $sfh, '>>', $skip_file) || $self->throw("Couldn't write to $skip_file");
            foreach my $bam (@skipped) {
                print $sfh $bam, "\n";
            }
            close($sfh);
        }
    }
}

1;

