=head1 NAME

VertRes::Pipelines::MergeUp - pipeline for creating sample bams from lane bams

=head1 SYNOPSIS

# make a file of absolute paths to the lane bams, eg:
lane_bams.pl --root /abs/META --db g1k_meta --all_samples
             --slx_mapper bwa --454_mapper ssaha2 --assembly_name NCBI37
             --ignore_platform SOLID --improved > lane_bams.fofn

# make a conf file with root pointing to where you want the merged bams built,
# and setting the lanes, mapper and assembly_name options. Optional settings
# also go here. Your mergeup.conf file may look like:
root    => '/path/to/merge_dir',
module  => 'VertRes::Pipelines::MergeUp',
prefix  => '_',
lane_bams => 'lane_bams.fofn',
previous_merge_root => '/path/to/previous/merge_root',
data => {
    tag_strip => { strip => [qw(OQ MD XM XG XO)] },
    simultaneous_merges => 200,
    do_cleanup => 1
}

# other options you could add to the mergeup.conf data {} section include:
# do_sample_merge => 1 (default is to platform level only)
# extract_intervals => {}

# outside the data {} section you can also specify:
# previous_merge_root => '/path/to/old/merge_root'
# the previous_merge_root will result in the pipeline checking if a particular
# bam needs to be (re)made depending on if lanes have been removed/added/altered
# since that previous merge. If a bam doesn't need to be made, the new
# merge directory will be a symlink to the old merge directory.

# make another config file that simply mentions the merge directory and the
# previous config file:
echo "__MERGEUP__ mergeup.conf" > mergeup.pipeline

# make your mergeup directory if it doesn't already exist:
mkdir /path/to/merge_dir

# run the pipeline:
run-pipeline -c mergeup.pipeline -v

# (and make sure it keeps running by adding that last to a regular cron job)

=head1 DESCRIPTION

A pipeline module for merging lane bams to the platform (optionaly, sample)
level, including marking of duplicates at the library level and optional
interval extraction for exome work. It also strips most tags from each record
to minimise final bam size.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Pipelines::MergeUp;

use strict;
use warnings;
use VertRes::Utils::Hierarchy;
use VertRes::IO;
use VertRes::Utils::Sam;
use VertRes::Utils::FileSystem;
use File::Basename;
use File::Spec;
use File::Copy;
use Cwd 'abs_path';
use LSF;

use base qw(VertRes::Pipeline);

our $actions = [{ name     => 'create_hierarchy',
                  action   => \&create_hierarchy,
                  requires => \&create_hierarchy_requires, 
                  provides => \&create_hierarchy_provides },
                { name     => 'library_merge',
                  action   => \&library_merge,
                  requires => \&library_merge_requires, 
                  provides => \&library_merge_provides },
                { name     => 'tag_strip',
                  action   => \&tag_strip,
                  requires => \&tag_strip_requires, 
                  provides => \&tag_strip_provides },
                { name     => 'lib_markdup',
                  action   => \&lib_markdup,
                  requires => \&lib_markdup_requires, 
                  provides => \&lib_markdup_provides },
                { name     => 'extract_intervals',
                  action   => \&extract_intervals,
                  requires => \&extract_intervals_requires, 
                  provides => \&extract_intervals_provides},
                { name     => 'platform_merge',
                  action   => \&platform_merge,
                  requires => \&platform_merge_requires, 
                  provides => \&platform_merge_provides },
                { name     => 'sample_merge',
                  action   => \&sample_merge,
                  requires => \&sample_merge_requires, 
                  provides => \&sample_merge_provides },
                { name     => 'cleanup',
                  action   => \&cleanup,
                  requires => \&cleanup_requires, 
                  provides => \&cleanup_provides }];

our %options = (do_cleanup => 0,
                do_sample_merge => 0,
                tag_strip => { strip => [qw(OQ MD XM XG XO)] },
                simultaneous_merges => 200,
                bsub_opts => '',
                dont_wait => 1);

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Pipelines::MergeUp->new(lane => '/path/to/lane');
 Function: Create a new VertRes::Pipelines::MergeUp object.
 Returns : VertRes::Pipelines::MergeUp object
 Args    : hierarchy_root => '/abs/path' (REQUIRED; the root of the hierarchy)
           lane_bams => fofn (REQUIRED; a file listing all the bams you want
                              to merge)
           
           tag_strip => { strip => [qw(OQ MD XM XG XO)] } (default as shown;
                          args as per VertRes::Utils::Sam::tag_strip)
           do_cleanup => boolean (default false: don't do the cleanup action)
           simultaneous_merges => int (default 200; the number of merge jobs to
                                       do at once - limited to avoid IO
                                       problems)
           do_sample_merge => boolean (default false: don't create sample-level
                                       bams)
           extract_intervals => {} (after marking duplicates, extract reads for
                                    only these intervals; default keep all
                                    reads)
           other optional args as per VertRes::Pipeline

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(%options, actions => $actions, @args);
    
    $self->{hierarchy_root} || $self->throw("hierarchy_root path not supplied, can't continue");
    $self->{lane_bams} || $self->throw("lane_bams fofn not supplied, can't continue");
    
    $self->{io} = VertRes::IO->new;
    $self->{fsu} = VertRes::Utils::FileSystem->new;
    
    return $self;
}

=head2 create_hierarchy_requires

 Title   : create_hierarchy_requires
 Usage   : my $required_files = $obj->create_hierarchy_requires('/path/to/lane');
 Function: Find out what files the create_hierarchy action needs before
           it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub create_hierarchy_requires {
    my $self = shift;
    return [];
}

=head2 create_hierarchy_provides

 Title   : create_hierarchy_provides
 Usage   : my $provided_files = $obj->create_hierarchy_provides('/path/to/lane');
 Function: Find out what files the create_hierarchy action generates on
           success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub create_hierarchy_provides {
    my ($self, $lane_path) = @_;
    return ['.hierarchy_made'];
}

=head2 create_hierarchy

 Title   : create_hierarchy
 Usage   : $obj->create_hierarchy('/path/to/lane', 'lock_filename');
 Function: Create a hierarchy with the desired lane-level bams symlinked in.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub create_hierarchy {
    my ($self, $sample_root, $action_lock) = @_;
    
    my @lane_bams = $self->{io}->parse_fofn($self->{lane_bams}, '/');
    
    my $hu = VertRes::Utils::Hierarchy->new;
    
    my @linked_bams;
    foreach my $bam (@lane_bams) {
        my $new_bam = $hu->lane_bam_link($bam, $self->{hierarchy_root});
        push(@linked_bams, $new_bam);
    }
    
    my $done_file = $self->{fsu}->catfile($sample_root, '.hierarchy_made');
    open(my $dfh, '>', $done_file) || $self->throw("Could not write to $done_file");
    foreach my $bam (@linked_bams) {
        print $dfh $bam, "\n";
    }
    my $expected = @linked_bams;
    print $dfh "# expecting $expected\n";
    close($dfh);
    
    # just check that the write was successful
    sleep(5);
    my @got = $self->{io}->parse_fofn($done_file);
    unless (@got == @linked_bams) {
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
    return ['.hierarchy_made'];
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
                              '.hierarchy_made',
                              undef,
                              'library_merge',
                              '.library_merge_expected',
                              'long');
    
    return $self->{No};
}


=head2 tag_strip_requires

 Title   : tag_strip_requires
 Usage   : my $required_files = $obj->tag_strip_requires('/path/to/lane');
 Function: Find out what files the tag_strip action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub tag_strip_requires {
    my $self = shift;
    return ['.library_merge_done'];
}

=head2 tag_strip_provides

 Title   : tag_strip_provides
 Usage   : my $provided_files = $obj->tag_strip_provides('/path/to/lane');
 Function: Find out what files the tag_strip action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub tag_strip_provides {
    my $self = shift;
    return ['.tag_strip_done'];
}

=head2 tag_strip

 Title   : tag_strip
 Usage   : $obj->tag_strip('/path/to/lane', 'lock_filename');
 Function: Strips unecessary tags from each record of the library level bam
           file.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub tag_strip {
    my ($self, $lane_path, $action_lock) = @_;
    
    my $fofn = $self->{fsu}->catfile($lane_path, '.library_merge_done');
    my @files = $self->{io}->parse_fofn($fofn, $lane_path);
    my $verbose = $self->verbose();
    
    my $tag_strip_args = '';
    foreach my $type ('strip', 'keep') {
        my $vals = $self->{tag_strip}->{$type} || next;
        $tag_strip_args .= " $type => [qw(@{$vals})]";
    }
    
    my @strip_bams;
    foreach my $merge_bam (@files) {
        $merge_bam = $self->{fsu}->catfile($lane_path, $merge_bam);
        my($basename, $path) = fileparse($merge_bam);
        
        my $strip_bam = $merge_bam;
        $strip_bam =~ s/\.bam$/.strip.bam/;
        push(@strip_bams, $strip_bam);
        next if -e $strip_bam;
        
        # if a higher-level bam already exists, don't repeat making this level
        # bam if we deleted it
        next if $self->_skip_if_higher_level_bam_present($path);
        next unless -s $merge_bam;
        
        my $job_name = $self->{prefix}.'tag_strip';
        $self->archive_bsub_files($path, $job_name);
        $job_name = $self->{fsu}->catfile($path, $job_name);
        
        LSF::run($action_lock, $lane_path, $job_name, $self,
                 qq{perl -MVertRes::Utils::Sam -Mstrict -e "VertRes::Utils::Sam->new(verbose => $verbose)->tag_strip(qq[$merge_bam], qq[$strip_bam], $tag_strip_args) || die qq[tag_strip failed for $merge_bam\n];"});
    }
    
    my $out_fofn = $self->{fsu}->catfile($lane_path, '.tag_strip_expected');
    open(my $ofh, '>', $out_fofn) || $self->throw("Couldn't write to $out_fofn");
    foreach my $out_bam (@strip_bams) {
        print $ofh $out_bam, "\n";
    }
    my $expected = @strip_bams;
    print $ofh "# expecting $expected\n";
    close($ofh);
    
    return $self->{No};
}

sub _skip_if_higher_level_bam_present {
    my ($self, $path) = @_;
    
    my @dirs = File::Spec->splitdir($path);
    pop(@dirs);
    pop(@dirs);
    my $higher_bam = File::Spec->catfile(@dirs, 'raw.bam');
    if ($self->{fsu}->file_exists($higher_bam) && ! -l $higher_bam) {
        return 1;
    }
    
    return 0;
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
    return ['.tag_strip_done'];
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
    
    my $fofn = $self->{fsu}->catfile($lane_path, '.tag_strip_done');
    my @files = $self->{io}->parse_fofn($fofn, $lane_path);
    my $verbose = $self->verbose();
    
    my $orig_bsub_opts = $self->{bsub_opts};
    $self->{bsub_opts} = '-q long -M6100000 -R \'select[mem>6100] rusage[mem=6100]\'';
    
    my @markdup_bams;
    foreach my $merge_bam (@files) {
        $merge_bam = $self->{fsu}->catfile($lane_path, $merge_bam);
        my($basename, $path) = fileparse($merge_bam);
        
        my $markdup_bam = $merge_bam;
        $markdup_bam =~ s/\.bam$/.markdup.bam/;
        push(@markdup_bams, $markdup_bam);
        next if -e $markdup_bam;
        
        # if a higher-level bam already exists, don't repeat making this level
        # bam if we deleted it
        next if $self->_skip_if_higher_level_bam_present($path);
        next unless -s $merge_bam;
        
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

=head2 extract_intervals_requires

 Title   : extract_intervals_requires
 Usage   : my $required_files = $obj->extract_intervals_requires('/path/to/lane');
 Function: Find out what files the extract_intervals action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub extract_intervals_requires {
    my $self = shift;
    return ['.lib_markdup_done'];
}

=head2 extract_intervals_provides

 Title   : extract_intervals_provides
 Usage   : my $provided_files = $obj->extract_intervals_provides('/path/to/lane');
 Function: Find out what files the extract_intervals action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub extract_intervals_provides {
    my $self = shift;
    return $self->{extract_intervals} ? ['.extract_intervals_done'] : ['.lib_markdup_done'];
}

=head2 extract_intervals

 Title   : extract_intervals
 Usage   : $obj->extract_intervals('/path/to/lane', 'lock_filename');
 Function: Extracts intervals at the library level, post markdup.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub extract_intervals {
    my ($self, $lane_path, $action_lock) = @_;
    
    $self->throw("extract_intervals not yet implemented");
    
    my $fofn = $self->{fsu}->catfile($lane_path, '.lib_markdup_done');
    my @files = $self->{io}->parse_fofn($fofn, $lane_path);
    my $verbose = $self->verbose();
    
    my @extract_bams;
    
    # ... makes .strip.markdup.intervals.bam
    
    my $out_fofn = $self->{fsu}->catfile($lane_path, '.extract_intervals_done');
    open(my $ofh, '>', $out_fofn) || $self->throw("Couldn't write to $out_fofn");
    foreach my $out_bam (@extract_bams) {
        print $ofh $out_bam, "\n";
    }
    my $expected = @extract_bams;
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
    return $self->{extract_intervals} ? ['.extract_intervals_done'] : ['.lib_markdup_done'];
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
                              $self->{extract_intervals} ? '.extract_intervals_done' : '.lib_markdup_done',
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
        
        # don't do more than desired merges at once, or we'll kill IO and jobs
        # will fail
        next if $jobs >= $self->{simultaneous_merges};
        
        # skip if we've already handled this group, or have created downstream
        # files from it and deleted the original
        next if $self->{fsu}->file_exists($out_bam);
        my $strip_bam = $out_bam;
        $strip_bam =~ s/\.bam$/.strip.bam/;
        next if $self->{fsu}->file_exists($strip_bam);
        my $markdup_bam = $strip_bam;
        $markdup_bam =~ s/\.bam$/.markdup.bam/;
        next if $self->{fsu}->file_exists($markdup_bam);
        my $intervals_bam = $markdup_bam;
        $intervals_bam =~ s/\.bam$/.intervals.bam/;
        if ($self->{extract_intervals}) {
            next if $self->{fsu}->file_exists($intervals_bam);
        }
        
        # if a higher-level bam already exists, don't repeat making this level
        # bam if we deleted it
        my (undef, $path) = fileparse($out_bam);
        next if $self->_skip_if_higher_level_bam_present($path);
	
        my $this_job_name = $self->{prefix}.$job_name;
        my $pathed_job_name = $self->{fsu}->catfile($path, $this_job_name);
        my $lock_file = $pathed_job_name.'.jids';
        
        # keep simultaneous_merges jobs running all the time
        my $is_running = LSF::is_job_running($lock_file);
        if ($is_running & $LSF::Error) {
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
            next if $jobs > $self->{simultaneous_merges};
            
            $self->archive_bsub_files($path, $this_job_name);
            
            LSF::run($lock_file, $lane_path, $pathed_job_name, $self,
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
    foreach my $job_base (qw(library_merge tag_strip lib_markdup extract_intervals platform_merge sample_merge)) {
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
           $action_name eq 'extract_intervals' || 
           $action_name eq 'tag_strip' || 
           $action_name eq 'library_merge' ||
           $action_name eq 'platform_merge' ||
           $action_name eq 'sample_merge') {
        my $expected_file = $self->{fsu}->catfile($lane_path, ".${action_name}_expected");
        my $not_yet_done = -e $expected_file;
        my $done_file = $self->{fsu}->catfile($lane_path, ".${action_name}_done");
        $self->_merge_check($expected_file, $done_file, $lane_path);
        
        if (($action_name eq 'lib_markdup' || $action_name eq 'extract_intervals' || $action_name eq 'tag_strip') && $not_yet_done && -s $done_file) {
            # we've just completed a non-merge operation; delete the original
            # bams to save space
            my $replace;
            if ($action_name eq 'lib_markdup') {
                $replace = 'markdup';
            }
            elsif ($action_name eq 'extract_intervals') {
                $replace = 'intervals';
            }
            else {
                $replace = 'strip';
            }
            
            my @files = $self->{io}->parse_fofn($done_file);
            foreach my $file (@files) {
                my $previous = $file;
                $previous =~ s/\.$replace.bam$/.bam/;
                if (-s $file && -s $previous) {
                    #unlink($previous);
                    move($previous, "$previous.deleted");
                }
            }
        }
    }
    
    return $self->SUPER::is_finished($lane_path, $action);
}

sub _merge_check {
    my ($self, $expected_file, $done_file, $lane_path) = @_;
    my $fsu = $self->{fsu};
    
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
            if ($fsu->file_exists($_)) {
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

