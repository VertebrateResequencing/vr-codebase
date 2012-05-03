=head1 NAME

VertRes::Pipelines::MergeUp - pipeline for creating sample bams from lane bams

=head1 SYNOPSIS

# make a file of absolute paths to the lane bams, eg:
lane_bams.pl --root /abs/META --db g1k_meta --all_samples \
             --slx_mapper bwa --454_mapper ssaha --assembly_name NCBI37 \
             --ignore_platform SOLID --improved --project_regex "low_coverage" \
             > lane_bams.fofn

# make a conf file with root pointing to where you want the merged bams built,
# and setting the lanes, mapper and assembly_name options. Optional settings
# also go here. Your mergeup.conf file may look like:
root    => '/path/to/merge_dir',
module  => 'VertRes::Pipelines::MergeUp',
prefix  => '_',
lane_bams => 'lane_bams.fofn',
previous_merge_root => '/path/to/previous/merge_root',
simultaneous_samples => 200,
data => {
    tag_strip => [qw(OQ XM XG XO)],
    do_cleanup => 1
}

# other options you could add to the mergeup.conf data {} section include:
do_sample_merge => 1 (default is to platform level only)
extract_intervals => {'intervals_file' => 'filename'}
do_index_bams => 1 (default is to not index)

# With the do_index_bams option set, platform bams will be indexed. If 
# do_sample_merge is also set, sample bams will also be indexed.

# if disk space is tight, then one or both of the following options may be set 
# to remove intermediate BAMs. Note that these BAMs may then have to be 
# recreated in a future release.
remove_tag_strip_bams => 1 (default is 0; false)
remove_library_level_bams => 1 (default is 0; false)

# outside the data {} section you can also specify:
previous_merge_root => '/path/to/old/merge_root'
# the previous_merge_root will result in the pipeline checking if a particular
# bam needs to be (re)made depending on if lanes have been removed/added/altered
# since that previous merge. If a bam doesn't need to be made, the new
# merge directory will be a symlink to the old merge directory.

# make another config file that gives the special __MERGEUP__ keyword and the
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
interval extraction for exome work. It can also strip tags from each record
to minimise final bam size.

NB: if you are working on many samples at once, you may run into problems with
too many tmp files being created, taking you over your quota and killing your
jobs. You can help avoid this by supplying the 'memory' option (eg. 12000)
in the data section of your config file.

# If you have manually deleted some files from the command line it may be 
# necessary to run the following command to clear the 
# VertRes::Utils::Filesystem->file_exists cache.
file-exists rm -r /path/to/merge_dir
# Likewise for previous merge directories, if linking to them.

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
                { name     => 'tag_strip',
                  action   => \&tag_strip,
                  requires => \&tag_strip_requires, 
                  provides => \&tag_strip_provides },
                { name     => 'library_merge',
                  action   => \&library_merge,
                  requires => \&library_merge_requires, 
                  provides => \&library_merge_provides },
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
                { name     => 'index_platform_bams',
                  action   => \&index_platform_bams,
                  requires => \&index_platform_bams_requires, 
                  provides => \&index_platform_bams_provides },
                { name     => 'index_sample_bam',
                  action   => \&index_sample_bam,
                  requires => \&index_sample_bam_requires, 
                  provides => \&index_sample_bam_provides },
                { name     => 'cleanup',
                  action   => \&cleanup,
                  requires => \&cleanup_requires, 
                  provides => \&cleanup_provides }];

our %options = (do_cleanup => 0,
                do_sample_merge => 0,
                bsub_opts => '',
                tag_strip => [],
                remove_tag_strip_bams => 0,
                remove_library_level_bams => 0,
                do_index_bams => 0,
                dont_wait => 1);

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Pipelines::MergeUp->new(lane => '/path/to/lane');
 Function: Create a new VertRes::Pipelines::MergeUp object.
 Returns : VertRes::Pipelines::MergeUp object
 Args    : hierarchy_root => '/abs/path' (REQUIRED; the root of the hierarchy)
           lane_bams => fofn (REQUIRED; a file listing all the bams you want
                              to merge)
           
           tag_strip => [qw(OQ XM XG XO)] (default do not strip any tags; supply
                                           an array ref of tags to remove if
                                           desired)
           remove_tag_strip_bams => boolean (default false: removes tag strip bams 
                                 after a library merge unless merged bams are 
                                 symlinks to the tag strip bams)
           remove_library_level_bams => boolean (default false: removes library 
                                 level bams after a platform merge unless merged 
                                 bams are symlinks to the library level bams)
           do_cleanup => boolean (default false: don't do the cleanup action)
           do_index_bams => boolean (default false: don't index bams)
           do_sample_merge => boolean (default false: don't create sample-level
                                       bams)
           extract_intervals => {intervals_file => 'filename'}
                                    (after marking duplicates, extract reads for
                                    only these intervals; default keep all
                                    reads.  if this option used,
                                    intervals_file is REQUIRED)
           memory => int (a default and minimum memory applies to each action;
                          if this is supplied and higher than the minimum it
                          will be used)
           other optional args as per VertRes::Pipeline

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(%options, actions => $actions, @args);
    
    $self->{hierarchy_root} || $self->throw("hierarchy_root path not supplied, can't continue");
    $self->{lane_bams} || $self->throw("lane_bams fofn not supplied, can't continue");
    
    if ($self->{extract_intervals} && !($self->{extract_intervals}->{intervals_file})) {
        $self->throw("extract_intervals must have intervals_file supplied, can't continue");
    }
    
    $self->{io} = VertRes::IO->new;
    $self->{fsu} = VertRes::Utils::FileSystem->new;
    
    $self->{do_tag_strip} = 0;
    if (defined $self->{tag_strip} && ref($self->{tag_strip}) eq 'ARRAY' && @{$self->{tag_strip}} > 0) {
        $self->{do_tag_strip} = 1;
    }
    unless ( $self->{do_tag_strip} ) { $self->{remove_tag_strip_bams} = 0 };
    
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

=head2 tag_strip_requires

 Title   : tag_strip_requires
 Usage   : my $required_files = $obj->tag_strip_requires('/path/to/lane');
 Function: Find out what files the tag_strip action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub tag_strip_requires {
    my $self = shift;
    return ['.hierarchy_made'];
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
    return $self->{do_tag_strip} ? ['.tag_strip_done'] : ['.hierarchy_made'];
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
    
    my $fofn = $self->{fsu}->catfile($lane_path, '.hierarchy_made');
    my @files = $self->{io}->parse_fofn($fofn, $lane_path);
    my $verbose = $self->verbose();
    
    my @strip_bams;
    foreach my $lane_bam (@files) {
        $lane_bam = $self->{fsu}->catfile($lane_path, $lane_bam);
        my($basename, $path) = fileparse($lane_bam);
        
        my $strip_bam = $lane_bam;
        $strip_bam =~ s/\.bam$/.strip.bam/;
        push(@strip_bams, $strip_bam);
        next if $self->{fsu}->file_exists($strip_bam);
        
        # if a higher-level bam already exists, don't repeat making this level
        # bam if we deleted it
        next if $self->_downstream_bam_present($strip_bam);
        next unless $self->{fsu}->file_exists($lane_bam);
        
        my $job_name = $self->{prefix}.'tag_strip';
        $self->archive_bsub_files($path, $job_name);
        $job_name = $self->{fsu}->catfile($path, $job_name);
        
        LSF::run($action_lock, $lane_path, $job_name, $self,
                 qq{perl -MVertRes::Utils::Sam -Mstrict -e "VertRes::Utils::Sam->new(verbose => $verbose)->tag_strip(qq[$lane_bam], qq[$strip_bam], qw(@{$self->{tag_strip}})) || die qq[tag_strip failed for $lane_bam\n];"});
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

sub _downstream_bam_present {
    my ($self, $bam) = @_;
    
	my ($basename, $path) = fileparse($bam);
	my $rel_path = File::Spec->abs2rel( $path, $self->{hierarchy_root} );
	$rel_path =~ s/\/$//;
	
    my @dirs = File::Spec->splitdir($rel_path);
	my @levels = ( 'root', 'pop', 'sample', 'platform', 'library', 'lane' );
    my $level = $levels[scalar @dirs];
	
	if ($level eq 'lane') {
		my $strip_bam = $bam;
        $strip_bam =~ s/\.bam$/.strip.bam/;
        return 1 if $self->{fsu}->file_exists($strip_bam);
		
		# No downstream bams exist at this level: move up one level.
		$path = dirname($path);
		$level = 'library';
	}
	
	if ($level eq 'library') {
		my $merged_bam = $bam;
		$merged_bam =~ s/lane\.//;
		$merged_bam =~ s/strip\.//;
        return 1 if ( $self->{fsu}->file_exists($merged_bam) && ! -l $merged_bam );
		
		my $markdup_bam = $merged_bam;
        $markdup_bam =~ s/\.bam$/.markdup.bam/;
        return 1 if $self->{fsu}->file_exists($markdup_bam);
        
        if ($self->{extract_intervals}) {
			my $intervals_bam = $markdup_bam;
	        $intervals_bam =~ s/\.bam$/.intervals.bam/;
            return 1 if $self->{fsu}->file_exists($intervals_bam);
        }
        
		# No downstream bams exist at this level: move up one level.
		$path = dirname($path);
		$level = 'platform';
	} 
	
	if ($level eq 'platform') {
		my $platform_bam = File::Spec->catfile($path, 'raw.bam');
		return 1 if ( $self->{fsu}->file_exists($platform_bam) && ! -l $platform_bam );

		# No downstream bams exist at this level: move up one level.
		$path = dirname($path);
		$level = 'sample';
	}
	
	if ($level eq 'sample') {
		my $sample_bam = File::Spec->catfile($path, 'raw.bam');
		return 1 if ( $self->{fsu}->file_exists($sample_bam) && ! -l $sample_bam );
	}
	
    return 0;
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
    return $self->{do_tag_strip} ? ['.tag_strip_done'] : ['.hierarchy_made'];
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
                              $self->{do_tag_strip} ? '.tag_strip_done' : '.hierarchy_made',
                              undef,
                              'library_merge',
                              '.library_merge_expected',
                              'long');
    
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
    
    my $memory = $self->{memory};
    my $java_mem;
    if (! defined $memory || $memory < 6100) {
        $memory = 6100;
        $java_mem = 5000;
    }
    $java_mem ||= int($memory * 0.9);
    my $queue = $memory >= 16000 ? "hugemem" : "long";
    
    my $orig_bsub_opts = $self->{bsub_opts};
    $self->{bsub_opts} = "-q $queue -M${memory}000 -R 'select[mem>$memory] rusage[mem=$memory]'";
    
    my @markdup_bams;
    foreach my $merge_bam (@files) {
        $merge_bam = $self->{fsu}->catfile($lane_path, $merge_bam);
        my($basename, $path) = fileparse($merge_bam);
        
        my $markdup_bam = $merge_bam;
        $markdup_bam =~ s/\.bam$/.markdup.bam/;
        push(@markdup_bams, $markdup_bam);
        next if $self->{fsu}->file_exists($markdup_bam);
        
        # if a higher-level bam already exists, don't repeat making this level
        # bam if we deleted it
        next if $self->_downstream_bam_present($markdup_bam);
        next unless $self->{fsu}->file_exists($merge_bam);
        
        my $job_name = $self->{prefix}.'lib_markdup_'.$basename;
        $self->archive_bsub_files($path, $job_name);
        $job_name = $self->{fsu}->catfile($path, $job_name);
        
        LSF::run($action_lock, $lane_path, $job_name, $self,
                 qq{perl -MVertRes::Utils::Sam -Mstrict -e "VertRes::Utils::Sam->new(verbose => $verbose)->markdup(qq[$merge_bam], qq[$markdup_bam], java_memory => $java_mem) || die qq[markdup failed for $merge_bam\n];"});
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
    
    my $fofn = $self->{fsu}->catfile($lane_path, '.lib_markdup_done');
    my @files = $self->{io}->parse_fofn($fofn, $lane_path);
    my $verbose = $self->verbose();
    
    my @extract_bams;
    foreach my $markdup_bam (@files){
        $markdup_bam = $self->{fsu}->catfile($lane_path, $markdup_bam); 
        my($basename, $path) = fileparse($markdup_bam);
        
        my $extract_bam = $markdup_bam;
        $extract_bam =~ s/\.bam$/.intervals.bam/;
        push(@extract_bams, $extract_bam);
        next if $self->{fsu}->file_exists($extract_bam);
        
        # if a higher-level bam already exists, don't repeat making this level
        # bam if we deleted it
        next if $self->_downstream_bam_present($extract_bam);
        next unless $self->{fsu}->file_exists($markdup_bam);

        my $job_name = $self->{prefix}.'extract_intervals_'.$basename;
        $self->archive_bsub_files($path, $job_name);
        $job_name = $self->{fsu}->catfile($path, $job_name);
        LSF::run($action_lock, $lane_path, $job_name, $self,
                 qq~perl -MVertRes::Utils::Sam -Mstrict -e "VertRes::Utils::Sam->new(verbose => $verbose)->extract_intervals_from_bam(qq[$markdup_bam], qq[$self->{extract_intervals}->{intervals_file}], qq[$extract_bam]) || die qq[extract_intervals failed for $markdup_bam\n];"~);
    }
    
    my $out_fofn = $self->{fsu}->catfile($lane_path, '.extract_intervals_expected');
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
    return ['.platform_merge_done'];
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
    
    my $memory = $self->{memory};
    if (! defined $memory || $memory < 5556) {
        $memory = 5556;
    }
    my $java_mem = int($memory * 0.9);
    $queue ||= 'normal';
    $queue = $memory >= 16000 ? "hugemem" : $queue;
    
    my $orig_bsub_opts = $self->{bsub_opts};
    $self->{bsub_opts} = "-q $queue -M${memory}000 -R 'select[mem>$memory] rusage[mem=$memory]'";
    
    my $group_by_basename = ! defined $output_basename;
    my %grouped_bams = $self->_fofn_to_bam_groups($lane_path, $fofn, $group_by_basename);
    
    my $verbose = $self->verbose;
    
    my @out_bams;
    my %lane_paths;
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
        
		my @rel_path = File::Spec->splitdir($out_dir);
		shift @rel_path;
        my $current_path = $self->{fsu}->catfile($lane_path, @rel_path);
        my $out_bam = $self->{fsu}->catfile($current_path, $out_bam_name);
        push(@out_bams, $out_bam);
        
        # skip if we've already handled this group, or have created downstream
        # files from it and deleted the original
        next if $self->{fsu}->file_exists($out_bam);
        next if $self->_downstream_bam_present($out_bam);
		
		my (undef, $path) = fileparse($out_bam);
        my $this_job_name = $self->{prefix}.$job_name.'_'.$out_bam_name;
        my $pathed_job_name = $self->{fsu}->catfile($path, $this_job_name);
        my $lock_file = $pathed_job_name.'.jids';
        
        # keep all the jobs running all the time
        my $is_running = LSF::is_job_running($lock_file);
        if ($is_running & $LSF::Error) {
            next;
        }
        elsif ($is_running & $LSF::Running) {
            next;
        }
        elsif ($is_running & $LSF::Done) {
            next;
        }
        else {
            $self->archive_bsub_files($path, $this_job_name);
            
            LSF::run($lock_file, $lane_path, $pathed_job_name, $self,
                 qq{perl -MVertRes::Utils::Sam -Mstrict -e "VertRes::Utils::Sam->new(verbose => $verbose, java_memory => $java_mem)->merge(qq[$out_bam], qw(@bams)) || die qq[merge failed for (@bams) -> $out_bam\n];"});
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
        my $parent = $self->{fsu}->catfile('SAMPLE_ROOT', @parts);
        
        my $group;
        if ($group_by_basename) {
            $basename =~ /^([^.]+).*\.bam$/ || $self->throw("Got a bad bamfile name in the fofn: $_");
            $group = "$parent:!:$1"; # :!: is a separator token
        }
        else {
            $group = $parent;
        }
        
        push(@{$bams_groups{$group}}, $self->{fsu}->catfile($lane_path, $path));
    }
    
    return %bams_groups;
}

=head2 index_platform_bams_requires

 Title   : index_platform_bams_requires
 Usage   : my $required_files = $obj->index_platform_bams_requires('/path/to/lane');
 Function: Find out what files the index_platform_bams action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub index_platform_bams_requires {
    my ($self, $lane_path) = @_;
	return $self->{do_index_bams} ? ['.platform_merge_done'] : [];
}

=head2 index_platform_bams_provides

 Title   : index_platform_bams_provides
 Usage   : my $provided_files = $obj->index_platform_bams_provides('/path/to/lane');
 Function: Find out what files the index_platform_bams action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub index_platform_bams_provides {
    my ($self, $lane_path) = @_;
	return $self->{do_index_bams} ? ['.index_platform_bams_done'] : [];
}

=head2 index_platform_bams

 Title   : index_platform_bams
 Usage   : $obj->index_platform_bams('/path/to/lane', 'lock_filename');
 Function: Index platform bams.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub index_platform_bams {
    my ($self, $lane_path, $action_lock) = @_;
	return $self->{Yes} unless $self->{do_index_bams};
    
    $self->_index_bams($lane_path,
                              $action_lock,
                              '.platform_merge_done',
                              'index_platform_bams',
                              '.index_platform_bams_expected');
    
    return $self->{No};
}


=head2 index_sample_bam_requires

 Title   : index_sample_bam_requires
 Usage   : my $required_files = $obj->index_sample_bam_requires('/path/to/lane');
 Function: Find out what files the index_sample_bam action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub index_sample_bam_requires {
    my ($self, $lane_path) = @_;
	return ( $self->{do_index_bams} && $self->{do_sample_merge} ) ? ['.index_platform_bams_done', '.sample_merge_done'] : [];
}

=head2 index_sample_bam_provides

 Title   : index_sample_bam_provides
 Usage   : my $provided_files = $obj->index_sample_bam_provides('/path/to/lane');
 Function: Find out what files the index_sample_bam action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub index_sample_bam_provides {
    my ($self, $lane_path) = @_;
	return ( $self->{do_index_bams} && $self->{do_sample_merge} ) ? ['.index_sample_bam_done'] : [];
}

=head2 index_sample_bam

 Title   : index_sample_bam
 Usage   : $obj->index_sample_bam('/path/to/lane', 'lock_filename');
 Function: Index sample bam.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub index_sample_bam {
    my ($self, $lane_path, $action_lock) = @_;
    return $self->{Yes} unless ( $self->{do_index_bams} && $self->{do_sample_merge} );
    
    $self->_index_bams($lane_path,
                              $action_lock,
                              '.sample_merge_done',
                              'index_sample_bam',
                              '.index_sample_bam_expected');
    
    return $self->{No};
}


sub _index_bams {
	my ($self, $lane_path, $action_lock, $fofn, $job_name, $out_fofn, $queue) = @_;
	
	my $this_fofn = $self->{fsu}->catfile($lane_path, $fofn);
    my @bams = $self->{io}->parse_fofn($this_fofn, '/');
	
	my $this_out_fofn = $self->{fsu}->catfile($lane_path, $out_fofn);
	
	my @bams_to_index;
	open(my $ofh, '>', $this_out_fofn) || $self->throw("Couldn't write to $this_out_fofn");
	foreach my $bam (@bams) {
		my $bai = "$bam.bai";
		if (-l $bam) {
			my $sym_bai = readlink($bam) . ".bai";
			my $path = dirname $bai;
			my $sym_bai_path = $self->{fsu}->catfile($path, $sym_bai);
			symlink $sym_bai, $bai if (-s $sym_bai_path); 
		}
	    print $ofh "$bai\n";
		push @bams_to_index, $bam unless (-s $bai);
	}
	my $expected = @bams;
	print $ofh "# expecting $expected\n";
	close($ofh);
	
	return unless @bams_to_index;
	
	my $verbose = $self->verbose();
    my $this_job_name = $self->{prefix}.$job_name;
    $self->archive_bsub_files($lane_path, $this_job_name);
    $job_name = $self->{fsu}->catfile($lane_path, $this_job_name);
	LSF::run($action_lock, $lane_path, $job_name, $self,
qq~perl -MVertRes::Utils::Sam -Mstrict -e "my \@bams = qw(@bams_to_index); VertRes::Utils::Sam->new(verbose => $verbose)->index_bams(files=>\\\@bams) || die qq[index_bams failed for $this_fofn\n];"~);

    return;
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
 Function: Unlink all the pipeline-related files (_*) in a lane.
           NB: do_cleanup => 1 must have been supplied to new();
 Returns : $VertRes::Pipeline::Yes
 Args    : lane path, name of lock file to use

=cut

sub cleanup {
    my ($self, $lane_path, $action_lock) = @_;
    return $self->{Yes} unless $self->{do_cleanup};
    
    my $prefix = $self->{prefix};
    my $fsu = $self->{fsu};
    
    foreach my $file (qw(log job_status)) {
        unlink($fsu->catfile($lane_path, $prefix.$file));
    }
    
	my @suffixes = qw(.o .e .jids .o.previous .e.previous);
	my @basenames = ('', '_raw.bam', '_pe.bam', '_se.bam', '_pe.markdup.bam', '_se.markdup.bam');
	my @jobs = qw(tag_strip library_merge lib_markdup extract_intervals platform_merge 
					sample_merge index_platform_bams index_sample_bam);
	
	my @dirs = ($lane_path);
    foreach my $job_base (@jobs) {
		my $done_file = $fsu->catfile($lane_path, ".$job_base\_done");
		next unless -e $done_file;
		my @files = $self->{io}->parse_fofn($done_file, '/');
		push @dirs, map( dirname($_), @files );
	}
	my %seen = ();
	my @uniq_dirs = grep { ! $seen{$_} ++ } @dirs;
	
	foreach my $job_base (@jobs) {
		foreach my $path (@uniq_dirs) {
		    my $file_base = $fsu->catfile($path, $prefix);
			foreach my $basename (@basenames) {
		        foreach my $suffix (@suffixes) {
					my $file = "$file_base$job_base$basename$suffix";
		            unlink($file) if -e $file;
		        }
			}
		}
    }
	
	if ( $self->{do_tag_strip} && $self->{remove_tag_strip_bams}) 
	{
		my $tag_done_file = $fsu->catfile($lane_path, ".tag_strip_done");
		my @files = $self->{io}->parse_fofn($tag_done_file);
		foreach my $tag_strip_bam (@files) {
			next unless $fsu->file_exists($tag_strip_bam);
			unlink $tag_strip_bam;
			$fsu->file_exists($tag_strip_bam, wipe_out => 1);
		}
	} 
	
	if ( $self->{remove_library_level_bams}) 
	{
		my $markdup_done_file = $fsu->catfile($lane_path, ".lib_markdup_done");
		my $extract_intervals_done_file = $fsu->catfile($lane_path, ".extract_intervals_done");
		my $done_file = $self->{extract_intervals} ? $extract_intervals_done_file : $markdup_done_file;
		my @files = $self->{io}->parse_fofn($done_file);
		foreach my $bam (@files) {
			next unless $fsu->file_exists($bam);
			my $path = dirname($bam);
			$path = dirname($path);
			my $platform_bam = File::Spec->catfile($path, 'raw.bam');
			if ($fsu->file_exists($platform_bam) && ! -l $platform_bam) {
				unlink $bam;
				$fsu->file_exists($bam, wipe_out => 1);
			}
		}
	}

    return $self->{Yes};
}

sub is_finished {
    my ($self, $lane_path, $action) = @_;
    
    my $action_name = $action->{name};
    my $fsu = $self->{fsu};
    
    if ($action_name eq 'cleanup') {
        return $self->{No};
    } elsif ($action_name eq 'lib_markdup' ||
           $action_name eq 'extract_intervals' || 
           $action_name eq 'tag_strip' || 
           $action_name eq 'library_merge' ||
           $action_name eq 'platform_merge' ||
           $action_name eq 'sample_merge' ||
           $action_name eq 'index_platform_bams' ||
           $action_name eq 'index_sample_bam') {
        my $expected_file = $fsu->catfile($lane_path, ".${action_name}_expected");
        my $not_yet_done = -e $expected_file;
        my $done_file = $fsu->catfile($lane_path, ".${action_name}_done");
        $self->_merge_check($expected_file, $done_file, $lane_path);
        
        if ($not_yet_done && -s $done_file) {
            if ($action_name eq 'lib_markdup' || $action_name eq 'extract_intervals') {
                # we've just completed a non-merge library operation; delete the
                # previous bam in the chain
                my $replace;
                if ($action_name eq 'lib_markdup') {
                    $replace = 'markdup';
                }
                elsif ($action_name eq 'extract_intervals') {
                    $replace = 'intervals';
                }
                
                my @files = $self->{io}->parse_fofn($done_file);
                foreach my $file (@files) {
					next unless ($fsu->file_exists($file));
                    my $previous = $file;
                    $previous =~ s/\.$replace.bam$/.bam/;
                    if ($fsu->file_exists($previous)) {
                        unlink($previous);
						$fsu->file_exists($previous, wipe_out => 1);
                    }
                }
            	
				if ( $self->{remove_tag_strip_bams} ) {
					# at this point it is now safe to unlink the tag_strip bams
					my $tag_done_file = $fsu->catfile($lane_path, ".tag_strip_done");
					my @files = $self->{io}->parse_fofn($tag_done_file);
					foreach my $tag_strip_bam (@files) {
						next unless $fsu->file_exists($tag_strip_bam);
						unlink $tag_strip_bam;
						$fsu->file_exists($tag_strip_bam, wipe_out => 1);
					}
				}
			} 
			elsif ($action_name eq 'platform_merge' && $self->{remove_library_level_bams}) {
				# we've just completed the platform merge; delete the library level
				# bams unless the platform level bams are just symlinks to those
				my $markdup_done_file = $fsu->catfile($lane_path, ".lib_markdup_done");
				my $extract_intervals_done_file = $fsu->catfile($lane_path, ".extract_intervals_done");
				my $done_file = $self->{extract_intervals} ? $extract_intervals_done_file : $markdup_done_file;
				my @files = $self->{io}->parse_fofn($done_file);
				foreach my $bam (@files) {
					next unless $fsu->file_exists($bam);
					my $path = dirname($bam);
					$path = dirname($path);
					my $platform_bam = File::Spec->catfile($path, 'raw.bam');
					if ($fsu->file_exists($platform_bam) && ! -l $platform_bam) {
						unlink $bam;
						$fsu->file_exists($bam, wipe_out => 1);
					}
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
			
			# If the file exists or one of the downstream bams exist, we can treat this 
			# level as complete
            if ( $fsu->file_exists($_) || $self->_downstream_bam_present($_) ) {
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

