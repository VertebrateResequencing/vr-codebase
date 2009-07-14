=head1 NAME

VertRes::Pipelines::Mapping - pipeline for mapping fastqs to reference

=head1 SYNOPSIS

# make a file of absolute paths to your desired lane directories, eg:
find $G1K/data -type d -name "*RR*" | grep -v "SOLID" > lanes.fofn

# make a config file:
echo '<lanes.fofn verbose=1;do_cleanup=1' > pipeline.config

# run the pipeline:
run-pipeline -c pipeline.config -v -l mapping_pipeline.log -m VertRes::Pipelines::Mapping

# (and make sure it keeps running by adding that last to a regular cron job)

=head1 DESCRIPTION

A module for carrying out mapping on the Vertebrate Resequencing Informatics 
mapping hierarchy. The hierarchy can consist of multiple different sequencing 
technologies.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Pipelines::Mapping;

use strict;
use warnings;
use VertRes::Utils::Mapping;
use HierarchyUtilities;
use VertRes::IO;
use LSF;

use base qw(VertRes::Pipeline);

our $actions = [ { name     => 'split',
                   action   => \&split,
                   requires => \&split_requires, 
                   provides => \&split_provides },
                { name     => 'map',
                   action   => \&map,
                   requires => \&map_requires, 
                   provides => \&map_provides },
                 { name     => 'merge_and_stat',
                   action   => \&merge_and_stat,
                   requires => \&merge_and_stat_requires, 
                   provides => \&merge_and_stat_provides },
                 { name     => 'cleanup',
                   action   => \&cleanup,
                   requires => \&cleanup_requires, 
                   provides => \&cleanup_provides }];

our %options = (bsub_opts => "-q normal -M5000000 -R 'select[mem>5000] rusage[mem=5000]'",
                sequence_index => '/nfs/sf8/G1K/misc/sequence.index',
                do_cleanup => 0,
                chunk_size => 1000000);

our $split_dir_name = 'split';

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Pipelines::Mapping->new(lane => '/path/to/lane');
 Function: Create a new VertRes::Pipelines::Mapping object;
 Returns : VertRes::Pipelines::Mapping object
 Args    : lane => '/path/to/lane'
           sequence_index => '/path/to/sequence.index' (there is a default)
           do_cleanup => boolean (default false: don't do the cleanup action)
           chunk_size => int (default 1000000)
           other optional args as per VertRes::Pipeline

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(%options, actions => $actions, @args);
    
    # we should have been supplied the option 'lane' which tells us which lane
    # we're in, which lets us choose which mapper module to use.
    my $lane = $self->{lane} || $self->throw("lane directory not supplied, can't continue");
    my $mapping_util = VertRes::Utils::Mapping->new();
    my $mapper_class = $mapping_util->lane_to_module($lane);
    $self->{mapper_class} = $mapper_class;
    
    $self->{io} = VertRes::IO->new;
    
    return $self;
}

=head2 split_requires

 Title   : split_requires
 Usage   : my $required_files = $obj->split_requires('/path/to/lane');
 Function: Find out what files the split action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub split_requires {
    my $self = shift;
    return [$self->_require_fastqs(@_)];
}

sub _require_fastqs {
    my ($self, $lane_path) = @_;
    my @requires;
    
    my @fastq_info = @{HierarchyUtilities::getFastqInfo($lane_path)};
    # *** currently HierarchyUtilities::getFastqInfo tells me what fastq files
    #     are actually there, not which ones are /supposed/ to be there...
    foreach my $info (@fastq_info) {
        push(@requires, $info->[0]);
    }
    
    return @requires;
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
    my @provides = ('.split_complete');
    return \@provides;
}

=head2 split

 Title   : split
 Usage   : $obj->split('/path/to/lane', 'lock_filename');
 Function: Split fastq files into smaller chunks in a subdirectory 'split'.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub split {
    my ($self, $lane_path, $action_lock) = @_;
    
    my @read_args = $self->_get_read_args($lane_path);
    my $mapper_class = $self->{mapper_class};
    my $verbose = $self->verbose;
    my $chunk_size = $self->{chunk_size};
    my $split_dir = $self->{io}->catfile($lane_path, $split_dir_name);
    my $complete_file = $self->{io}->catfile($lane_path, '.split_complete');
    my $script_name = $self->{io}->catfile($lane_path, $self->{prefix}.'split.pl');
    
    # run mapping in an LSF call to a temp script (escaping a perl -e is a bitch)
    open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
    print $scriptfh qq{
use strict;
use $mapper_class;
use VertRes::IO;

my \$mapper = $mapper_class->new(verbose => $verbose);

# split
my \$splits = \$mapper->split_fastq(@read_args,
                                    split_dir => '$split_dir',
                                    chunk_size => $chunk_size);

# (it will only return >0 if all splits created and checked out fine)
\$mapper->throw("split failed - try again?") unless \$splits;

# note how many splits there are in our complete file
my \$io = VertRes::IO->new(file => ">$complete_file");
my \$fh = \$io->fh;
print \$fh \$splits, "\n";
\$io->close;

exit;
    };
    close $scriptfh;
    
    LSF::run($action_lock, $lane_path, $self->{prefix}.'split', $self, qq{perl -w $script_name});
    
    # we've only submitted to LSF, so it won't have finished; we always return
    # that we didn't complete
    return $self->{No};
}

sub _get_read_args {
    my ($self, $lane_path) = @_;
    
    # get the fastq filenames
    # *** currently HierarchyUtilities::getFastqInfo tells me what fastq files
    #     are actually there, not which ones are /supposed/ to be there...
    my @fastq_info = @{HierarchyUtilities::getFastqInfo($lane_path)};
    my @read_args;
    if ($fastq_info[0]->[0]) {
        @read_args = ("read0 => '".$self->{io}->catfile($lane_path, $fastq_info[0]->[0])."'");
    }
    if ($fastq_info[1]->[0]) {
        push(@read_args, ("read1 => '".$self->{io}->catfile($lane_path, $fastq_info[1]->[0])."', ",
                          "read2 => '".$self->{io}->catfile($lane_path, $fastq_info[2]->[0])."'"));
    }
    
    return @read_args;
}

=head2 map_requires

 Title   : map_requires
 Usage   : my $required_files = $obj->map_requires('/path/to/lane');
 Function: Find out what files the map action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub map_requires {
    my ($self, $lane_path) = @_;
    my @requires = $self->_require_fastqs($lane_path);
    
    my %lane_info = %{HierarchyUtilities::lane_info($lane_path)};
    # we require the .bwt so that multiple lanes during 'map' action don't all
    # try and create the .bwt at once automatically if it is missing; could add
    # an index action before map action to solve this...
    push(@requires, '.split_complete', $lane_info{fa_ref}, $lane_info{fa_ref}.'.bwt');
    
    return \@requires;
}

=head2 map_provides

 Title   : map_provides
 Usage   : my $provided_files = $obj->map_provides('/path/to/lane');
 Function: Find out what files the map action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub map_provides {
    my ($self, $lane_path) = @_;
    
    my @provides = ($self->{io}->catfile($lane_path, '.mapping_complete'));
    
    return \@provides;
}

=head2 map

 Title   : map
 Usage   : $obj->map('/path/to/lane', 'lock_filename');
 Function: Carry out mapping of the split fastq files in a lane directory to the
           reference genome. Always generates a bam file ('[splitnum].raw.sorted.bam')
           per split, regardless of the technology/ mapping tool.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub map {
    my ($self, $lane_path, $action_lock) = @_;
    
    # get the reference filename
    my %lane_info = %{HierarchyUtilities::lane_info($lane_path)};
    my $ref_fa = $lane_info{fa_ref} || $self->throw("the reference fasta wasn't known for $lane_path");
    
    my @read_args = $self->_get_read_args($lane_path);
    
    # get the number of splits
    my $num_of_splits = $self->_num_of_splits($lane_path);
    my $split_dir = $self->{io}->catfile($lane_path, $split_dir_name);
    
    # run mapping of each split in LSF calls to temp scripts
    my $mapper_class = $self->{mapper_class};
    my $verbose = $self->verbose;
    my $sequence_index = $self->{sequence_index};
    
    foreach my $split (1..$num_of_splits) {
        my $sam_file = $self->{io}->catfile($lane_path, $split_dir_name, $split.'.raw.sam');
        my $bam_file = $self->{io}->catfile($lane_path, $split_dir_name, $split.'.raw.sorted.bam');
        my $script_name = $self->{io}->catfile($lane_path, $split_dir_name, $self->{prefix}.$split.'.map.pl');
        
        my @split_read_args = ();
        foreach my $read_arg (@read_args) {
            my $split_read_arg = $read_arg;
            $split_read_arg =~ s/\.f[^.]+(?:\.gz)?('(?:, )?)$/.$split.fastq$1/;
            $split_read_arg =~ s/\/([^\/]+)$/\/$split_dir_name\/$1/;
            push(@split_read_args, $split_read_arg);
        }
        
        open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
        print $scriptfh qq{
use strict;
use $mapper_class;

my \$mapper = $mapper_class->new(verbose => $verbose);

# mapping won't get repeated if mapping works the first time but subsequent
# steps fail
my \$ok = \$mapper->do_mapping(ref => '$ref_fa',
                               @split_read_args,
                               output => '$sam_file',
                               insert_size => 2000);

# (it will only return ok and create output if the sam file was created and not
#  truncated)
\$mapper->throw("mapping failed - try again?") unless \$ok;

# add sam header
open(my \$samfh, '$sam_file') || \$mapper->throw("Unable to open sam file '$sam_file'");
my \$head = <\$samfh>;
close(\$samfh);
unless (\$head =~ /^\\\@HD/) {
    my \$ok = \$mapper->add_sam_header('$sam_file', '$sequence_index',
                                       lane_path => '$lane_path');
    \$mapper->throw("Failed to add sam header!") unless \$ok;
}

# convert to mate-fixed sorted bam
unless (-s '$bam_file') {
    my \$ok = \$mapper->sam_to_fixed_sorted_bam('$sam_file', '$bam_file');
    
    unless (\$ok) {
        # (will only return ok and create output bam file if bam was created and
        #  not truncted)
        \$mapper->throw("Failed to create sorted bam file!");
    }
    else {
        unlink('$sam_file');
    }
}

exit;
        };
        close $scriptfh;
        
        LSF::run($action_lock, $lane_path, $self->{prefix}.'map', $self, qq{perl -w $script_name});
    }
    
    # we've only submitted to LSF, so it won't have finished; we always return
    # that we didn't complete
    return $self->{No};
}

sub _num_of_splits {
    my ($self, $lane_path) = @_;
    
    # get the number of splits
    my $split_file = $self->{io}->catfile($lane_path, '.split_complete');
    my $io = VertRes::IO->new(file => $split_file);
    my $split_fh = $io->fh;
    my $splits = <$split_fh>;
    chomp($splits);
    $io->close;
    
    return $splits;
}

=head2 merge_and_stat_requires

 Title   : merge_and_stat_requires
 Usage   : my $required_files = $obj->merge_and_stat_requires('/path/to/lane');
 Function: Find out what files the merge_and_stat action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub merge_and_stat_requires {
    my ($self, $lane_path) = @_;
    
    # get the number of splits
    my $num_of_splits = $self->_num_of_splits($lane_path);
    my $split_dir = $self->{io}->catfile($lane_path, $split_dir_name);
    
    my @requires;
    foreach my $split (1..$num_of_splits) {
        push(@requires, $self->{io}->catfile($lane_path, $split_dir_name, $split.'.raw.sorted.bam'));
    }
    
    return \@requires;
}

=head2 merge_and_stat_provides

 Title   : merge_and_stat_provides
 Usage   : my $provided_files = $obj->merge_and_stat_provides('/path/to/lane');
 Function: Find out what files the merge_and_stat action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub merge_and_stat_provides {
    my ($self, $lane_path) = @_;
    my @provides;
    foreach my $file (qw(raw.sorted.bam rmdup.bam unmapped.bam)) {
        push(@provides, $file);
        
        foreach my $suffix ('bamstat', 'flagstat') {
            push(@provides, $file.'.'.$suffix);
        }
    }
    return \@provides;
}

=head2 merge_and_stat

 Title   : merge_and_stat
 Usage   : $obj->merge_and_stat('/path/to/lane', 'lock_filename');
 Function: Takes the bam files generated during map(), and creates a merged,
           sorted bam file from them. Also removes duplicates and creates
           statistic files. Finally, creates a bam file of unmapped reads, along
           with stat files for that too.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub merge_and_stat {
    my ($self, $lane_path, $action_lock) = @_;
    
    # setup filenames etc. we'll use within our temp script
    my $mapper_class = $self->{mapper_class};
    my $script_name = $self->{io}->catfile($lane_path, $self->{prefix}.'merge_and_stat.pl');
    my $verbose = $self->verbose;
    
    my @bams = @{$self->merge_and_stat_requires($lane_path)};
    
    my $bam_file = $self->{io}->catfile($lane_path, 'raw.sorted.bam');
    my $rmdup_file = $self->{io}->catfile($lane_path, 'rmdup.bam');
    my $unmapped_out = $self->{io}->catfile($lane_path, 'unmapped.bam');
    my @stat_files;
    foreach my $file (qw(raw.sorted.bam rmdup.bam unmapped.bam)) {
        foreach my $suffix ('bamstat', 'flagstat') {
            push(@stat_files, $self->{io}->catfile($lane_path, $file.'.'.$suffix));
        }
    }
    
    # run the multiple steps required for this in an LSF call to a temp script
    open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
    print $scriptfh qq{
use strict;
use $mapper_class;
use VertRes::Wrapper::samtools;

my \$mapper = $mapper_class->new(verbose => $verbose);
my \$samtools = VertRes::Wrapper::samtools->new(verbose => $verbose);

# merge bams
unless (-s '$bam_file') {
    \$samtools->merge_and_check('$bam_file', [qw(@bams)]);
    \$mapper->throw("merging bam failed - try again?") unless \$samtools->run_status == 2;
}

my \$ok = 0;

# rmdup
unless (-s '$rmdup_file') {
    \$ok = \$mapper->rmdup('$bam_file', '$rmdup_file',
                          lane_path => '$lane_path');
    
    unless (\$ok) {
        unlink('$rmdup_file');
        \$mapper->throw("Failed to rmdup the bam file!");
    }
}

# make an unmapped bam
unless (-s '$unmapped_out') {
    \$ok = \$mapper->make_unmapped_bam('$bam_file', '$unmapped_out');
    \$mapper->throw("making unmapped bam failed - try again?") unless \$ok;
}

# make stat files
my \$num_present = 0;
foreach my \$stat_file (qw(@stat_files)) {
    \$num_present++ if -s \$stat_file;
}
unless (\$num_present == ($#stat_files + 1)) {
    my \$ok = \$mapper->stats('$bam_file', '$rmdup_file', '$unmapped_out');
    
    unless (\$ok) {
        foreach my \$stat_file (qw(@stat_files)) {
            unlink(\$stat_file);
        }
        \$mapper->throw("Failed to get stats for the bam files!");
    }
}

exit;
    };
    close $scriptfh;
    
    LSF::run($action_lock, $lane_path, $self->{prefix}.'merge_and_stat', $self, qq{perl -w $script_name});
    
    # we've only submitted to LSF, so it won't have finished; we always return
    # that we didn't complete
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
    my ($self, $lane_path) = @_;
    my @requires = ();
    return \@requires;
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
    my @provides = ();
    return \@provides;
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
    
    my @sais = $self->{io}->get_filepaths($lane_path, suffix => 'sai');
    foreach my $sai (@sais) {
        unlink($sai);
    }
    
    foreach my $file (qw(_log
                         _split.e _split.o _split.pl
                         _map.e _map.o _map.pl
                         _merge_and_stat.e _merge_and_stat.o _merge_and_stat.pl)) {
        unlink($self->{io}->catfile($lane_path, $file));
    }
    
    $self->{io}->rmtree($self->{io}->catfile($lane_path, 'split'));
    
    return $self->{Yes};
}

sub is_finished {
    my ($self, $lane_path, $action) = @_;
    
    # so that we can delete the split dir at the end of the pipeline, and not
    # redo the map action when it sees there are no split bams
    if ($action->{name} eq 'map') {
        my $done_file = $self->{io}->catfile($lane_path, '.mapping_complete');
        
        my $done_bams = 0;
        my $num_of_splits = $self->_num_of_splits($lane_path);
        my $split_dir = $self->{io}->catfile($lane_path, $split_dir_name);
        foreach my $split (1..$num_of_splits) {
            $done_bams += (-s $self->{io}->catfile($lane_path, $split_dir_name, $split.'.raw.sorted.bam')) ? 1 : 0;
        }
        
        if (! -e $done_file && $done_bams == $num_of_splits) {
            system("touch $done_file");
        }
    }
    elsif ($action->{name} eq 'cleanup') {
        return $self->{No};
    }
    
    return $self->SUPER::is_finished($lane_path, $action);
}

1;

