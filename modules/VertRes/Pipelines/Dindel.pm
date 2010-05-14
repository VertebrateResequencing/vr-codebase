=head1 NAME

VertRes::Pipelines::Dindel - pipeline for running Dindel on bams

=head1 SYNOPSIS

# make a file of absolute paths to your release bam-containing directories:
find /path/to/REL -maxdepth 3 -mindepth 3 -type d > dindel.fofn
# (each directoy should contain 'release.bam' and 'release.bam.bai' files)

# make a conf file with root pointing to the root of the release directory.
# Optional settings also go here. Your dindel.conf file may look like:
root    => '/path/to/REL',
module  => 'VertRes::Pipelines::Dindel',
prefix  => '_',
data => {
    ref => '/path/to/ref.fa',
    simultaneous_jobs => 100
}

# make a pipeline file that associates the lanes with the conf file:
echo "<dindel.fofn dindel.conf" > dindel.pipeline

# run the pipeline:
run-pipeline -c dindel.pipeline -v

# (and make sure it keeps running by adding that last to a regular cron job)

=head1 DESCRIPTION

A module for running Dindel on a set of release bams. Dindel output appears in
a dindel subdir of each dir containing a release bam.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Pipelines::Dindel;

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

our $actions = [{ name     => 'get_variants',
                  action   => \&get_variants,
                  requires => \&get_variants_requires, 
                  provides => \&get_variants_provides },
                { name     => 'split_variants',
                  action   => \&split_variants,
                  requires => \&split_variants_requires, 
                  provides => \&split_variants_provides },
                { name     => 'call',
                  action   => \&call,
                  requires => \&call_requires, 
                  provides => \&call_provides },
                { name     => 'merge',
                  action   => \&merge,
                  requires => \&merge_requires, 
                  provides => \&merge_provides }];
                
my $dindel_base = '/lustre/scratch102/user/sb10/dindel/';
our %options = (simultaneous_jobs => 100,
                vars_per_file => 1000,
                files_per_block => 200,
                dindelPars => '--maxRead 4000 --pMut 1e-6 --priorSNP 0.001 --priorIndel 0.0001 --maxLengthIndel 5 --slower --outputGLF --doDiploid',
                addDindelOpt => '',
                dindel_scripts => $dindel_base.'PythonCode/dindel',
                dindel_bin => $dindel_base.'dindel',
                bsub_opts => '');

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Pipelines::Dindel->new(lane => '/path/to/lane');
 Function: Create a new VertRes::Pipelines::Dindel object.
 Returns : VertRes::Pipelines::Dindel object
 Args    : lane_path => '/path/to/dir/containing/release.bam' (REQUIRED, set by
                         run-pipeline automatically)
           ref => '/path/to/ref.fa' (REQUIRED)

           chrom => string (to only do 1 chr, provide the chromosome name, eg.
                            20)
           simultaneous_jobs => int (default 200; the number of jobs to
                                     do at once - limited to avoid IO
                                     problems)
           other optional args as per VertRes::Pipeline

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(%options, actions => $actions, @args);
    
    $self->{lane_path} || $self->throw("lane_path (misnomer; actually releasese bam) directory not supplied, can't continue");
    $self->{ref} || $self->throw("ref not supplied, can't continue");
    -x $self->{dindel_bin} || $self->throw("bad dindel_bin $self->{dindel_bin}");
    -d $self->{dindel_scripts} || $self->throw("bad dindel_scripts $self->{dindel_scripts}");
    
    $self->{io} = VertRes::IO->new;
    $self->{fsu} = VertRes::Utils::FileSystem->new;
    
    $self->{outdir_basename} = 'dindel';
    if ($self->{chrom}) {
        my $chrom_bam = `ls $self->{lane_path}/*.chrom$self->{chrom}.*.bam`;
        chomp($chrom_bam);
        if ($chrom_bam) {
            $self->{bam_file_name} = basename($chrom_bam);
            $self->{outdir_basename} .= '_chrom'.$self->{chrom};
        }
        else {
            $self->throw("Could not find a chrom$self->{chrom} bam in $self->{lane_path}");
        }
    }
    else {
        $self->{bam_file_name} = 'release.bam';
    }
    
    # setup output and log dirs
    my $out_dir = $self->{fsu}->catfile($self->{lane_path}, $self->{outdir_basename});
    $self->{outdir} = $out_dir;
    unless (-d $out_dir) {
        mkdir($out_dir) || $self->throw("Could not create directory '$out_dir'");
    }
    
    my $logs_dir = $self->{fsu}->catfile($out_dir, 'logs');
    $self->{logs} = $logs_dir;
    unless (-d $logs_dir) {
        mkdir($logs_dir) || $self->throw("Could not create directory '$logs_dir'");
    }
    
    $self->{bamfiles_fofn} = $self->{fsu}->catfile($out_dir, 'bamfiles.fofn');
    unless (-s $self->{bamfiles_fofn}) {
        open(my $bfh, '>', $self->{bamfiles_fofn}) || $self->throw("Could not write to $self->{bamfiles_fofn}");
        my $release_bam = $self->{fsu}->catfile($self->{lane_path}, $self->{bam_file_name});
        print $bfh $release_bam, "\n";
        close($bfh);
    }
    
    return $self;
}

=head2 get_variants_requires

 Title   : get_variants_requires
 Usage   : my $required_files = $obj->get_variants_requires('/path/to/lane');
 Function: Find out what files the get_variants action needs before
           it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub get_variants_requires {
    my $self = shift;
    return [$self->{ref}, $self->{bam_file_name}, $self->{bam_file_name}.'.bai'];
}

=head2 get_variants_provides

 Title   : get_variants_provides
 Usage   : my $provided_files = $obj->get_variants_provides('/path/to/lane');
 Function: Find out what files the get_variants action generates on
           success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub get_variants_provides {
    my $self = shift;
    return [$self->{fsu}->catfile($self->{outdir_basename}, 'variants', 'variants.variants.txt'),
            $self->{fsu}->catfile($self->{outdir_basename}, 'variants', 'variants.libraries.txt')];
}

=head2 get_variants

 Title   : get_variants
 Usage   : $obj->get_variants('/path/to/lane', 'lock_filename');
 Function: Gets variants from the bam file.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub get_variants {
    my ($self, $lane_path, $action_lock) = @_;
    
    my $out_dir = $self->{fsu}->catfile($self->{outdir}, 'variants');
    unless (-d $out_dir) {
        mkdir($out_dir) || $self->throw("Could not create directory '$out_dir'");
    }
    
    my $bam_file = $self->{fsu}->catfile($lane_path, $self->{bam_file_name});
    my $out_file = $self->{fsu}->catfile($out_dir, 'running.variants');
    
    my $job_name = $self->{fsu}->catfile($self->{logs}, 'get_variants');
    $self->archive_bsub_files($self->{logs}, 'get_variants');
    
    LSF::run($action_lock, $lane_path, $job_name, $self,
             qq{$self->{dindel_bin} --analysis getCIGARindels --bamFile $bam_file --ref $self->{ref} --outputFile $out_file});
    
    return $self->{No};
}

=head2 split_variants_requires

 Title   : split_variants_requires
 Usage   : my $required_files = $obj->split_variants_requires('/path/to/lane');
 Function: Find out what files the split_variants action needs before
           it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub split_variants_requires {
    my $self = shift;
    return $self->get_variants_provides;
}

=head2 split_variants_provides

 Title   : split_variants_provides
 Usage   : my $provided_files = $obj->split_variants_provides('/path/to/lane');
 Function: Find out what files the split_variants action generates on
           success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub split_variants_provides {
    my ($self, $lane_path) = @_;
    
    my $split_dir = $self->{fsu}->catfile($self->{outdir_basename}, 'variants', 'split');
    my $abs_split_dir = $self->{fsu}->catfile($lane_path, $split_dir);
    my @provides;
    
    # don't know at the moment how many splits we're supposed to have, so cheat
    # and see how many we do have and assume that's correct!
    if (-d $abs_split_dir) {
        opendir(my $splitfh, $abs_split_dir) || $self->throw("Could not open dir $abs_split_dir");
        foreach my $file (readdir($splitfh)) {
            if ($file =~ /\.variants\.txt$/) {
                push(@provides, $self->{fsu}->catfile($split_dir, $file));
            }
        }
        @provides > 1 || $self->throw("Something went wrong; not enough variant files in @provides");
    }
    else {
        @provides = ($split_dir);
    }
    
    return \@provides;
}

=head2 split_variants

 Title   : split_variants
 Usage   : $obj->split_variants('/path/to/lane', 'lock_filename');
 Function: Split variants.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub split_variants {
    my ($self, $lane_path, $action_lock) = @_;
    
    # split the variants per chromosome
    # (group variants into windows and split into individual files)
    my $split_dir = $self->{fsu}->catfile($self->{outdir}, 'variants', 'split');
    unless (-d $split_dir) {
        mkdir($split_dir) || $self->throw("Could not create directory '$split_dir'");
    }
    
    my $job_name = $self->{fsu}->catfile($self->{logs}, 'split_variants');
    $self->archive_bsub_files($self->{logs}, 'split_variants');
    
    my $var_file = $self->{fsu}->catfile($self->{outdir}, 'variants', 'variants.variants.txt');
    
    #*** we need a way of checking if it completed successfully; need
    #    to know how many files it is supposed to create in $split_dir
    my $orig_bsub_opts = $self->{bsub_opts};
    $self->{bsub_opts} = ' -M4000000 -R \'select[mem>4000] rusage[mem=4000]\'';

    LSF::run($action_lock, $lane_path, $job_name, $self,
             qq{cat $var_file | python $self->{dindel_scripts}/MergeCandidates.py --outputDir $split_dir --minDist 25 --varPerFile $self->{vars_per_file}});
    
    $self->{bsub_opts} = $orig_bsub_opts;
    
    return $self->{No};
}

=head2 call_requires

 Title   : call_requires
 Usage   : my $required_files = $obj->call_requires('/path/to/lane');
 Function: Find out what files the call action needs before
           it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub call_requires {
    my ($self, $lane_path) = @_;
    return $self->split_variants_provides($lane_path);
}

=head2 call_provides

 Title   : call_provides
 Usage   : my $provided_files = $obj->call_provides('/path/to/lane');
 Function: Find out what files the call action generates on
           success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub call_provides {
    my ($self, $lane_path) = @_;
    
    my @provides = ($self->{fsu}->catfile($self->{outdir_basename}, 'libraries', 'libraries.txt'));
    
    my @split_var_files = @{$self->split_variants_provides($lane_path)};
    my $splits = 0;
    my $blocks = 0;
    my $blocks_dir = $self->{fsu}->catfile($self->{outdir_basename}, 'blocks');
    my $block_dir;
    foreach my $var_file (@split_var_files) {
        $splits++;
        if ($splits % $self->{files_per_block} == 1) {
            $blocks++;
            $block_dir = $self->{fsu}->catfile($blocks_dir, 'block.'.$blocks);
        }
        
        foreach my $type ('calls', 'glf', 'log') {
            push(@provides, $self->{fsu}->catfile($block_dir, "dindel.$splits.$type.txt.gz"));
        }
    }
    
    return \@provides;
}

=head2 call

 Title   : call
 Usage   : $obj->call('/path/to/lane', 'lock_filename');
 Function: Calls variants in blocks.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub call {
    my ($self, $lane_path, $action_lock) = @_;
    
    # merge libraries
    #*** not necessary? We only ever run on one bam at a time which is already
    #    a merge of multiple libraries, so $lib_out_file and $lib_file end up
    #    being the same thing
    my $lib_dir = $self->{fsu}->catfile($self->{outdir}, 'libraries');
    unless (-d $lib_dir) {
        mkdir($lib_dir) || $self->throw("Could not create directory '$lib_dir'");
    }
    my $lib_out_file = $self->{fsu}->catfile($lib_dir, 'libraries.txt');
    unless (-s $lib_out_file) {
        #*** also need to properly check that this was created ok
        open(my $lofh, '>', $lib_out_file) || $self->throw("Could not write to $lib_out_file");
        my $lib_file = $self->{fsu}->catfile($self->{outdir}, 'variants', 'variants.libraries.txt');
        open(my $lfh, $lib_file) || $self->throw("Could not open $lib_file");
        while (<$lfh>) {
            if (/^#/) {
                chomp;
                my @a = split;
                shift(@a);
                print $lofh '#LIB ', join('_', @a), "\n";
            }
            else {
                print $lofh $_;
            }
        }
        close($lfh);
        close($lofh);
    }
    
    # create dindel jobs - divide jobs over blocks/subdirectories with at most
    # files_per_block files
    my $blocks_dir = $self->{fsu}->catfile($self->{outdir}, 'blocks');
    unless (-d $blocks_dir) {
        mkdir($blocks_dir) || $self->throw("Could not create directory '$blocks_dir'");
    }
    my @split_var_files = @{$self->split_variants_provides($lane_path)};
    my $splits = 0;
    my $blocks = 0;
    my $block_dir;
    my $jobs = 0;
    foreach my $var_file (@split_var_files) {
        $var_file = $self->{fsu}->catfile($lane_path, $var_file);
        
        $splits++;
        if ($splits % $self->{files_per_block} == 1) {
            $blocks++;
            $block_dir = $self->{fsu}->catfile($blocks_dir, 'block.'.$blocks);
            unless (-d $block_dir) {
                mkdir($block_dir) || $self->throw("Could not create directory '$block_dir'");
            }
        }
        
        my $done = 0;
        foreach my $type ('calls', 'glf', 'log') {
            my $out_file = $self->{fsu}->catfile($block_dir, "dindel.$splits.$type.txt.gz");
            if (-s $out_file) {
                $done++;
            }
        }
        next if $done == 3;
        
        my $job_base_name = "varfile.$splits";
        my $job_name = $self->{fsu}->catfile($block_dir, $job_base_name);
        my $lock_file = $job_name.'.jids';
        
        my $is_running = LSF::is_job_running($lock_file);
        if ($is_running & $LSF::Error) {
            warn "$job_name failed!\n";
            #unlink($lock_file);
            
            # for some reason rare jobs (1-3 per strain) are failing because
            # the gzip complains it could not find the file to compress, yet
            # the files are actually created, compressed!...
            
            next;
        }
        elsif ($is_running & $LSF::Running) {
            $jobs++;
            next;
        }
        elsif ($is_running & $LSF::Done) {
            #*** $LSF::Done means that the .o file said we were successful, but
            #    that could just mean the gzip worked... don't know how to
            #    check this properly yet
            foreach my $type ('calls', 'glf', 'log') {
                my $out_file = $self->{fsu}->catfile($block_dir, "dindel.$splits.$type.txt.gz");
                my $running_file = $self->{fsu}->catfile($block_dir, "running.dindel.$splits.$type.txt.gz");
                move($running_file, $out_file) || $self->throw("failed to move $running_file to $out_file");
            }
            next;
        }
        else {
            $jobs++;
            last if $jobs > $self->{simultaneous_jobs};
            
            $self->archive_bsub_files($block_dir, $job_base_name);
            
            # jobs can fail because the .gz already exists, causing gzip to exit
            # with code 1, so always delete the .gz files first
            foreach my $type ('log', 'glf', 'calls') {
                unlink("$block_dir/running.dindel.$splits.$type.txt.gz");
            }
            
            LSF::run($lock_file, $block_dir, $job_base_name, {bsub_opts => '-q long'},
                     qq{$self->{dindel_bin} --analysis indels --bamFiles $self->{bamfiles_fofn} --varFile $var_file --ref $self->{ref} --outputFile $block_dir/running.dindel.$splits --mapUnmapped --libFile $lib_out_file $self->{dindelPars} $self->{addDindelOpt} > $block_dir/running.dindel.$splits.log.txt; gzip $block_dir/running.dindel.$splits.*.txt});
        }
        
        last if $jobs >= $self->{simultaneous_jobs};
    }
    
    return $self->{No};
}

=head2 merge_requires

 Title   : merge_requires
 Usage   : my $required_files = $obj->merge_requires('/path/to/lane');
 Function: Find out what files the merge action needs before
           it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub merge_requires {
    my ($self, $lane_path) = @_;
    return $self->call_provides($lane_path);
}

=head2 merge_provides

 Title   : merge_provides
 Usage   : my $provided_files = $obj->merge_provides('/path/to/lane');
 Function: Find out what files the merge action generates on
           success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub merge_provides {
    my $self = shift;
    return [$self->{fsu}->catfile($self->{outdir_basename}, 'calls.txt'),
            $self->{fsu}->catfile($self->{outdir_basename}, 'indels.txt')];
}

=head2 merge

 Title   : merge
 Usage   : $obj->merge('/path/to/lane', 'lock_filename');
 Function: Merge results and output final result file.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub merge {
    my ($self, $lane_path, $action_lock) = @_;
    
    # filter output from .glf files into single file, and call genotypes from
    # these
    my $blocks_dir = $self->{fsu}->catfile($self->{outdir}, 'blocks');
    
    my $job_name = $self->{fsu}->catfile($self->{logs}, 'merge');
    $self->archive_bsub_files($self->{logs}, 'merge');
    
    my $orig_bsub_opts = $self->{bsub_opts};
    $self->{bsub_opts} = ' -M5500000 -R \'select[mem>5500] rusage[mem=5500]\'';
    
    LSF::run($action_lock, $lane_path, $job_name, $self,
             qq{python $self->{dindel_scripts}/MergeOutput.py -w merge -d $blocks_dir -e glf.txt.gz -t newdindel -o results.merged.glf.txt.gz; python $self->{dindel_scripts}/ParseDindelGLF.py -w callDiploidGLF -g $blocks_dir/results.merged.glf.txt.sorted.gz -o $self->{outdir}/running.calls.txt});
    
    $self->{bsub_opts} = $orig_bsub_opts;
    
    return $self->{No};
}

sub is_finished {
    my ($self, $lane_path, $action) = @_;
    
    my $action_name = $action->{name};
    
    if ($action_name eq 'get_variants') {
        my $provided = $self->get_variants_provides($lane_path);
        foreach my $file (@{$provided}) {
            my $finished = $self->{fsu}->catfile($lane_path, $file);
            next if -s $finished;
            
            my $basename = basename($finished);
            my $running_basename = 'running.'.$basename;
            my $running = $finished;
            $running =~ s/$basename$/$running_basename/;
            
            if (-s $running) {
                # check the bsub output file to see if we completed successfully
                # (could do seperate tests for writing all indels for each
                #  ref sequence, and for writing the lib file, but just use
                #  the same test for the 'done!' claim for now
                my $bsub_o = $self->{fsu}->catfile($self->{logs}, 'get_variants.o');
                if (-s $bsub_o) {
                    open(my $bofh, $bsub_o) || $self->throw("Could not open $bsub_o");
                    my $done = 0;
                    while (<$bofh>) {
                        if (/^done!$/) {
                            $done++;
                        }
                        if (/^Successfully completed./) {
                            $done++;
                        }
                    }
                    close($bofh);
                    
                    if ($done == 2) {
                        move($running, $finished);
                    }
                }
            }
        }
    }
    if ($action_name eq 'merge') {
        my ($finished_calls, $finished_indels) = @{$self->merge_provides($lane_path)};
        $finished_calls = $self->{fsu}->catfile($lane_path, $finished_calls);
        $finished_indels = $self->{fsu}->catfile($lane_path, $finished_indels);
        
        unless (-s $finished_calls) {
            my $basename = basename($finished_calls);
            my $running_basename = 'running.'.$basename;
            my $running = $finished_calls;
            $running =~ s/$basename$/$running_basename/;
            
            if (-s $running) {
                # check the bsub output file to see if we completed successfully
                my $bsub_o = $self->{fsu}->catfile($self->{logs}, 'merge.o');
                if (-s $bsub_o) {
                    open(my $bofh, $bsub_o) || $self->throw("Could not open $bsub_o");
                    my $done = 0;
                    while (<$bofh>) {
                        if (/^Successfully completed./) {
                            $done++;
                        }
                    }
                    close($bofh);
                    
                    if ($done == 1) {
                        move($running, $finished_calls);
                    }
                }
            }
        }
        
        if (-s $finished_calls && ! -e $finished_indels) {
            # the third column contains the phred-scaled posterior probability,
            # I generally use a minimum of 20 for this, corresponding to 99 %
            # confidence.
            #
            # Columns 10 and 11 correspond to number of reads covering the indel
            # site on the forward strand, and number of reads covering the indel
            # on the reverse strand. (Even requiring a minimum of 5 on both
            # strands gives 700K calls)
            #
            # The files also contains SNP calls. THe $6 ~ /[-+]/ requires that
            # the most likely non-reference allele is an indel.
            #
            # The fourth column gives the phred-scaled confidence in the
            # genotype. The difference with the third column is that you may be
            # very confident about the statement that a non-reference allele is
            # present, but not about the precise genotype.
            #
            # The 7th column ('genotype') contains the most likely genotype (the
            # one with the highest posterior), where */* means ref/ref.
            my $problem = system("less -S $finished_calls | awk '\$3>20 && \$10>5 && \$11>5 && \$6 ~ /[-+]/' > $finished_indels.running");
            if ($problem) {
                $self->throw("Failed to filter calls to $finished_indels: $!");
            }
            else {
                move("$finished_indels.running", $finished_indels);
            }
        }
    }
    
    return $self->SUPER::is_finished($lane_path, $action);
}

# cleanup _get_variants.jids _job_status _log in $lane_path


# we override running_status to allow calls() to be called even while some jobs
# are still running, because it is doing 200 jobs at a time and we don't want
# to wait for all 200 to finish before the next 200 are scheduled.
sub running_status {
    my ($self, $jids_file) = @_;
    if ($jids_file =~ /call/) {
        return $LSF::No;
    }
    
    return $self->SUPER::running_status($jids_file);
}

1;
