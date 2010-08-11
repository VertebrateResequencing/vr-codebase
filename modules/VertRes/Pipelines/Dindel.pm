=head1 NAME

VertRes::Pipelines::Dindel - pipeline for running Dindel on bams

=head1 SYNOPSIS

# Make multiple files for each of your populations/strains/other_sample_grouping
# containing absolute paths to all the bams for all the samples in that group.
# Even if you want to ultimately analysis multiple populations together, you
# should create seperate files for each population. Dindel is only known to work
# properly with SLX (Illumina) bams, so it's best to limit to those as well.
# Eg. for a 1000 genomes directory structure:
find /path/to/REL -maxdepth 3 -mindepth 3 -type d -name SLX | grep low_coverage > samples.fod
cat samples.fod | perl -ne 'chomp; $path = $_; ($pop, $na) = $path =~ qr{/(\w\w\w)_low_coverage/([^/]+)/SLX}; $na || next; open($fh, ">>$pop.bams.fofn"); print $fh "$path/release.bam\n"; close($fh);'

# make a conf file with root pointing to anywhere you'd like the dindel results
# to appear, and that specifies the group => fofn mapping.
# Optional settings also go here. If you wanted to analyses multiple groups
# together, make use of the group_groups option.
# Your dindel.conf file may look like:
root    => '/path/to/REL/dindel',
module  => 'VertRes::Pipelines::Dindel',
prefix  => '_',
sample_groups => {
    population_A => 'A.fofn',
    population_B => 'B.fofn',
    population_C => 'C.fofn',
    population_D => 'D.fofn'
},
group_groups =>  {
    group_AB => ['population_A', 'population_B'],
    group_CD => ['population_C', 'population_D']
},
data => {
    ref => '/path/to/ref.fa',
    simultaneous_jobs => 100,
    type => 'diploid',
    dindel_args => '--maxRead 5000',
}

# make a pipeline file:
echo "__DINDEL__ dindel.conf" > dindel.pipeline

# run the pipeline:
run-pipeline -c dindel.pipeline -v

# (and make sure it keeps running by adding that last to a regular cron job)

=head1 DESCRIPTION

A module for running Dindel on groups of bams. If a 'sample_group' specifies
a fofn with only 1 bam, dindel will be run in diploid mode. If it has multiple
bams, dindel will be run in pool mode.

NB: Dindel is hard-coded to assume that different bams are for different
samples (individuals), so if you have bams split by chr you should run multiple
independant pipelines for each chr and manually merge the final VCFs together
once all the pipelines complete.

For 1000 genomes there is a script 'dindel_pipeline_prepare.pl' that will
generate all the fofn, conf and the pipeline files for you.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Pipelines::Dindel;

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

use base qw(VertRes::Pipeline);

our $actions = [{ name     => 'extract_indels',
                  action   => \&extract_indels,
                  requires => \&extract_indels_requires, 
                  provides => \&extract_indels_provides },
                { name     => 'concat_libvar',
                  action   => \&concat_libvar,
                  requires => \&extract_indels_provides, 
                  provides => \&concat_libvar_provides },
                { name     => 'make_windows',
                  action   => \&make_windows,
                  requires => \&concat_libvar_provides, 
                  provides => \&make_windows_provides },
                { name     => 'realign_windows',
                  action   => \&realign_windows,
                  requires => \&realign_windows_requires, 
                  provides => \&realign_windows_provides },
                { name     => 'merge',
                  action   => \&merge,
                  requires => \&realign_windows_provides, 
                  provides => \&merge_provides }];
                
my $dindel_base = $ENV{DINDEL_SCRIPTS} || die "DINDEL_SCRIPTS environment variable not set\n";
our %options = (simultaneous_jobs => 100,
                dindel_scripts => $dindel_base,
                dindel_bin => 'dindel',
                dindel_args => '--maxRead 5000',
                bsub_opts => '');

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Pipelines::Dindel->new(lane => '/path/to/lane');
 Function: Create a new VertRes::Pipelines::Dindel object.
 Returns : VertRes::Pipelines::Dindel object
 Args    : lane_path => '/path/to/dindel_group_dir' (REQUIRED, set by
                         run-pipeline automatically)
           bam_fofn => '/path/to/bam.fofn' (REQUIRED, set by run-pipeline
                                            automatically)
           ref => '/path/to/ref.fa' (REQUIRED)
           type => 'diploid'|'haploid' (default diploid)

           simultaneous_jobs => int (default 200; the number of jobs to
                                     do at once - limited to avoid IO
                                     problems)
           dindel_args => '--maxRead 5000' (specify special non-default
                                            settings to the dindel window
                                            realignment call)
           other optional args as per VertRes::Pipeline

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(%options, actions => $actions, @args);
    
    $self->{lane_path} || $self->throw("lane_path (misnomer; actually dindel group output) directory not supplied, can't continue");
    $self->{ref} || $self->throw("ref not supplied, can't continue");
    #-x $self->{dindel_bin} || $self->throw("bad dindel_bin $self->{dindel_bin}");
    -d $self->{dindel_scripts} || $self->throw("bad dindel_scripts $self->{dindel_scripts}");
    $self->{bam_fofn} || $self->throw("bam_fofn not supplied, can't continue");
    
    $self->{io} = VertRes::IO->new;
    $self->{fsu} = VertRes::Utils::FileSystem->new;
    
    my @bams = $self->{io}->parse_fofn($self->{bam_fofn});
    $self->{bam_files} = \@bams;
    my @bais;
    my %bam_to_sample;
    foreach my $bam (@bams) {
        push(@bais, $bam.'.bai');
    }
    $self->{bai_files} = \@bais;
    
    $self->{type} ||= 'diploid';
    
    return $self;
}

sub _get_bam_to_sample {
    my $self = shift;
    
    return if defined $self->{bam_to_sample};
    
    my %bam_to_sample;
    foreach my $bam (@{$self->{bam_files}}) {
        my $pars = VertRes::Parser::sam->new(file => $bam);
        my @samples = $pars->samples();
        $pars->close;
        
        if (@samples > 1) {
            $self->throw("The bam '$bam' was for more than 1 sample (@samples)");
        }
        elsif (@samples < 1) {
            $self->throw("Could not detect what sample '$bam' was for");
        }
        
        $bam_to_sample{$bam} = $samples[0];
    }
    $self->{bam_to_sample} = \%bam_to_sample;
    return;
}

=head2 extract_indels_requires

 Title   : extract_indels_requires
 Usage   : my $required_files = $obj->extract_indels_requires('/path/to/lane');
 Function: Find out what files the extract_indels action needs before
           it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub extract_indels_requires {
    my $self = shift;
    return [$self->{ref}, @{$self->{bam_files}}, @{$self->{bai_files}}];
}

=head2 extract_indels_provides

 Title   : extract_indels_provides
 Usage   : my $provided_files = $obj->extract_indels_provides('/path/to/lane');
 Function: Find out what files the extract_indels action generates on
           success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub extract_indels_provides {
    my ($self, $lane_path) = @_;
    
    my $out_dir = $self->{fsu}->catfile($lane_path, 'extract_indels');
    my @bam_files = @{$self->{bam_files}};
    $self->_get_bam_to_sample;
    
    my @lib_var_files;
    foreach my $bam (@bam_files) {
        my $sample = $self->{bam_to_sample}->{$bam} || $self->throw("No sample for bam '$bam'!");
        my $out_base = $self->{fsu}->catfile($out_dir, $sample.'.dindel_extract_indels');
        
        foreach my $type ('libraries', 'variants') {
            my $out_file = $out_base.".$type.txt";
            push(@lib_var_files, $out_file);
        }
    }
    
    return \@lib_var_files;
}

=head2 extract_indels

 Title   : extract_indels
 Usage   : $obj->extract_indels('/path/to/lane', 'lock_filename');
 Function: Extracts candidate indels and infer insert-size distribution on each
           bam seperately, then merge.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub extract_indels {
    my ($self, $lane_path, $action_lock) = @_;
    
    my $out_dir = $self->{fsu}->catfile($lane_path, 'extract_indels');
    unless (-d $out_dir) {
        mkdir($out_dir) || $self->throw("Could not create directory '$out_dir'");
    }
    
    my @bam_files = @{$self->{bam_files}};
    $self->_get_bam_to_sample;
    
    foreach my $bam (@bam_files) {
        my $sample = $self->{bam_to_sample}->{$bam} || $self->throw("No sample for bam '$bam'!");
        my $out_base = $self->{fsu}->catfile($out_dir, $sample.'.dindel_extract_indels');
        my $running_base = $out_base.'.running';
        
        my $done = 0;
        foreach my $type ('libraries', 'variants') {
            my $out_file = $out_base.".$type.txt";
            if ($self->{fsu}->file_exists($out_file)) {
                $done++;
            }
        }
        next if $done == 2;
        
        my $job_base_name = $sample.'.dindel_extract_indels';
        my $job_name = $self->{fsu}->catfile($out_dir, $job_base_name);
        my $lock_file = $job_name.'.jids';
        
        my $is_running = LSF::is_job_running($lock_file);
        if ($is_running & $LSF::Error) {
            warn "$job_name failed!\n";
            unlink($lock_file);
            next;
        }
        elsif ($is_running & $LSF::Running) {
            next;
        }
        elsif ($is_running & $LSF::Done) {
            #*** $LSF::Done means that the .o file said we were successful
            foreach my $type ('libraries', 'variants') {
                my $out_file = $out_base.".$type.txt";
                my $running_file = $running_base.".$type.txt";
                move($running_file, $out_file) || $self->throw("failed to move $running_file to $out_file");
            }
            next;
        }
        else {
            $self->archive_bsub_files($out_dir, $job_base_name);
            
            LSF::run($lock_file, $out_dir, $job_base_name, $self,
                     qq{$self->{dindel_bin} --analysis getCIGARindels --bamFile $bam --ref $self->{ref} --outputFile $running_base});
        }
    }
    
    return $self->{No};
}

sub concat_libvar_provides {
    my ($self, $lane_dir) = @_;
    return [$self->{fsu}->catfile($lane_dir, 'variants.txt'), $self->{fsu}->catfile($lane_dir, 'libraries.txt')];
}

=head2 concat_libvar

 Title   : concat_libvar
 Usage   : $obj->concat_libvar('/path/to/lane', 'lock_filename');
 Function: Concatenates the individual sample library and variant files made
           by extract_indels() to give one file of each type per this sample
           group.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub concat_libvar {
    my ($self, $lane_path, $action_lock) = @_;
    
    my $var_out = $self->{fsu}->catfile($lane_path, 'variants.txt.running');
    my $lib_out = $self->{fsu}->catfile($lane_path, 'libraries.txt.running');
    open(my $vofh, '>', $var_out) || $self->throw("Could not write to $var_out");
    open(my $lofh, '>', $lib_out) || $self->throw("Could not write to $lib_out");
    
    my ($v_exp, $l_exp) = (0, 0);
    foreach my $file (@{$self->extract_indels_provides($lane_path)}) {
        my ($ofh, $counter);
        if ($file =~ /variants\.txt$/) {
            $ofh = $vofh;
            $counter = \$v_exp;
        }
        else {
            $ofh = $lofh;
            $counter = \$l_exp;
        }
        
        open(my $fh, $file) || $self->throw("Could not open $file");
        while (<$fh>) {
            $$counter++;
            print $ofh $_;
        }
        close($fh);
    }
    close($vofh);
    close($lofh);
    
    # check the cat files are complete
    my $v_lines = VertRes::IO->new(file => $var_out)->num_lines;
    my $l_lines = VertRes::IO->new(file => $lib_out)->num_lines;
    if ($v_lines != $v_exp) {
        $self->throw("$var_out ended up with $v_lines instead of $v_exp lines");
    }
    if ($l_lines != $l_exp) {
        $self->throw("$lib_out ended up with $l_lines instead of $l_exp lines");
    }
    
    move($var_out, $self->{fsu}->catfile($lane_path, 'variants.txt'));
    move($lib_out, $self->{fsu}->catfile($lane_path, 'libraries.txt'));
    
    return $self->{Yes};
}

sub make_windows_provides {
    my ($self, $lane_path) = @_;
    
    # we don't know how many windows files there will be, but we can work out
    # how many lines should be found across all the window files, and once we
    # pass that check we'll make a .made_windows file listing all the window
    # files
    
    return ['.made_windows'];
}

sub _check_all_windows_files {
    my ($self, $lane_path) = @_;
    
    my $done_file = $self->{fsu}->catfile($lane_path, '.made_windows');
    if ($self->{fsu}->file_exists($done_file)) {
        return 1;
    }
    
    # the number of lines in variants.txt should match the number of column 4+
    # entries in all window files...
    # actually, when there are multiple samples, they may share indels, and so
    # this won't be true.
    # Instead we'll compare the number of variants in the variants file
    # (> lines) to the STDERR output of makeWindows.py, which says how many
    # variants it read. If that matches, then we'll trust the number of
    # candidates it reports as the expected, and make sure the window files
    # match that.
    
    my $var_file = $self->{fsu}->catfile($lane_path, 'variants.txt');
    my $window_dir = $self->{fsu}->catfile($lane_path, 'windows');
    my $e_file = $self->{fsu}->catfile($window_dir, 'dindel_make_windows.e');
    
    if (-d $window_dir && -s $e_file) {
        my ($read_variants, $num_candidates) = (0, 0);
        open(my $efh, $e_file) || $self->throw("Could not open error file '$e_file'");
        while (<$efh>) {
            if (/Total variants read: (\d+)/) {
                $read_variants = $1;
            }
            elsif (/Number of candidates: (\d+)/) {
                $num_candidates = $1;
            }
        }
        close($efh);
        unless ($read_variants && $num_candidates && $num_candidates <= $read_variants) {
            return 0;
        }
        
        my $num_variants = 0;
        open(my $vfh, $var_file) || $self->throw("Could not open variants file '$var_file'");
        while (<$vfh>) {
            my ($after_hash) = $_ =~ /# (.+)$/;
            my @n = split(" ", $after_hash);
            $num_variants += @n;
        }
        close($vfh);
        unless ($num_variants == $read_variants) {
            $self->throw("read $read_variants variants, but there are $num_variants variants in total!");
        }
        
        my $actual = 0;
        opendir(my $windowfh, $window_dir) || $self->throw("Could not open dir $window_dir");
        my @window_files = ();
        foreach my $file (readdir($windowfh)) {
            if ($file =~ /window\./) {
                my $this_window_file = $self->{fsu}->catfile($window_dir, $file);
                push(@window_files, $this_window_file);
                open(my $wfh, $this_window_file) || $self->throw("Could not open window file '$this_window_file'");
                while (<$wfh>) {
                    my @a = split;
                    $actual += scalar(@a) - 3;
                }
            }
        }
        
        
        if ($actual >= $num_candidates) {
            open(my $dfh, '>', $done_file) || $self->throw("Could not write to $done_file");
            foreach my $window_file (@window_files) {
                print $dfh $window_file, "\n";
            }
            close($dfh);
            $self->{windows_files} = \@window_files;
            return 1;
        }
        else {
            warn "$actual entries in window files < $num_candidates\n";
        }
    }
    
    return 0;
}

sub _get_windows_files {
    my ($self, $lane_path) = @_;
    
    if (defined $self->{windows_files}) {
        return @{$self->{windows_files}};
    }
    
    my $done_file = $self->{fsu}->catfile($lane_path, '.made_windows');
    if ($self->{fsu}->file_exists($done_file)) {
        open(my $fh, $done_file) || $self->throw("Could not open $done_file");
        my @window_files;
        while (<$fh>) {
            chomp;
            push(@window_files, $_);
        }
        close($fh);
        $self->{windows_files} = \@window_files;
        return @window_files;
    }
    else {
        $self->throw(".made_windows not present, this method should not have been called");
    }
}

=head2 make_windows

 Title   : make_windows
 Usage   : $obj->make_windows('/path/to/lane', 'lock_filename');
 Function: Create realignment windows.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub make_windows {
    my ($self, $lane_path, $action_lock) = @_;
    
    my $var_file = $self->{fsu}->catfile($lane_path, 'variants.txt');
    my $window_dir = $self->{fsu}->catfile($lane_path, 'windows');
    unless (-d $window_dir) {
        mkdir($window_dir) || $self->throw("Could not create directory '$window_dir'");
    }
    
    return $self->{Yes} if $self->_check_all_windows_files($lane_path);
    
    my $job_name = $self->{fsu}->catfile($window_dir, 'dindel_make_windows');
    $self->archive_bsub_files($window_dir, 'dindel_make_windows');
    
    my $orig_bsub_opts = $self->{bsub_opts};
    $self->{bsub_opts} = ' -M7900000 -R \'select[mem>7900] rusage[mem=7900]\'';

    LSF::run($action_lock, $window_dir, $job_name, $self,
             qq{python $self->{dindel_scripts}/makeWindows.py --inputVarFile $lane_path/variants.txt --windowFilePrefix $window_dir/window --numWindowsPerFile 1000});
    
    $self->{bsub_opts} = $orig_bsub_opts;
    
    return $self->{No};
}

sub realign_windows_requires {
    my ($self, $lane_path) = @_;
    my @window_files = $self->_get_windows_files($lane_path);
    $self->throw("Surely too few window files?!") if @window_files < 2;
    return [$self->{fsu}->catfile($lane_path, 'libraries.txt'), @window_files];
}

sub realign_windows_provides {
    my ($self, $lane_path) = @_;
    
    my @window_files = $self->_get_windows_files($lane_path);
    
    my @glfs;
    foreach my $window_file (@window_files) {
        my $glf = $window_file;
        $glf =~ s/\.txt$/.glf.txt/;
        push(@glfs, $glf);
    }
    
    @glfs > 1 || $self->throw("Surely too few glfs?!");
    
    return \@glfs;
}

=head2 realign_windows

 Title   : realign_windows
 Usage   : $obj->realign_windows('/path/to/lane', 'lock_filename');
 Function: Realign the reads in each window (the computationally intensive bit).
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub realign_windows {
    my ($self, $lane_path, $action_lock) = @_;
    
    my $lib_file = $self->{fsu}->catfile($lane_path, 'libraries.txt');
    my $window_dir = $self->{fsu}->catfile($lane_path, 'windows');
    my @window_files = $self->_get_windows_files($lane_path);
    
    my $orig_bsub_opts = $self->{bsub_opts};
    $self->{bsub_opts} = ' -q long -M7900000 -R \'select[mem>7900] rusage[mem=7900]\'';
    
    my $jobs = 0;
    foreach my $window_file (@window_files) {
        my $glf = $window_file;
        $glf =~ s/\.txt$/.glf.txt/;
        
        next if $self->{fsu}->file_exists($glf);
        
        my $running_base = $glf;
        $running_base =~ s/\.glf\.txt/.running/;
        
        my $job_base_name = basename($running_base);
        $job_base_name =~ s/\.running//;
        my $job_name = $self->{fsu}->catfile($window_dir, $job_base_name);
        my $lock_file = $job_name.'.jids';
        
        my $is_running = LSF::is_job_running($lock_file);
        if ($is_running & $LSF::Error) {
            warn "$job_name failed!\n";
            unlink($lock_file);
            next;
        }
        elsif ($is_running & $LSF::Running) {
            $jobs++;
            next;
        }
        elsif ($is_running & $LSF::Done) {
            my $running_file = $running_base.".glf.txt";
            
            # to check the glf file for completion properly, the number of lines
            # should match the number of lines in the window file * number of
            # variants on that line (cols 4+) * number of samples, or something
            # like that. However, when a window is skipped only 1 line will
            # appear for that window line. Cheat and just check the glf file
            # ends with a line corresponding to the number of lines in the
            # window
            my $expected_index = $self->{io}->new(file => $window_file)->num_lines;
            my $last_index = 0;
            open(my $rfh, "tail -1 $running_file |") || $self->throw("Could not open tail to $running_file");
            while (<$rfh>) {
                (undef, $last_index) = split;
            }
            close($rfh);
            
            if ($last_index == $expected_index) {
                move($running_file, $glf) || $self->throw("failed to move $running_file to $glf");
            }
            else {
                $self->warn("Made a glf file $glf, but it ended on window file index $last_index instead of $expected_index; moving it to .bad");
                move($running_file, "$running_file.bad");
            }
            
            unlink($lock_file);
            next;
        }
        else {
            $jobs++;
            last if $jobs > $self->{simultaneous_jobs};
            
            $self->archive_bsub_files($window_dir, $job_base_name);
            
            my $bam_mode_args = '';
            my @bam_files = @{$self->{bam_files}};
            if (@bam_files > 1) {
                $bam_mode_args = "--bamFiles $self->{bam_fofn} --doEM";
            }
            else {
                $bam_mode_args = "--bamFile @bam_files --doDiploid";
            }
            
            LSF::run($lock_file, $window_dir, $job_base_name, $self,
                     qq{$self->{dindel_bin} --analysis indels $bam_mode_args $self->{dindel_args} --ref $self->{ref} --varFile $window_file --libFile $lib_file --outputFile $running_base});
        }
    }
    
    $self->{bsub_opts} = $orig_bsub_opts;
    
    return $self->{No};
}

sub merge_provides {
    my ($self, $lane_path) = @_;
    return [$self->{fsu}->catfile($lane_path, 'calls.vcf')];
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
    
    my $glf_fofn = $self->{fsu}->catfile($lane_path, 'glf.fofn');
    unless ($self->{fsu}->file_exists($glf_fofn)) {
        my @glf_files = @{$self->realign_windows_provides($lane_path)};
        open(my $ofh, '>', $glf_fofn) || $self->throw("Could not write to $glf_fofn");
        print $ofh join("\n", @glf_files), "\n";
        close($ofh);
    }
    
    my $job_name = $self->{fsu}->catfile($lane_path, 'dindel_merge');
    $self->archive_bsub_files($lane_path, 'dindel_merge');
    
    my $bam_mode_args = '';
    my @bam_files = @{$self->{bam_files}};
    if (@bam_files > 1) {
        $bam_mode_args = "--bamFiles $self->{bam_fofn} --type pool";
    }
    else {
        $bam_mode_args = "--bamFile @bam_files --type $self->{type}";
    }
    
    my $running_out = $self->{fsu}->catfile($lane_path, 'calls.vcf.running');
    
    LSF::run($action_lock, $lane_path, $job_name, $self,
             qq{python $self->{dindel_scripts}/mergeOutput.py $bam_mode_args --inputFiles $glf_fofn --outputFile $running_out --ref $self->{ref}});
    
    return $self->{No};
}

sub is_finished {
    my ($self, $lane_path, $action) = @_;
    
    my $action_name = $action->{name};
    
    if ($action_name eq 'make_windows') {
        $self->_check_all_windows_files($lane_path);
    }
    if ($action_name eq 'merge') {
        my $vcf = $self->{fsu}->catfile($lane_path, 'calls.vcf');
        my $running = $vcf.'.running';
        
        if (! $self->{fsu}->file_exists($vcf) && -s $running) {
            my $lock_file = $self->{fsu}->catfile($lane_path, $self->{prefix}.'merge.jids');
            
            #*** checking is_running didn't work for some reason...
            #my $is_running = LSF::is_job_running($lock_file);
            #if ($is_running & $LSF::Error) {
            #    warn "$lock_file indicates failure\n";
            #    unlink($lock_file);
            #}
            #elsif ($is_running & $LSF::Done) {
                # *** check for truncation? we could check the .o file which
                #     lists all the glf files calls were made from...
                move($running, $vcf) || $self->throw("failed to move $running to $vcf");
            #}
        }
    }
    
    return $self->SUPER::is_finished($lane_path, $action);
}

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

# we override this to use file_exists, which could be dangerous for general
# use, but is pretty much required here
sub what_files_are_missing {
    my ($self, $path, $files) = @_;

    my @missing = ();
    for my $file (@$files) {
        my $file_path = index($file, '/') == 0 ? $file : "$path/$file";
        if (! $self->{fsu}->file_exists($file_path)) {
            push @missing, $file_path;
        }
    }
    return \@missing;
}

1;
