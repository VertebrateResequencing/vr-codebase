=head1 NAME

VertRes::Pipelines::Dindel - pipeline for running Dindel on bams

=head1 SYNOPSIS

# Make multiple files for each of your populations/strains/other_sample_grouping
# containing absolute paths to all the bams for all the samples in that group.
# Dindel is only known to work properly with SLX (Illumina) bams, so it's best
# to limit to those as well. Eg. for a 1000 genomes directory structure:
find /path/to/REL -maxdepth 3 -mindepth 3 -type d -name SLX | grep low_coverage > samples.fod
cat samples.fod | perl -ne 'chomp; $path = $_; ($pop, $na) = $path =~ qr{/(\w\w\w)_low_coverage/([^/]+)/SLX}; $na || next; open($fh, ">>$pop.bams.fofn"); print $fh "$path/release.bam\n"; close($fh);'

# make a conf file with root pointing to anywhere you'd like the dindel results
# to appear, and that specifies the group => fofn mapping.
# Optional settings also go here.
# Your dindel.conf file may look like:
root    => '/path/to/REL/dindel',
module  => 'VertRes::Pipelines::Dindel',
prefix  => '_',
sample_groups => {
    panel_A => 'A.fofn',
    panel_B => 'B.fofn',
    panel_C => 'C.fofn',
},
data => {
    ref => '/path/to/ref.fa',
    simultaneous_jobs => 100,
    type => 'pooled',
    dindel_args => '--maxRead 5000',
}

# make a pipeline file:
echo "__DINDEL__ dindel.conf" > dindel.pipeline

# run the pipeline:
run-pipeline -c dindel.pipeline -v

# (and make sure it keeps running by adding that last to a regular cron job)

=head1 DESCRIPTION

A module for running Dindel on groups of bams. A 'sample_group' may specify
a fofn with only 1 bam or many, and in either case the 'type' option may be set
to 'diploid' if you want to treat all bams as being for the same sample. If your
bams are all for different samples, type should be set to 'pooled'.

NB: Dindel in pooled mode is hard-coded to assume that different bams are for
different samples (individuals), so if you have bams split by chr you should run
multiple independant pipelines for each chr and manually merge the final VCFs
together once all the pipelines complete.

If you want to analyse across multiple populations it is recommended that your
sample_groups consist of those multiple populations. However, if you want to
analyse populations seperately but still increase sensitivity and gain
'completeness' by calling indels in each population seen in other populations,
make use of the group_groups option:
sample_groups => {
    population_A => 'A.fofn',
    population_B => 'B.fofn',
    population_C => 'C.fofn',
    population_D => 'D.fofn'
},
group_groups =>  {
    panel_AB => ['population_A', 'population_B'],
    panel_CD => ['population_C', 'population_D']
}
NB: group_groups is not yet fully implemented!

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
use VertRes::LSF;
use Utils;

use base qw(VertRes::Pipeline);

our $actions = [{ name     => 'extract_indels',
                  action   => \&extract_indels,
                  requires => \&extract_indels_requires,
                  provides => \&extract_indels_provides },
                { name     => 'concat_libvar',
                  action   => \&concat_libvar,
                  requires => \&extract_indels_provides,
                  provides => \&concat_libvar_provides },
                { name     => 'select_candidates',
                  action   => \&select_candidates,
                  requires => \&concat_libvar_provides,
                  provides => \&select_candidates_provides },
                { name     => 'filter_candidates',
                  action   => \&filter_candidates,
                  requires => \&concat_libvar_provides,
                  provides => \&filter_candidates_provides },
                { name     => 'make_windows',
                  action   => \&make_windows,
                  requires => \&select_candidates_provides,
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
                make_windows_bsub_opts => " -M1400000 -R 'select[mem>1400] rusage[mem=1400]'",
                realign_windows_bsub_opts => " -q long -M2800000 -R 'select[mem>2800] rusage[mem=2800]'",
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
           type => 'diploid'|'pooled' (default diploid if 1 bam, pooled if many
                                       bams)
           retry_from_others => \@vcfs (in group_groups mode, when calls are
                                        complete for this sample group and the
                                        other sample groups in the group_group,
                                        the pipeline is repeated using
                                        candidates from the other sample groups
                                        that were not attempted for this one)

           simultaneous_jobs => int (default 200; the number of jobs to
                                     do at once - limited to avoid IO
                                     problems)
           dindel_args => '--maxRead 5000' (specify special non-default
                                            settings to the dindel window
                                            realignment call)
           realign_windows_bsub_opts => ' -M7900000 -R \'select[mem>7900] rusage[mem=7900]\'',
           make_windows_bsub_opts => ' -M7900000 -R \'select[mem>7900] rusage[mem=7900]\'',
                                      (specify bsub options to make_windows)
           min_count => int (default 2 in pooled mode (multiple bams supplied),
                             default 0==off in other modes; the number of times
                             a candidate indel needs to have been seen before
                             being considered)
           other optional args as per VertRes::Pipeline

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(%options, actions => $actions, @args);
    
    $self->{lane_path} || $self->throw("lane_path (misnomer; actually dindel group output) directory not supplied, can't continue");
    $self->{ref} || $self->throw("ref not supplied, can't continue");
    -d $self->{dindel_scripts} || $self->throw("bad dindel_scripts $self->{dindel_scripts}");
    $self->{bam_fofn} || $self->throw("bam_fofn not supplied, can't continue");
    
    $self->{io} = VertRes::IO->new;
    $self->{fsu} = VertRes::Utils::FileSystem->new;
    
    my @bams = $self->{io}->parse_fofn($self->{bam_fofn}, "/");
    $self->{bam_files} = \@bams;
    my @bais;
    my %seen_bams;
    foreach my $bam (@bams) {
        if (exists $seen_bams{$bam}) {
            $self->throw("bam_fofn $self->{bam_fofn} contains the bam '$bam' more than once, which dindel will fail on");
        }
        $seen_bams{$bam} = 1;
        push(@bais, $bam.'.bai');
    }
    $self->{bai_files} = \@bais;

    $self->{type} ||= @bams > 1 ? 'pooled' : 'diploid';
    $self->{min_count} ||= @bams > 1 ? 2 : 0;
    if ($self->{retry_from_others}) {
        $self->{min_count} = 0;
    }
    
    return $self;
}

sub _get_bam_to_sample {
    my $self = shift;
    
    return if defined $self->{bam_to_sample};
    
    my %bam_to_sample;
    my $bts_file = File::Spec->catfile($self->{lane_path}, 'bam_to_sample.txt');
    if (-s $bts_file) {
        open(my $fh, $bts_file) || $self->throw("Could not open $bts_file");
        while (<$fh>) {
            my ($bam, $sample) = split;
            $bam_to_sample{$bam} = $sample;
        }
    }
    else {
        open(my $ofh, '>', $bts_file) || $self->throw("Could not write to $bts_file");
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
            
            print $ofh "$bam\t$samples[0]\n";
            $bam_to_sample{$bam} = $samples[0];
        }
        close($ofh);
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
    
    warn "extract_indels_requires called!\n";
    my @reqs = ($self->{ref}, @{$self->{bam_files}}, @{$self->{bai_files}});
    
    if ($self->{retry_from_others}) {
        # we need the final calls.vcf from the other sample groups and ourselves
        push(@reqs, @{$self->{retry_from_others}});
        push(@reqs, File::Spec->catfile($self->{lane_path}, 'calls.vcf'));
    }
    
    return \@reqs;
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
    
    my @lib_var_files;
    if ($self->{retry_from_others}) {
        # we provide a single variants file and we already made the single
        # library file in the previous run of this pipeline
        push(@lib_var_files, File::Spec->catfile($lane_path, 'extra_variants.txt'),
                             File::Spec->catfile($lane_path, 'libraries.txt'));
    }
    else {
        # we provide a library and variant file per bam, listed in the
        # .extracts_done file
        @lib_var_files = (File::Spec->catfile($lane_path, '.extracts_done'));
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
        
        my $is_running = VertRes::LSF::is_job_running($lock_file);
        if ($is_running & $VertRes::LSF::Error) {
            warn "$job_name failed!\n";
            unlink($lock_file);
            next;
        }
        elsif ($is_running & $VertRes::LSF::Running) {
            next;
        }
        elsif ($is_running & $VertRes::LSF::Done) {
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
            
            VertRes::LSF::run($lock_file, $out_dir, $job_base_name, $self,
                     qq{$self->{dindel_bin} --analysis getCIGARindels --bamFile $bam --ref $self->{ref} --outputFile $running_base});
        }
    }
    
    return $self->{No};
}

sub concat_libvar_provides {
    my ($self, $lane_dir) = @_;
    my $varfile = $self->{retry_from_others} ? 'extra_variants.txt' : 'variants.txt';
    return [$self->{fsu}->catfile($lane_dir, $varfile), $self->{fsu}->catfile($lane_dir, 'libraries.txt')];
}

=head2 concat_libvar

 Title   : concat_libvar
 Usage   : $obj->concat_libvar('/path/to/lane', 'lock_filename');
 Function: Concatenates the individual sample library and variant files made
           by extract_indels() to give one file of each type per this sample
           group. Also sorts the variants file by chr,pos.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub concat_libvar {
    my ($self, $lane_path, $action_lock) = @_;
    
    my $var_out = $self->{fsu}->catfile($lane_path, 'variants.txt.running');
    my $lib_out = $self->{fsu}->catfile($lane_path, 'libraries.txt.running');
    open(my $vofh, '>', $var_out) || $self->throw("Could not write to $var_out");
    open(my $lofh, '>', $lib_out) || $self->throw("Could not write to $lib_out");
    
    # concat
    my ($v_exp, $l_exp) = (0, 0);
    my $ex_done_file = File::Spec->catfile($lane_path, '.extracts_done');
    open(my $edffh, $ex_done_file) || $self->throw("Could not open $ex_done_file");
    while (<$edffh>) {
        chomp;
        my $file = $_;
        $file || next;
        
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
    close($edffh);
    
    # sort variants file
    system("sort -k1,1n -k2,2n $var_out > $var_out.sorted; mv $var_out.sorted $var_out");
    
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

sub select_candidates_provides {
    my ($self, $lane_dir) = @_;
    my $varfile = $self->{retry_from_others} ? 'extra_variants.txt' : 'selected_variants.txt';
    return [$self->{fsu}->catfile($lane_dir, $varfile), $self->{fsu}->catfile($lane_dir, 'libraries.txt')];
}

=head2 select_candidates

 Title   : select_candidates
 Usage   : $obj->select_candidates('/path/to/lane', 'lock_filename');
 Function: Select candidate variants to work on in subsequent steps.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub select_candidates {
    my ($self, $lane_path, $action_lock) = @_;
    
    my $var_file = $self->{fsu}->catfile($lane_path, 'variants.txt');
    my $sel_file = $self->{fsu}->catfile($lane_path, 'selected_variants.txt');
    my $min_count = $self->{min_count};
    
    if ($min_count == 0) {
        copy($var_file, $sel_file) || $self->throw("Failed to cp $var_file to $sel_file");
        unless (-s $var_file == -s $sel_file) {
            unlink($sel_file);
            $self->throw("cp of $var_file to $sel_file was bad; $sel_file unlinked");
        }
        return $self->{Yes};
    }
    
    my $job_basename = 'select_candidates';
    my $job_name = $self->{fsu}->catfile($lane_path, $job_basename);
    $self->archive_bsub_files($lane_path, $job_basename);
    
    VertRes::LSF::run($action_lock, $lane_path, $job_basename, $self,
             qq{python $self->{dindel_scripts}/selectCandidates.py --minCount $min_count -i $var_file -o $sel_file.running});
    
    return $self->{No};
}

sub filter_candidates_provides {
    my ($self, $lane_dir) = @_;
    my $varfile = '.filter_candidates_done';
    return [$self->{fsu}->catfile($lane_dir, $varfile)];
}

=head2 filter_candidates

 Title   : filter_candidates
 Usage   : $obj->filter_candidates('/path/to/lane', 'lock_filename');
 Function: Filter candidate variants from the previous step, exclude variants not present in the 'candidate_filter' VCF.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub filter_candidates {
    my ($self, $lane_path, $action_lock) = @_;

    my $done_file = File::Spec->catfile($lane_path, '.filter_candidates_done');
    if ( !exists($$self{filter_candidates}) ) 
    { 
        system("touch $done_file"); 
        return $$self{Yes};
    }

    my $candidates = $self->{fsu}->catfile($lane_path, 'selected_variants.txt');
    my $filter = $$self{filter_candidates};
    my $win    = exists($$self{filter_window}) ? $$self{filter_window} : 0;

    my $job_basename = 'filter_candidates';
    my $job_name = $self->{fsu}->catfile($lane_path, $job_basename);
    $self->archive_bsub_files($lane_path, $job_basename);

    VertRes::LSF::run($action_lock, $lane_path, $job_basename, $self,
            qq{perl -MVertRes::Pipelines::Dindel -e '\\''VertRes::Pipelines::Dindel->filter_candidates_run(q[$win],q[$candidates],q[$filter],q[$done_file])'\\''});
    
    return $self->{No};
}

sub filter_candidates_run
{
    my ($self,$win,$candidates,$filter,$done_file) = @_;

    use Utils;
    Utils::CMD("cat $candidates | sed 's,\\s\\s*,\\t,g' | bgzip -c > $candidates.unfiltered.tab.gz",{verbose=>1});
    Utils::CMD("tabix -f -s 1 -b 2 -e 2 $candidates.unfiltered.tab.gz",{verbose=>1});
    Utils::CMD("vcf-isec -o -w $win -n =2 -t 1:2:$candidates.unfiltered.tab.gz -t 1:2:$filter | sed 's,\\t, ,g' > $candidates.out",{verbose=>1});
    Utils::CMD("vcf-isec -o -w $win -c -t 1:2:$candidates.unfiltered.tab.gz -t 1:2:$filter > $candidates.dindel-only",{verbose=>1});
    Utils::CMD("vcf-isec -o -w $win -c -t 1:2:$filter -t 1:2:$candidates.unfiltered.tab.gz > $candidates.filter-only",{verbose=>1});
    rename($candidates,"$candidates.unfiltered") or $self->throw("rename $candidates $candidates.unfiltered");
    rename("$candidates.out",$candidates) or $self->throw("rename $candidates.out $candidates");
    Utils::CMD("touch $done_file",{verbose=>1});
}

sub make_windows_provides {
    my ($self, $lane_path) = @_;
    
    # we don't know how many windows files there will be, but we can work out
    # how many lines should be found across all the window files, and once we
    # pass that check we'll make a .made_windows file listing all the window
    # files
    
    my $complete = $self->{retry_from_others} ? '.made_windows_retry' : '.made_windows';
    return ['.made_windows'];
}

sub _check_all_windows_files {
    my ($self, $lane_path) = @_;
    
    my $done_base = $self->{retry_from_others} ? '.made_windows_retry' : '.made_windows';
    my $done_file = $self->{fsu}->catfile($lane_path, $done_base);
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
    
    my $var_base = $self->{retry_from_others} ? 'extra_variants.txt' : 'selected_variants.txt';;
    my $var_file = $self->{fsu}->catfile($lane_path, $var_base);
    my $window_base = $self->{retry_from_others} ? 'windows_retry' : 'windows';
    my $window_dir = $self->{fsu}->catfile($lane_path, $window_base);
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
    
    my $done_base = $self->{retry_from_others} ? '.made_windows_retry' : '.made_windows';
    my $done_file = $self->{fsu}->catfile($lane_path, $done_base);
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
        $self->throw("$done_base not present, this method should not have been called");
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
    
    my $var_base = $self->{retry_from_others} ? 'extra_variants.txt' : 'selected_variants.txt';;
    my $var_file = $self->{fsu}->catfile($lane_path, $var_base);
    my $window_base = $self->{retry_from_others} ? 'windows_retry' : 'windows';
    my $window_dir = $self->{fsu}->catfile($lane_path, $window_base);
    unless (-d $window_dir) {
        mkdir($window_dir) || $self->throw("Could not create directory '$window_dir'");
    }
    
    return $self->{Yes} if $self->_check_all_windows_files($lane_path);
    
    my $job_basename = 'dindel_make_windows';
    my $job_name = $self->{fsu}->catfile($window_dir, $job_basename);
    $self->archive_bsub_files($window_dir, $job_basename);

    my $orig_bsub_opts = $self->{bsub_opts};
    $self->{bsub_opts} = $self->{make_windows_bsub_opts};
    
    VertRes::LSF::run($action_lock, $window_dir, $job_basename, $self,
             qq{python $self->{dindel_scripts}/makeWindows.py --inputVarFile $var_file --windowFilePrefix $window_dir/window --numWindowsPerFile 1000});
    
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
    my $window_base = $self->{retry_from_others} ? 'windows_retry' : 'windows';
    my $window_dir = $self->{fsu}->catfile($lane_path, $window_base);
    my @window_files = sort { my ($n1) = $a =~ /(\d+)\.txt/; my ($n2) = $b =~ /(\d+)\.txt/; $n1 <=> $n2 } $self->_get_windows_files($lane_path);
    
    my $orig_bsub_opts = $self->{bsub_opts};
    $self->{bsub_opts} = $self->{realign_windows_bsub_opts};
    
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
        
        my $is_running = VertRes::LSF::is_job_running($lock_file);
        if ($is_running & $VertRes::LSF::Error) {
            warn "$job_name failed!\n";
            unlink($lock_file);
            next;
        }
        elsif ($is_running & $VertRes::LSF::Running) {
            $jobs++;
            next;
        }
        elsif ($is_running & $VertRes::LSF::Done) {
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
            
            if ("$last_index" ne "index" && $last_index == $expected_index) {
                move($running_file, $glf) || $self->throw("failed to move $running_file to $glf");
            }
            else {
                if ("$last_index" eq "index") {
                    $self->warn("Made a glf file $glf, but it is just the header instead of $expected_index records; moving it to .bad");
                }
                else {
                    $self->warn("Made a glf file $glf, but it ended on window file index $last_index instead of $expected_index; moving it to .bad");
                }
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
            my $mode = $self->{type} eq 'diploid' ? '--doDiploid' : '--doPooled';
            my $files = @bam_files > 1 ? "--bamFiles $self->{bam_fofn}" : "--bamFile @bam_files";
            
            VertRes::LSF::run($lock_file, $window_dir, $job_base_name, $self,
                     qq{$self->{dindel_bin} --analysis indels $files $mode $self->{dindel_args} --ref $self->{ref} --varFile $window_file --libFile $lib_file --outputFile $running_base});
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
    
    my $job_basename = 'dindel_merge';
    my $job_name = $self->{fsu}->catfile($lane_path, $job_basename);
    $self->archive_bsub_files($lane_path, $job_basename);
    
    my $script = '';
    if ($self->{type} eq 'pooled') {
        $script = "mergeOutputPooled.py --bamFiles $self->{bam_fofn}";
    }
    else {
        $script = "mergeOutputDiploid.py";
    }
    
    my $running_out = $self->{fsu}->catfile($lane_path, 'calls.vcf.running');
    
    VertRes::LSF::run($action_lock, $lane_path, $job_basename, $self,
             qq{python $self->{dindel_scripts}/$script --inputFiles $glf_fofn --outputFile $running_out --ref $self->{ref}});
    
    return $self->{No};
}

sub is_finished {
    my ($self, $lane_path, $action) = @_;
    
    my $action_name = $action->{name};
    
    if ($action_name eq 'extract_indels') {
        my $done_file = $self->{fsu}->catfile($lane_path, '.extracts_done');
        
        unless ($self->{fsu}->file_exists($done_file)) {
            my @bam_files = @{$self->{bam_files}};
            $self->_get_bam_to_sample;
            
            my @lib_var_files;
            my $have = 0;
            foreach my $bam (@bam_files) {
                my $sample = $self->{bam_to_sample}->{$bam} || $self->throw("No sample for bam '$bam'!");
                my $out_dir = $self->{fsu}->catfile($lane_path, 'extract_indels');
                my $out_base = $self->{fsu}->catfile($out_dir, $sample.'.dindel_extract_indels');
                
                foreach my $type ('libraries', 'variants') {
                    my $out_file = $out_base.".$type.txt";
                    $have++ if -e $out_file;
                    push(@lib_var_files, $out_file);
                }
            }
            
            if ($have == @lib_var_files) {
                open(my $ofh, '>', $done_file) || $self->throw("Could not write to $done_file");
                print $ofh join("\n", @lib_var_files), "\n";
                close($ofh);
            }
        }
    }
    elsif ($action_name eq 'select_candidates') {
        my $sel_file = $self->{fsu}->catfile($lane_path, 'selected_variants.txt');
        my $running = $sel_file.'.running';
        
        if (-s $running) {
            my $lock_file = $self->{fsu}->catfile($lane_path, $self->{prefix}.'select_candidates.jids');
            
            my $is_running = VertRes::LSF::is_job_running($lock_file);
            if ($is_running & $VertRes::LSF::Error) {
                warn "$lock_file indicates failure\n";
                unlink($lock_file);
            }
            elsif ($is_running & $VertRes::LSF::Done) {
                move($running, $sel_file) || $self->throw("failed to move $running to $sel_file");
            }
        }
    }
    elsif ($action_name eq 'make_windows') {
        $self->_check_all_windows_files($lane_path);
    }
    elsif ($action_name eq 'merge') {
        my $vcf = $self->{fsu}->catfile($lane_path, 'calls.vcf');
        my $running = $vcf.'.running';
        
        if (! $self->{fsu}->file_exists($vcf) && -s $running) {
            my $lock_file = $self->{fsu}->catfile($lane_path, $self->{prefix}.'merge.jids');
            
            my $is_running = VertRes::LSF::is_job_running($lock_file);
            if ($is_running & $VertRes::LSF::Error) {
                warn "$lock_file indicates failure\n";
                unlink($lock_file);
            }
            elsif ($is_running & $VertRes::LSF::Done) {
                 #*** check for truncation? we could check the .o file which
                 #    lists all the glf files calls were made from...
                move($running, $vcf) || $self->throw("failed to move $running to $vcf");
            }
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
        return $VertRes::LSF::No;
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
