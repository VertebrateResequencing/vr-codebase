=head1 NAME

VertRes::Pipelines::SNPs - pipeline for calling SNPs from BAM files

=head1 SYNOPSIS

# Can either provide a list of BAM files, or provide a database from
# which improved BAM files will be found and have SNPs called.
# If using a database, the calls are made independently on
# each BAM file.
# Either way, make a .conf file and .pipeline file.

# The .pipeline file depends on how the BAM files are provided.
# If providing a list of BAMs:
echo '/path/to/root/dir/of/output/files/ snps.conf' > snps.pipeline
# If using a database:
echo '__VRTrack_SNPs__ snps.conf' > snps.pipeline

# The .conf file (for running off a database) may look like:
root    => '/abs/path/to/root/data/dir',
module => 'VertRes::Pipelines::SNPs',
log     => '/abs/path/to/logfile/foo.log',
max_lanes => 10,

db  => 
{
    database => 'db_name',
    host     => 'db_host',
    port     => 42,
    user     => 'db_user',
    password => 'db_password',
},

data =>
{
    task => 'mpileup,gatk,update_db,cleanup,pseudo_genome',
   
    bsub_opts       => "-q normal -M3500000 -R 'select[type==X86_64 && mem>3500] rusage[mem=3500,thouio=1,tmp=16000]'",
    bsub_opts_long  => "-q normal -M3500000 -R 'select[type==X86_64 && mem>3500] rusage[mem=3500,thouio=1,tmp=16000]'",
    bsub_opts_qcall  => "-q normal -M3500000 -R 'select[type==X86_64 && mem>3500 && tmp>17000] rusage[mem=3500,thouio=1, tmp=17000]'",
    bsub_opts_mpileup    => "-q normal -R 'select[type==X86_64] rusage[thouio=1]'",
 
    split_size_mpileup   => 300_000_000,
    split_size_gatk      => 100_000_000,         
    split_size_qcall     => 100_000_000,         
 
    tmp_dir => '/abs/path/to/tmp/dir',
 
    mpileup_cmd => 'samtools mpileup -d 500 -C50 -m3 -F0.0002 -aug',
    qcall_cmd       => 'QCALL -ct 0.01 -snpcan',
    gatk_opts =>
    {
        all =>
        {
            verbose => 1,
            #_extras => [ '-U ALLOW_UNSET_BAM_SORT_ORDER --platform Solexa' ],
            _extras => [ '-U ALLOW_UNSET_BAM_SORT_ORDER ' ],
        },
        unified_genotyper => 
        {
            genotype_likelihoods_model => 'BOTH',         # Call both SNPs and INDELs
        },
        variant_filtration =>
        {
            filters => { HARD_TO_VALIDATE => 'MQ0 >= 4 && (MQ0 / (1.0 * DP)) > 0.1' },
            maskName => 'InDel',
            clusterWindowSize => 11,
        },
        variant_recalibrator =>
        {
            mode => 'BOTH',         # recalibrate SNPs and INDELs
            target_titv => 2.08,    # only for plotting
            ignore_filter => 'HARD_TO_VALIDATE',
            percentBadVariants => 0.05,
            maxGaussians => 4,
            training_sets => ['dbsnp,VCF,known=true,training=false,truth=false,prior=8.0,/path/to/dbsnp_123.b37.vcf.gz',
                              'hapmap,VCF,known=false,training=true,truth=true,prior=15.0,/path/to/hapmap_3.3.b37.sites.vcf.gz',
                              'omni,VCF,known=false,training=true,truth=false,prior=12.0,/path/to/1000G_omni2.5.b37.sites.vcf.gz'],
        },
        apply_variant_cuts =>
        {
            fs_filter_level => 99.0,
            mode => 'BOTH',
        },
    },
 
    # max number of running jobs per task
    max_jobs => 25,
 
    dbSNP_rod   => '/abs/path/to/dbsnp/rod/file',
    [or dbSNP_vcf   => '/abs/path/to/dbsnp/vcf/file',]
    indel_mask  => '/abs/path/to/indel/mask/file',
    fa_ref => '/abs/path/to/ref.fasta',
    fai_ref => '/abs/path/to/ref.fasta.fai',
}


# max_lanes outside the data section limits the number of lanes
# which can have SNP calling run on them simultaneously, and
# is only applicable when running off a database.
# max_jobs inside the data section is the max number of jobs allowed
# to run per caller.
# The max grand total of jobs which could be running is:
# max_lanes * (number of callers) * max_jobs.

# The possible SNP callers are mpileup, gatk, varFilter and qcall.
# Those which are specified in task => '...' will be used.
#
# Add 'cleanup' to the list of tasks to cleanup all files made
# by the pipeline, leaving you with:
# caller.vcf.gz, caller.vcf.gz.stats, caller.vcf.gz.tbi
# for each caller.  In the case of mpileup, you also get
# mpileup.unfilt.vcf.gz*, which are unfiltered calls.
# The mpileup.vcf.gz* files are the filtered calls.

# If you want to supply a list of BAM files instead of running off
# a database, then
#  1) don't provide db => {...}.  Instead, provide a file of
#     filenames in the data section with
#     file_list => '/abs/path/to/file.list',
#     The extension must be .list to stop gatk complaining.
#  2) Don't include update_db in the 'task' option. e.g. use
#     task => 'mpileup,qcall,gatk',
#     instead.
#  3) max_lanes outside the data section will be ignored
# If calling on many SNP files, it is sensible to make the
# split_size_* opions smaller, in the region of
# 1_000_000 to 5_000_000.  The numbers in the above exmaple
# work well for a single lane of exome data.

# Note that QCALL is not designed to be run on one sample,
# so there's no point in adding it to the task list when
# getting the BAMs from a database.  The pipeline will complain.

# If you don't want to use the database, but also want to call 
# separately on each of bams, can also use the config file above
# without the db section, and in the pipeline file, this instead:
echo '__SNPs__ snps.conf' > snps.pipeline

# Additional options in the data section are new_root to rebase
# output location by replacing root in the bam files with new_root.
# Also, outdir may be set to put output into a common subdir at its
# destination.

# for GATK, can mow supply the option dbSNP_vcf instead of dbSNP_rod
# if your snp sites are in vcf format rather than rod format.

=head1 DESCRIPTION

A module for calling SNPs from BAM files using any of
the callers mpileup, GATK, QCALL and varFilter.

=cut

package VertRes::Pipelines::SNPs;
use base qw(VertRes::Pipeline);

use strict;
use warnings;
use LSF;
use Utils;
use VertRes::Parser::sam;
use VertRes::Utils::FileSystem;
use VertRes::Wrapper::GATK;
use File::Copy;

our @actions =
(
    # Samtools varFilter
    {
        'name'     => 'varFilter',
        'action'   => \&varFilter,
        'requires' => \&varFilter_requires, 
        'provides' => \&varFilter_provides,
    },
    # Samtools mpileup
    {
        'name'     => 'mpileup',
        'action'   => \&mpileup,
        'requires' => \&mpileup_requires, 
        'provides' => \&mpileup_provides,
    },
    # GATK SNP calling
    {
        'name'     => 'gatk',
        'action'   => \&gatk,
        'requires' => \&gatk_requires, 
        'provides' => \&gatk_provides,
    },
    # Qcall SNP calling
    {
        'name'     => 'qcall',
        'action'   => \&qcall,
        'requires' => \&qcall_requires, 
        'provides' => \&qcall_provides,
    },
    # create pseudo-genomes for the mapped organisms
    {
        'name'     => 'pseudo_genome',
        'action'   => \&pseudo_genome,
        'requires' => \&pseudo_genome_requires, 
        'provides' => \&pseudo_genome_provides,
    },    
    # clean up all files after SNP calling
    {
        'name'     => 'cleanup',
        'action'   => \&cleanup,
        'requires' => \&cleanup_requires, 
        'provides' => \&cleanup_provides,
    },
    # when getting jobs from a database, update the finished lanes
    {
        'name'     => 'update_db',
        'action'   => \&update_db,
        'requires' => \&update_db_requires, 
        'provides' => \&update_db_provides,
    },
);

our $options = 
{
    bcftools        => 'bcftools',
    bsub_opts_varfilter  => "-q normal -R 'select[type==X86_64] rusage[thouio=1]'",
    bsub_opts_mpileup    => "-q normal -R 'select[type==X86_64] rusage[thouio=1]'",
    bsub_opts_qcall      => "-q normal -R 'select[type==X86_64] rusage[thouio=1]'",
    bsub_opts_gatk       => "-q normal -M3000000 -R 'select[type==X86_64 && mem>3000] rusage[mem=3000,thouio=1]'",
    chunks_overlap  => 250,
    fai_chr_regex   => '\d+|x|y',
    gatk_opts       => { all=>{verbose=>1} },
    dbSNP_rod       => undef,
    dbSNP_vcf       => undef,
    max_jobs        => undef,
    merge_vcf       => 'merge-vcf -d',
    mpileup_cmd     => 'samtools mpileup -C50 -aug',
    qcall_cmd       => 'QCALL -ct 0.01 -snpcan', # mouse: -pphet 0',
    sort_cmd        => 'sort',
    sam2vcf         => 'sam2vcf.pl',
    split_size_varfilter => 1_000_000,
    split_size_mpileup   => 1_000_000,
    split_size_qcall     => 1_000_000,
    split_size_gatk      => 1_000_000,
    varfilter       => 'samtools.pl varFilter -S 20 -i 20',
    pileup_rmdup    => 'pileup-rmdup',
    samtools_pileup_params => '-d 500',   # homozygous: '-r 0.0001 -d 500'
    vcf_rmdup       => 'vcf-rmdup',
    vcf_stats       => 'vcf-stats',
    filter4vcf      => qq[awk '/^#/||\\\$6>=3' | vcfutils.pl filter4vcf],
    bam_suffix      => 'intervals.bam',
};


# --------- OO stuff --------------

=head2 new

        Options    : See Pipeline.pm for general options and the code for default values.

                    bam_suffix      .. when getting jobs from a database (db =>... used outside the 
                                       data section in config) which have been improved by
                                       the bamImprovement pipeline, only consider bams whose name
                                       end in the given suffix.  This is likely to be 'calmd.bam'
                                       or 'intervals.bam', depending on how BamImprovement was run.
                                       Default: 'intervals.bam'.
                    bcftools        .. The bcftools binary for mpileup task
                    dbSNP_rod       .. The dbSNP file for GATK
                    dbSNP_vcf       .. The dbSNP vcf for GATK.  This or dbSNP_rod.
                    chunks_overlap  .. Allow small overlap between chunks
                    file_list       .. File name containing a list of bam files (e.g. the 17 mouse strains). 
                    fa_ref          .. The reference sequence in fasta format
                    fai_ref         .. The reference fai file to read the chromosomes and lengths.
                    fai_chr_regex   .. The chromosomes to be processed.
                    gatk_opts       .. Options to pass to the GATK wrapper
                    indel_mask      .. GATK precomputed indel mask for filtering SNPs [optional]
                    max_jobs        .. The maximum number of running jobs per task
                    merge_vcf       .. The merge-vcf script.
                    mpileup_cmd     .. The mpileup command.
                    pileup_rmdup    .. The script to remove duplicate positions.
                    qcall_cmd       .. The qcall command
                    sam2vcf         .. The convertor from samtools pileup format to VCF.
                    samtools_pileup_params .. The options to samtools.pl varFilter (Used by Qcall and varFilter.)
                    sort_cmd        .. Change e.g. to 'sort -T /big/space'
                    split_size_gatk       .. The size of the chunks.
                    split_size_mpileup    .. The size of the chunks.
                    split_size_qcall      .. The size of the chunks.
                    split_size_varfilter  .. The size of the chunks.
                    tmp_dir         .. Big space for QCall pileup sorting.
                    varfilter       .. The samtools varFilter command (samtools.pl varFilter).
                    vcf_rmdup       .. The script to remove duplicate positions.
                    filter4vcf      .. vcfutils.pl filter4vcf from bcftools package

=cut

sub new 
{
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(%$options,'actions'=>\@actions,@args);
    $self->write_logs(1);

    if ( !$$self{file_list} ) { $self->throw("Missing the option file_list.\n"); }
    $$self{fsu} = VertRes::Utils::FileSystem->new(reconnect_db=>1);

    return $self;
}


#---------- utilities ---------------------


# Read a list of files and dies if some of them does not exist. Returns arrayref.
sub read_files
{
    my ($self,$fname) = @_;
    open(my $fh,'<',$fname) or $self->throw("$fname: $!");
    my @files;
    while (my $line=<$fh>)
    {
        chomp($line);
        if ( !$$self{fsu}->file_exists($line) ) { $self->throw("The file does not exist: [$line]\n"); }
        push @files,$line;
    }
    close($fh);
    return \@files;
}


# Returns the common prefix of the files.
sub common_prefix
{
    my ($self,$files) = @_;
    my @paths;
    my $len = -1;
    for my $file (@$files)
    {
        my @path = split(m{/+},$file);
        if ( $len<0 || $len>scalar @path ) { $len=scalar @path; }
        push @paths, \@path;
    }
    my @common;
    for (my $i=0; $i<$len-1; $i++)
    {
        my $identical=1;
        for (my $ifile=1; $ifile<scalar @paths; $ifile++)
        {
            if ( $paths[$ifile]->[$i] ne $paths[0]->[$i] ) { $identical=0; last; }
        }
        if ( !$identical ) { last; }
        push @common, $paths[0]->[$i];
    }
    return join('/+',@common);
}


# Multiple BAM files can belong to the same sample and should be merged.
sub bam_file_groups
{
    my ($self,$files) = @_;
    my %groups;
    for my $file (@$files)
    {
        my $pars = VertRes::Parser::sam->new(file => $file);
        my @samples = $pars->samples();
        $pars->close;

        if ( !@samples )
        {
            # With the merged BAM files, this will have to be handled differently.
            $self->throw("TODO: [$file]\n");
        }

        for my $sample (@samples)
        {
            push @{$groups{$sample}}, $file
        }
    }
    return \%groups;
}


sub read_chr_lengths
{
    my ($self,$fai) = @_;

    if ( $$self{fai_chr_regex}=~/^\^/ or $$self{fai_chr_regex}=~/\$$/ )
    {
        $self->throw("The regex must not contain [^\$], this will be appended automagically.\n");
    }

    # Determine the chromosomes and their lengths
    open(my $fh,'<',$fai) or $self->throw("$fai: $!"); 
    my @chr_lengths;
    while (my $line=<$fh>)
    {
        if ( !($line=~/^($$self{fai_chr_regex})\t(\d+)/i) ) { next; }
        push @chr_lengths, { name=>$1, from=>1, to=>$2 };
    }
    close($fh);
    return \@chr_lengths;
}

# tab-separated file CHR\tFROM\tPOS
sub read_regions
{
    my ($self,$fname) = @_;

    open(my $fh,'<',$fname) or $self->throw("$fname: $!"); 
    my @regions;
    while (my $line=<$fh>)
    {
        if ( !($line=~/^(\S+)\t(\d+)\t(\d+)$/) ) { next; }
        push @regions, { name=>$1, from=>$2, to=>$3 };
    }
    close($fh);
    return \@regions;
}


# Create a list of chromosomes splitted into slightly overlapping chunks in the format chr:start-end.
sub chr_chunks
{
    my ($self,$split_size) = @_;

    # It takes too long to call pileup on big bam files - split them into chunks. Here, create
    #   the list of chunks.
    my @chunks;
    my $regions = exists($$self{regions}) ? $self->read_regions($$self{regions}) : $self->read_chr_lengths($$self{fai_ref});
    for my $region (@$regions)
    {
        my $chr     = $$region{name};
        my $pos     = $$region{from};
        my $end_pos = $$region{to};
        while ($pos<$end_pos)
        {
            my $from = $pos;
            my $to   = $from+$split_size-1;

            # GATK will fail if the region goes beyond the end of the chromosome
            if ( $to>$end_pos ) { $to=$end_pos; }

            push @chunks, "$chr:$from-$to";

            # Create some overlap between the chunks, duplicate records will be removed
            #   later. (There would be overlap SNPs anyway, at least for samtools pileup,
            #   and it would be hard to decide which of the SNPs should be kept. This way,
            #   all covering reads should be there and the SNPs should have similar/identical 
            #   depth.)
            $pos += $split_size - $$self{chunks_overlap};
            if ( $pos<1 ) { $self->throw("The split size too small [$split_size]?\n"); }
        
            # Eearly exit for debugging: do two chunks for one chromosome only for testing.
            if ( $$self{debug_chunks} && scalar @chunks>$$self{debug_chunks} ) { return \@chunks; }
        }
    }
    return \@chunks;
}


# Clean .pl .e .o and .jid files.
sub clean_files
{
    my ($self,$dir,$name,$outdir) = @_;
    if ( !-e $outdir ) { Utils::CMD("mkdir -p $outdir"); }

    if ( -e "$dir/_$name.pl" ) { rename("$dir/_$name.pl","$outdir/_$name.pl"); }
    $$self{fsu}->file_exists("$dir/_$name.pl",wipe_out=>1);

    if ( -e "$dir/_$name.o" ) { rename("$dir/_$name.o","$outdir/_$name.o"); }
    $$self{fsu}->file_exists("$dir/_$name.o",wipe_out=>1);

    if ( -e "$dir/_$name.e" ) { rename("$dir/_$name.e","$outdir/_$name.e"); }
    $$self{fsu}->file_exists("$dir/_$name.e",wipe_out=>1);

    if ( -e "$dir/_$name.jid" ) { unlink("$dir/_$name.jid"); }
    $$self{fsu}->file_exists("$dir/_$name.jid",wipe_out=>1);
}


sub merge_vcf_files
{
    my ($self,$vcfs,$name) = @_;

    my $args = join(' ',@$vcfs);
    return qq[
use strict;
use warnings;
use Utils;
Utils::CMD(qq[$$self{merge_vcf} $args | bgzip -c > $name.vcf.gz.part]);
Utils::CMD(qq[zcat $name.vcf.gz.part | $$self{vcf_stats} > $name.vcf.gz.stats]);
rename('$name.vcf.gz.part','$name.vcf.gz') or Utils::error("rename $name.vcf.gz.part $name.vcf.gz: \$!");
Utils::CMD(qq[tabix -f -p vcf $name.vcf.gz]);
    ];
}


# Runs jobs in parallel for each BAM file, splitting each chromosome into chunks.
#   1) For each bam file, get a list of chromosomes and split them into chunks (calls $$opts{split_chunks}).
#   2) Then merge the parts into one VCF file and delete the intermediate chunks (calls $$opts{merge_chunks}).
#   3) Finally, create one big VCF file with all the data, each column corresponds to one bam file (calls $$opts{merge_vcf_files}).
#
# GATK and QCall run on populations and produce single VCF file straight away. Therefore the merge_chunks 
#   step is the final one and merge_vcf_files is not called on them.
#
#
# The parameters:
#   work_dir .. the working directory, should be of the form dir/(varFilter|mpileup|qcall|gatk)
#
# split_chunks function:
#   parameters .. $bam,$chunk_name,$chunk
#   must produce the file _$chunk_name.done, otherwise it will be run over and over again.
#
# merge_chunks function:
#   parameters .. $chunks,$name,$dir
#   must produce the file $name.vcf.gz out of the files from the $dir
#
# merge_vcf_files:
#   parameters .. $vcfs,$name
#   must produce the file $name.vcf.gz from the $vcfs 
#
sub run_in_parallel
{
    my ($self,$opts) = @_;

    if ( !$$self{fa_ref} ) { $self->throw("Missing the option fa_ref.\n"); }
    if ( !$$self{split_size} && !$$opts{split_size} ) { $self->throw("Missing the option split_size\n"); }

    if ( !$$opts{work_dir} ) { $self->throw("Missing the option work_dir\n"); }

    my $split_size = $$opts{split_size} ? $$opts{split_size} : $$self{split_size};
    my $chunks = $self->chr_chunks($split_size);
    if ( !@$chunks ) { $self->throw("No chunks to run? Please check your fai_chr_regex settings.\n"); }

    my $work_dir = $$opts{work_dir};
    Utils::CMD("mkdir -p $work_dir") unless -e $work_dir;

    my $is_finished = 1;

    my $prefix = $self->common_prefix($$opts{bams});
    my %names;
    my @to_be_merged;
    my $ntasks = 0;
    my $max_jobs = $$self{max_jobs} ? $$self{max_jobs} : 0;
    my $done = $LSF::Done;

    # Step 1, call the split_chunks command.
    for my $bam (@{$$opts{bams}})
    {
        # Take from each bam file name a unique part 
        my ($outdir,$name,$suffix) = Utils::basename($bam);
        $outdir =~ s{^/*$prefix/*}{};

        # Sanity check - is the result really unique? Should be, but does not hurt to check.
        if ( $names{"$outdir/$name"} ) { $self->throw("FIXME, not ready for this: Conflicting bam file names: $outdir/$name\n"); }
        $names{"$outdir/$name"} = 1;

        # Add to the list of VCF files to be merged in step 3.
        push @to_be_merged, "$work_dir/$outdir/$name.vcf.gz";

        # If the merged file from step 2 is already in place, no need for splitting. Check however, if there 
        #   are any temp files left and move them into the _done directory. (There can be thousands of files.)
        if ( $$self{fsu}->file_exists("$work_dir/$outdir/$name.vcf.gz") ) 
        {
            $self->clean_files("$work_dir/$outdir",$name,"$work_dir/$outdir/_done"); 
            for my $chunk (@$chunks)
            {
                my $chunk_name = $name.'_'.$chunk;
                unlink("$work_dir/$outdir/_$chunk_name.done") unless !$$self{fsu}->file_exists("$work_dir/$outdir/_$chunk_name.done");
            }
            next; 
        }

        # If we are here, the step 1 has not finished yet.
        $is_finished = 0;

        # Run split_chunks command for each chunk. (E.g. pileup for each chunks and samtools varfilter.)
        my $chunks_finished = 1;
        for my $chunk (@$chunks)
        {
            my $chunk_name = $name.'_'.$chunk;
            if ( $$self{fsu}->file_exists("$work_dir/$outdir/_$chunk_name.done") || $$self{fsu}->file_exists("$work_dir/$outdir/$name.vcf.gz") ) 
            { 
                $self->clean_files("$work_dir/$outdir",$chunk_name,"$work_dir/$outdir/_done");
                next; 
            }
            $chunks_finished = 0;
            Utils::CMD("mkdir -p $work_dir/$outdir") unless $$self{fsu}->file_exists("$work_dir/$outdir");

            my $jids_file = "$work_dir/$outdir/_$chunk_name.jid";
            my $status = LSF::is_job_running($jids_file);

            $ntasks++;
            if ( $max_jobs && $ntasks>$max_jobs )
            {
                $done |= $LSF::Running;
                $self->debug("The limit reached, $max_jobs jobs are running.\n");
                last;
            }

            if ( $status&$LSF::Running ) { $done |= $LSF::Running; next; }
            if ( $status&$LSF::Error ) 
            { 
                $done |= $LSF::Error;
                $self->warn("The command failed: $work_dir/$outdir .. perl -w $$self{prefix}$chunk_name.pl\n");
            }

            # Now get the splitting command
            my $cmd = &{$$opts{split_chunks}}($self,$bam,$chunk_name,$chunk);

            open(my $fh,'>',"$work_dir/$outdir/$$self{prefix}$chunk_name.pl") or $self->throw("$work_dir/$outdir/$$self{prefix}$chunk_name.pl: $!");
            print $fh $cmd;
            close($fh);
            LSF::run($jids_file,"$work_dir/$outdir","$$self{prefix}$chunk_name",$$opts{bsub_opts},qq[perl -w $$self{prefix}$chunk_name.pl]);
            $self->debug("Submitting $work_dir/$outdir .. perl -w $$self{prefix}$chunk_name.pl\n");
            $done |= $LSF::Running;
        }

        # Some of the chunks is not ready
        if ( !$chunks_finished ) { next; }

        # Now merge the chunks for one bam file into one VCF file.
        my $jids_file = "$work_dir/$outdir/_$name.jid";
        my $status = LSF::is_job_running($jids_file);
        if ( $status&$LSF::Running ) { $done |= $LSF::Running; next; }
        if ( $status&$LSF::Error ) 
        { 
            $done |= $LSF::Error;
            $self->warn("The command failed: $work_dir/$outdir .. perl -w $$self{prefix}$name.pl\n");
        }

        # Get the merging command
        my $cmd = &{$$opts{merge_chunks}}($self,$chunks,$name,"$work_dir/$outdir");

        open(my $fh,'>',"$work_dir/$outdir/_$name.pl") or $self->throw("$work_dir/$outdir/$$self{prefix}$name.pl: $!");
        print $fh $cmd;
        close($fh);
        LSF::run($jids_file,"$work_dir/$outdir","$$self{prefix}$name",$$opts{bsub_opts},qq[perl -w $$self{prefix}$name.pl]);
        $self->debug("Submitting $work_dir/$outdir .. $$self{prefix}$name\n");
    }

    if ( $done&$LSF::Error ) { $self->throw("$$opts{job_type} failed for some of the files.\n"); }
    if ( $done&$LSF::Running ) { $self->debug("Some files not finished...\n"); return $$self{No}; }


    # Some chunks still not done?
    if ( !$is_finished ) { return; }

    # The big VCF file from the step 3. already exists.
    if ( $$self{fsu}->file_exists("$work_dir/$$opts{job_type}.vcf.gz") ) 
    { 
        Utils::CMD("touch $$opts{dir}/$$opts{job_type}.done");
        return; 
    }

    if ( !$$opts{merge_vcf_files} ) 
    {
        # GATK and QCall run on multiple BAMs and create a single VCF output file, no need to merge. 
        if ( @to_be_merged != 1 ) { $self->throw("FIXME: [",join(' ',@to_be_merged), "] --> $work_dir/$$opts{job_type}.vcf.gz"); }
        rename($to_be_merged[0],"$work_dir/$$opts{job_type}.vcf.gz") 
            or $self->throw("rename $to_be_merged[0] $work_dir/$$opts{job_type}.vcf.gz: $!");
        rename("$to_be_merged[0].tbi","$work_dir/$$opts{job_type}.vcf.gz.tbi") 
            or $self->throw("rename $to_be_merged[0].tbi $work_dir/$$opts{job_type}.vcf.gz.tbi: $!");
        rename("$to_be_merged[0].stats","$work_dir/$$opts{job_type}.vcf.gz.stats") 
            or $self->throw("rename $to_be_merged[0].stats $work_dir/$$opts{job_type}.vcf.gz.stats: $!");
        Utils::CMD("touch $$opts{dir}/$$opts{job_type}.done");
        return;
    }


    # Because this subroutine returns as if it has already finished, a custom jids_file must
    #   be used: Pipeline.pm will delete the $lock_file.
    my $jids_file = "$work_dir/_$$opts{job_type}_merge.jid";
    my $status = LSF::is_job_running($jids_file);
    if ( $status&$LSF::Running ) { return; }
    if ( $status&$LSF::Error ) { $self->warn("Some jobs failed: $jids_file\n"); }

    # Get the VCF merging command
    my $cmd = &{$$opts{merge_vcf_files}}($self,\@to_be_merged,$$opts{job_type});
    open(my $fh,'>',"$work_dir/_merge.pl") or $self->throw("$work_dir/_merge.pl: $!");
    print $fh $cmd;
    close($fh);
    LSF::run($jids_file,$work_dir,"_merge",$$opts{bsub_opts},qq[perl -w _merge.pl]);
    $self->debug("Submitting $work_dir .. _merge\n");
    
    return;
}

sub dump_opts
{
    my ($self,@keys) = @_;
    my %opts;
    for my $key (@keys)
    {
        $opts{$key} = exists($$self{$key}) ? $$self{$key} : undef;
    }
    return Data::Dumper->Dump([\%opts],["opts"]);
}


sub glue_vcf_chunks
{
    my ($self,$chunks,$name,$work_dir) = @_;

    open(my $fh,'>',"$work_dir/$name.chunks.list") or $self->throw("$work_dir/$name.chunks.list: $!");
    for my $chunk (@$chunks)
    {
        print $fh "${name}_$chunk.vcf.gz\n";
    }
    close($fh);

    my $out = qq[
use strict;
use warnings;
use Utils;

Utils::CMD("rm -f $name.vcf-tmp.gz.part");
Utils::CMD("vcf-concat -s 3 -f $name.chunks.list | $$self{vcf_rmdup} | bgzip -c > $name.vcf.gz.part");
Utils::CMD(qq[zcat $name.vcf.gz.part | $$self{vcf_stats} > $name.vcf.gz.stats]);
Utils::CMD(qq[tabix -f -p vcf $name.vcf.gz.part]);
rename('$name.vcf.gz.part.tbi','$name.vcf.gz.tbi') or Utils::error("rename $name.vcf.gz.part.tbi $name.vcf.gz.tbi: \$!");
rename('$name.vcf.gz.part','$name.vcf.gz') or Utils::error("rename $name.vcf.gz.part $name.vcf.gz: \$!");

];

    for my $chunk (@$chunks)
    {
        $out .= qq[unlink('${name}_$chunk.vcf.gz');\n]; 
        $out .= qq[unlink('${name}_$chunk.vcf.gz.tbi');\n]; 
        $out .= qq[unlink('$$self{prefix}$chunk.names');\n]; 
        $$self{fsu}->file_exists("$work_dir/${name}_$chunk.vcf.gz",wipe_out=>1);
    }

    return $out;
}



#---------- varFilter ---------------------


# Requires the bam files listed in the file
sub varFilter_requires
{
    my ($self,$dir) = @_;
    return [$$self{file_list}];
}


sub varFilter_provides
{
    my ($self,$dir) = @_;
    my @provides = ('varFilter.done');
    return \@provides;
}

sub varFilter
{
    my ($self,$dir,$lock_file) = @_;

    if ( !$$self{varfilter} ) { $self->throw("Missing the option varfilter.\n"); }
    if ( !$$self{sam2vcf} ) { $self->throw("Missing the option sam2vcf.\n"); }

    my $bams = $self->read_files($$self{file_list});
    my %opts = 
    (
        dir        => $dir,
        work_dir   => "$dir/varFilter",
        job_type   => 'varFilter',
        bsub_opts  => {bsub_opts=>$$self{bsub_opts_varfilter},dont_wait=>1,append=>0},
        bams       => $bams,
        split_size => $$self{split_size_varfilter},

        split_chunks    => \&varfilter_split_chunks,    
        merge_chunks    => \&varfilter_merge_chunks,
        merge_vcf_files => \&merge_vcf_files,
    );

    $self->run_in_parallel(\%opts);

    return $$self{Yes};
}

sub run_varfilter
{
    my ($self,$bam,$name,$chunk) = @_;
    if ( !-e "$name.pileup" )
    {
        Utils::CMD(qq[samtools view -bh $bam $chunk | samtools pileup $$self{samtools_pileup_params} -c -f $$self{fa_ref} - | $$self{varfilter} > $name.pileup.part],{verbose=>1});
        rename("$name.pileup.part","$name.pileup") or $self->throw("rename $name.pileup.part $name.pileup: $!");
    }
    Utils::CMD("touch _$name.done",{verbose=>1});
}


sub varfilter_split_chunks
{
    my ($self,$bam,$chunk_name,$chunk) = @_;

    return qq[
use strict;
use warnings;
use VertRes::Pipelines::SNPs;

my \$opts = {
    file_list  => q[$$self{file_list}],
    fa_ref     => q[$$self{fa_ref}],
    fai_ref    => q[$$self{fai_ref}],
    sam2vcf    => q[$$self{sam2vcf}],
    varfilter  => q[$$self{varfilter}],
    samtools_pileup_params => q[$$self{samtools_pileup_params}],
};
my \$snps = VertRes::Pipelines::SNPs->new(%\$opts);
\$snps->run_varfilter(q[$bam],q[$chunk_name],q[$chunk]);

    ];
}

sub varfilter_merge_chunks
{
    my ($self,$chunks,$name,$work_dir) = @_;

    open(my $fh,'>',"$work_dir/$name.chunks.list") or $self->throw("$work_dir/$name.chunks.list: $!");
    for my $chunk (@$chunks)
    {
        my $chunk_name = $name.'_'.$chunk;
        print $fh "$chunk_name.pileup\n";
    }
    close($fh);

    return qq[
use strict;
use warnings;
use Utils;

if ( !-e "$name.pileup.gz" )
{
    # Make sure xargs does not split the arguments to multiple commands: check using `sort -c`
    Utils::CMD("cat $name.chunks.list | xargs sort -k1,1 -k2,2n -m | sort -c -k1,1 -k2,2n | bgzip -c > $name.pileup.gz.part");
    rename("$name.pileup.gz.part","$name.pileup.gz") or Utils::error("rename $name.pileup.gz.part $name.pileup.gz: \$!");
    Utils::CMD("cat $name.chunks.list | xargs rm -f");
    Utils::CMD("rm -f $name.chunks.list");
}
if ( !-e "$name.vcf.gz" )
{
    Utils::CMD("zcat $name.pileup.gz | $$self{pileup_rmdup} | $$self{sam2vcf} -s -t $name | bgzip -c > $name.vcf.gz.part",{verbose=>1});
    Utils::CMD(qq[zcat $name.vcf.gz.part | $$self{vcf_stats} > $name.vcf.gz.stats]);
    rename("$name.vcf.gz.part","$name.vcf.gz") or Utils::error("rename $name.vcf.gz.part $name.vcf.gz: \$!");
    Utils::CMD(qq[tabix -f -p vcf $name.vcf.gz]);
}

    ];
}

#---------- pseudo_genome ---------------------


# Requires the bam files listed in the file
sub pseudo_genome_requires
{
    my ($self,$dir) = @_;
    return [$$self{file_list}];
}


sub pseudo_genome_provides
{
    my ($self,$dir) = @_;
    my @provides = ('.pseudo_genome.done');
    return \@provides;
}

# Create a pseudo genome for the organism
sub pseudo_genome {
    my ($self,$lane_path,$action_lock) = @_;

    #$self->{vrtrack} = VRTrack::VRTrack->new($self->{db}) or $self->throw("Could not connect to the database\n");
    #my $vlane = VRTrack::Lane->new_by_name($self->{vrtrack}, $self->{lane}); 
    
    if ( !$$self{fa_ref} ) { $self->throw("Missing the option fa_ref.\n"); }
    
    #There will be a single BAM file in the list file ( $$self{file_list} )  
    my $bam_in_array_ref = $self->read_files($$self{file_list});
    my $single_bam_file = $bam_in_array_ref->[0];
    
    my $memory = $self->{memory};
    if (! defined $memory || $memory < 2000) {
        $memory = 2000;
    }
    
    
    # Script to be run by LSF
    my $script_name = $self->{fsu}->catfile($lane_path, $self->{prefix}."pseudo_genome.pl");
    open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
    
    print $scriptfh qq[
use strict;
use warnings;
use Log::Log4perl qw(get_logger);
use Pathogens::Variant::Utils::PseudoReferenceMaker;
use Pathogens::Variant::Exception qw(:try);

Log::Log4perl->init(\\ qq{
    log4perl.logger = INFO, AppInfo, AppError

    # Filter to match level WARN
    log4perl.filter.MatchError = Log::Log4perl::Filter::LevelMatch
    log4perl.filter.MatchError.LevelToMatch  = WARN
    log4perl.filter.MatchError.AcceptOnMatch = true

    # Filter to match level INFO
    log4perl.filter.MatchInfo  = Log::Log4perl::Filter::LevelMatch
    log4perl.filter.MatchInfo.LevelToMatch  = INFO
    log4perl.filter.MatchInfo.AcceptOnMatch = true

    # Error appender
    log4perl.appender.AppError = Log::Log4perl::Appender::Screen
    log4perl.appender.AppError.stderr   = 1
    log4perl.appender.AppError.layout   = SimpleLayout
    log4perl.appender.AppError.Filter   = MatchError

    # Info appender
    log4perl.appender.AppInfo = Log::Log4perl::Appender::Screen
    log4perl.appender.AppInfo.stderr   = 0
    log4perl.appender.AppInfo.layout   = SimpleLayout
    log4perl.appender.AppInfo.Filter   = MatchInfo
});
my \$logger = get_logger();

my \%args = (
     out       => 'pseudo_genome.fasta'
   , bam       => '$single_bam_file'
   , reference => '$$self{fa_ref}'
   , lane_name => '$self->{lane}'
);

my \$pseudo_reference_maker;
try {
   \$pseudo_reference_maker = Pathogens::Variant::Utils::PseudoReferenceMaker->new(arguments => \\%args);
   \$pseudo_reference_maker->create_pseudo_reference();
   \$logger->info(\$pseudo_reference_maker->get_statistics_dump);

} catch Error with {
    my \$E = shift;
    \$logger->error(\$E->stacktrace);

    \$pseudo_reference_maker->_clean_up_the_temporary_files;
    unlink 'pseudo_genome.fasta' if -e 'pseudo_genome.fasta';
};

system('touch .pseudo_genome.done');

exit 0;
];
    close $scriptfh;
    
    my $job_name = $self->{prefix}.'pseudo_genome';
    $self->archive_bsub_files($lane_path, $job_name);
    LSF::run($action_lock, $lane_path, $job_name, {bsub_opts => "-q normal -M${memory}000 -R 'select[mem>$memory] rusage[mem=$memory]'", dont_wait=>1 }, qq{perl -w $script_name});
    
    # we've only submitted to LSF, so it won't have finished; we always return
    # that we didn't complete
    return $self->{No};
}


#---------- mpileup ---------------------


# Requires the bam files listed in the file
sub mpileup_requires
{
    my ($self,$dir) = @_;
    return [$$self{file_list}];
}


sub mpileup_provides
{
    my ($self,$dir) = @_;
    my @provides = ('mpileup.done');
    return \@provides;
}

sub mpileup
{
    my ($self,$dir,$lock_file) = @_;

    if ( !$$self{mpileup_cmd} ) { $self->throw("Missing the option mpileup_cmd.\n"); }
    if ( !$$self{bcftools} ) { $self->throw("Missing the option bcftools.\n"); }

    my $bams = [ $$self{file_list} ];
    my %opts = 
    (
        dir        => $dir,
        work_dir   => "$dir/mpileup",
        job_type   => 'mpileup',
        bsub_opts  => {bsub_opts=>$$self{bsub_opts_mpileup},dont_wait=>1,append=>0},
        bams       => $bams,
        split_size => $$self{split_size_mpileup},
        bcf_based  => $$self{bcf_based},

        split_chunks    => \&mpileup_split_chunks,    
        merge_chunks    => \&glue_vcf_chunks,
        merge_vcf_files => $$self{filter4vcf} ? \&mpileup_postprocess : undef,
    );
    if ( exists($$self{mpileup_samples}) ) { $opts{mpileup_samples} = $$self{mpileup_samples}; }
    if ( $$self{bcf_based} )
    {
        $opts{merge_chunks}    = \&glue_bcf_chunks;
        $opts{merge_vcf_files} = \&merge_bcf_files;
    }

    $self->run_in_parallel(\%opts);

    return $$self{Yes};
}

sub glue_bcf_chunks
{
    my ($self,$chunks,$name,$work_dir) = @_;
    return qq[ `touch $name.vcf.gz` ];
}

sub merge_bcf_files
{
    my ($self,$vcfs,$name) = @_;
    return qq[ `touch $name.vcf.gz` ];
}

sub run_mpileup
{
    my ($self,$file_list,$name,$chunk) = @_;

    if ( $$self{bcf_based} )
    {
        # Produce BCFs instead of VCFs
        if ( ! -e "$name.bcf" )
        {
            Utils::CMD(qq[$$self{mpileup_cmd} -g -b $file_list -r $chunk -f $$self{fa_ref} > $name.bcf.part],{verbose=>1});
            rename("$name.bcf.part","$name.bcf") or $self->throw("rename $name.bcf.part $name.bcf: $!");
        }

        Utils::CMD("touch _$name.done",{verbose=>1});
        return;
    }
    
    if ( -e "$name.bcf" && exists($$self{mpileup_samples}) )
    {
        # Take existing BCF and create a population subsample
        if ( ! -e "$name.pop.bcf" ) 
        {
            Utils::CMD(qq[$$self{bcftools} view -p 0.99 -bvcg -s $$self{mpileup_samples} $name.bcf > $name.pop.bcf.part],{verbose=>1});
            rename("$name.pop.bcf.part","$name.pop.bcf") or $self->throw("$name.pop.bcf.part $name.pop.bcf: $!");
        }
        if ( ! -e "$name.vcf.gz" )
        {
            Utils::CMD(qq[$$self{bcftools} view $name.pop.bcf | bgzip -c > $name.vcf.gz.part],{verbose=>1});
            Utils::CMD(qq[tabix -f -p vcf $name.vcf.gz.part],{verbose=>1});
            rename("$name.vcf.gz.part.tbi","$name.vcf.gz.tbi") or $self->throw("rename $name.vcf.gz.part.tbi $name.vcf.gz.tbi: $!");
            rename("$name.vcf.gz.part","$name.vcf.gz") or $self->throw("rename $name.vcf.gz.part $name.vcf.gz: $!");
        }
    }
    else
    {
        # Original code which does not save BCFs and pipes mpileup directly through bcftools to produce VCFs 
        if ( ! -e "$name.vcf.gz" )
        {
            Utils::CMD(qq[$$self{mpileup_cmd} -b $file_list -r $chunk -f $$self{fa_ref} | $$self{bcftools} view -p 0.99 -gcv - | bgzip -c > $name.vcf.gz.part],{verbose=>1});
            Utils::CMD(qq[tabix -f -p vcf $name.vcf.gz.part],{verbose=>1});
            rename("$name.vcf.gz.part.tbi","$name.vcf.gz.tbi") or $self->throw("rename $name.vcf.gz.part.tbi $name.vcf.gz.tbi: $!");
            rename("$name.vcf.gz.part","$name.vcf.gz") or $self->throw("rename $name.vcf.gz.part $name.vcf.gz: $!");
        }
    }
    Utils::CMD("touch _$name.done",{verbose=>1});
}


sub mpileup_split_chunks
{
    my ($self,$bam,$chunk_name,$chunk) = @_;

    my $opts = $self->dump_opts(qw(file_list fa_ref fai_ref mpileup_cmd mpileup_samples bcftools bcf_based filter4vcf));

    return qq[
use strict;
use warnings;
use VertRes::Pipelines::SNPs;

my $opts

my \$snps = VertRes::Pipelines::SNPs->new(%\$opts);
\$snps->run_mpileup(q[$bam],q[$chunk_name],q[$chunk]);

    ];
}


# Apply the vcfutils.pl filter4vcf
sub mpileup_postprocess
{
    my ($self,$vcfs,$name) = @_;

    if ( scalar @$vcfs > 1 ) { $self->throw("FIXME: expected one file only\n"); }

    my $basename = $$vcfs[0];
    $basename =~ s/\.gz$//i;
    $basename =~ s/\.vcf$//i;

    return qq[
use strict;
use warnings;
use Utils;

if ( -e "$basename.vcf.gz" )
{
    rename("$basename.vcf.gz.tbi","$name.unfilt.vcf.gz.tbi") or Utils::error("rename $basename.vcf.gz.tbi $name.unfilt.vcf.gz.tbi: \$!");
    rename("$basename.vcf.gz.stats","$name.unfilt.vcf.gz.stats") or Utils::error("rename $basename.vcf.gz.stats $name.unfilt.vcf.gz.stats: \$!");
    rename("$basename.vcf.gz","$name.unfilt.vcf.gz") or Utils::error("rename $basename.vcf.gz $name.unfilt.vcf.gz: \$!");
}

Utils::CMD("zcat $name.unfilt.vcf.gz | $$self{filter4vcf} | bgzip -c > $basename.filt.vcf.gz");
Utils::CMD("zcat $basename.filt.vcf.gz | $$self{vcf_stats} > $name.vcf.gz.stats");
Utils::CMD("tabix -f -p vcf $basename.filt.vcf.gz");
rename("$basename.filt.vcf.gz.tbi","$name.vcf.gz.tbi") or Utils::error("rename $basename.filt.vcf.gz.tbi $name.vcf.gz.tbi: \$!");
rename("$basename.filt.vcf.gz","$name.vcf.gz") or Utils::error("rename $basename.filt.vcf.gz $name.vcf.gz: \$!");

];
}


#---------- gatk ---------------------


# Requires the bam files listed in the file
sub gatk_requires
{
    my ($self,$dir) = @_;
    return [$$self{file_list}];
}


sub gatk_provides
{
    my ($self,$dir) = @_;
    my @provides = ('gatk.done');
    return \@provides;
}

sub gatk
{
    my ($self,$dir,$lock_file) = @_;

    if ( !$$self{fa_ref} ) { $self->throw("Missing the option fa_ref.\n"); }
    if ( !$$self{fai_ref} ) { $self->throw("Missing the option fai_ref.\n"); }
    if ( !$$self{gatk_opts} ) { $self->throw("Missing the option gatk_opts.\n"); }

    if ( !($$self{file_list}=~/\.list$/) ) { $self->throw("GATK requires .list extension for the file of file names: $$self{file_list}\n"); }

    if ( !exists($$self{gatk_opts}{all}) ) { $$self{gatk_opts}{all}={}; }
    my $gopts = $$self{gatk_opts}{all};

    # dbSNP options - bit of a mess.
    # if gatk_opts{all}{dbsnp} is set, use that (dbSNP rod file)
    # if dbSNP_rod is set, set dbsnp to that (dbSNP rod file)
    # if dbSNP_vcf is set, set_b for that (dbSNP vcf file)
    # --DBSNP is deprecated, but keep dbsnp & dbSNP_rod around for those rod files.
    # should now be using vcf format dbsnp files, and specifying with dbSNP_vcf
    #  
    # Use $self->{have_dbSNP} as flag to pass around that we're using dbSNP (used to use dbSNP_rod)

    $$self{have_dbSNP} = 0;
    # old style dbSNP rod file
    if ( exists($$gopts{dbsnp})){
        $$self{have_dbSNP} = 1;
        if ($$gopts{dbsnp} ne $$self{dbSNP_rod} ) {
            $self->throw("Conflicting options for dbSNP, which one to use: $$gopts{dbsnp} or $$self{dbSNP_rod}?\n");
        }
    }
    if ( !exists($$gopts{dbsnp}) )
    {
        if ( exists($$self{dbSNP_rod}) ) { 
            $$gopts{dbsnp} = $$self{dbSNP_rod}; 
            $$self{have_dbSNP} = 1;
        }
        else { 
            $$gopts{dbsnp} = undef; # undef the default dbSNP rod file explicitly
        }  
    }

    # new style vcf dbSNP file
    if ( exists($$self{dbSNP_vcf}) && $$self{dbSNP_vcf}){
        push @{$gopts->{bs}}, "dbsnp,VCF $$self{dbSNP_vcf}";
        # override any old rod file
        $$gopts{dbsnp} = undef;
        $$self{have_dbSNP} = 1;
    }

    if ( exists($$gopts{reference}) && $$gopts{reference} ne $$self{fa_ref} )
    {
        $self->throw("Conflicting options for reference, which one to use: $$gopts{reference} or $$self{fa_ref}?\n");
    }
    if ( !exists($$gopts{reference}) ) { $$gopts{reference} = $$self{fa_ref}; }
    my $bams = [ $$self{file_list} ];
    my %opts =
    (
        dir        => $dir,
        work_dir   => "$dir/gatk",
        job_type   => 'gatk',
        bsub_opts  => {bsub_opts=>$$self{bsub_opts_gatk},dont_wait=>1,append=>0},
        split_size => $$self{split_size_gatk},
        bams       => $bams,
        gatk_opts  => $$self{gatk_opts},

        split_chunks    => \&gatk_split_chunks,
        merge_chunks    => \&glue_vcf_chunks,
        merge_vcf_files => \&gatk_postprocess,
    );

    $self->run_in_parallel(\%opts);

    return $$self{Yes};
}


sub gatk_split_chunks
{
    my ($self,$bam,$chunk_name,$chunk) = @_;

    my $opts = $self->dump_opts(qw(file_list fa_ref fai_ref sam2vcf gatk_opts have_dbSNP indel_mask));

    return qq[
use strict;
use warnings;
use VertRes::Pipelines::SNPs;

my $opts

my \$snps = VertRes::Pipelines::SNPs->new(%\$opts);
\$snps->run_gatk_chunk(q[$bam],q[$chunk_name],q[$chunk]);
    ];
}


# Set the hash for GATK wrapper. $$self{gatk_opts} is set in the config file and
#   $keys is an array which contains a list of subkeys to be used, such as 
#   qw(all unified_genotyper).
#
# and
#   gatk_opts =>
#   {
#       all =>
#       {
#           _extras => [ '-U ALLOW_UNSET_BAM_SORT_ORDER' ],
#       }
#       unified_genotyper => 
#       {
#           standard_min_confidence_threshold_for_calling => 10,
#       }
#   }
#
sub get_gatk_opts
{
    my ($self,@keys) = @_;
    my %out = ();
    for my $key (@keys)
    {
        if ( !exists($$self{gatk_opts}) or !exists($$self{gatk_opts}{$key}) ) { next; }
        %out = ( %out, %{$$self{gatk_opts}{$key}} );
    }
    return %out;
}

sub run_gatk_chunk
{
    my ($self,$bam,$name,$chunk) = @_;
    
    my $gatk;
    my %opts = $self->get_gatk_opts(qw(all unified_genotyper));
    
    $gatk = VertRes::Wrapper::GATK->new(%opts);
    $gatk->unified_genotyper($bam,"$name.vcf.gz",L=>$chunk,%opts);
    Utils::CMD("tabix -f -p vcf $name.vcf.gz");
    
    if ( !$$self{indel_mask} or !-e $$self{indel_mask} )
    {
        %opts = $self->get_gatk_opts(qw(all indel_genotyper));
        $gatk = VertRes::Wrapper::GATK->new(%opts);
        $gatk->somatic_indel_detector($bam, "$name.indels.raw.bed", "$name.indels.detailed.bed",L=>$chunk,%opts);
        Utils::CMD("make_indel_mask.pl $name.indels.raw.bed 10 $name.indels.mask.bed",{verbose=>1});
    
        %opts = $self->get_gatk_opts(qw(all variant_filtration));
        $gatk = VertRes::Wrapper::GATK->new(%opts);
        $gatk->set_b("variant,VCF,$name.vcf.gz","mask,Bed,$name.indels.mask.bed");
        $gatk->variant_filtration("$name.filtered.vcf.gz",%opts);
    
        unlink("$name.vcf.gz");
        unlink("$name.vcf.gz.tbi");
        unlink("$name.indels.detailed.bed");
        unlink("$name.indels.detailed.bed.idx");
        unlink("$name.indels.mask.bed");
        unlink("$name.indels.mask.bed.idx");
        unlink("$name.indels.raw.bed");
    
        Utils::CMD("tabix -f -p vcf $name.filtered.vcf.gz");
        rename("$name.filtered.vcf.gz.tbi","$name.vcf.gz.tbi") or $self->throw("rename $name.filtered.vcf.gz.tbi $name.vcf.gz.tbi: $!");
        rename("$name.filtered.vcf.gz","$name.vcf.gz") or $self->throw("rename $name.filtered.vcf.gz $name.vcf.gz: $!");
    }
    
    Utils::CMD("touch _$name.done",{verbose=>1});
}


sub gatk_postprocess
{
    my ($self,$vcfs,$name) = @_;

    if ( scalar @$vcfs > 1 ) { $self->throw("FIXME: expected one file only\n"); }

    my $opts = $self->dump_opts(qw(file_list fa_ref fai_ref sam2vcf gatk_opts have_dbSNP indel_mask));

    return qq[
use strict;
use warnings;
use VertRes::Pipelines::SNPs;

my $opts

my \$snps = VertRes::Pipelines::SNPs->new(%\$opts);
\$snps->run_gatk_postprocess(q[$$vcfs[0]],q[$name]);
];
}


sub run_gatk_postprocess
{
    my ($self,$vcf,$name) = @_;

    my $basename = $vcf;
    $basename =~ s/\.gz$//i;
    $basename =~ s/\.vcf$//i;

    my (%opts,$gatk);

    # Filter
    if ( $$self{indel_mask} && -e $$self{indel_mask} )
    {
        %opts = $self->get_gatk_opts(qw(all variant_filtration));
        $gatk = VertRes::Wrapper::GATK->new(%opts);
        $gatk->set_b("variant,VCF,$vcf", "mask,Bed,$$self{indel_mask}");
        $gatk->variant_filtration("$basename.filtered.vcf.gz", %opts);
        rename("$basename.filtered.vcf.gz","$basename.tmp.vcf.gz") or $self->throw("rename $basename.filtered.vcf.gz $basename.tmp.vcf.gz: $!");
        Utils::CMD("tabix -f -p vcf $basename.tmp.vcf.gz");
    }
    else
    {
        Utils::CMD("cp $basename.vcf.gz $basename.tmp.vcf.gz");
        Utils::CMD("tabix -f -p vcf $basename.tmp.vcf.gz");
    }
    
    if ( $$self{have_dbSNP} )
    {
        # Recalibrate
        %opts = $self->get_gatk_opts(qw(all variant_recalibrator));
        $gatk = VertRes::Wrapper::GATK->new(%opts);
        $gatk->set_b("input,VCF,$basename.tmp.vcf.gz");
        $gatk->add_b(@{$opts{training_sets}});
        $gatk->set_annotations('HaplotypeScore', 'MQRankSum', 'ReadPosRankSum', 'HRun');
        $gatk->variant_recalibrator("$basename.recal", "$basename.tranches", %opts);
        
        # Filter the calls down to an implied 10.0% novel false discovery rate
        %opts = $self->get_gatk_opts(qw(all apply_recalibration));
        $gatk = VertRes::Wrapper::GATK->new(%opts);
        $gatk->set_b("input,VCF,$basename.tmp.vcf.gz");
        Utils::CMD("tabix -f -p vcf $basename.tmp.vcf.gz");
        $gatk->apply_recalibration("$basename.recal", "$basename.tranches","$basename.recalibrated.vcf.gz", %opts);
        rename("$basename.recalibrated.vcf.gz","$basename.tmp.vcf.gz") or $self->throw("rename $basename.recalibrated.vcf.gz $basename.tmp.vcf.gz: $!");
        Utils::CMD("tabix -f -p vcf $basename.tmp.vcf.gz");
    }

    Utils::CMD(qq[zcat $basename.tmp.vcf.gz | $$self{vcf_stats} > $name.vcf.gz.stats]);
    rename("$basename.tmp.vcf.gz.tbi","$name.vcf.gz.tbi") or Utils::error("rename $basename.tmp.vcf.gz.tbi $name.vcf.gz.tbi: \$!");
    rename("$basename.tmp.vcf.gz","$name.vcf.gz") or Utils::error("rename $basename.tmp.vcf.gz $name.vcf.gz: \$!");
}




#---------- qcall ---------------------


# Requires the bam files listed in the file
sub qcall_requires
{
    my ($self,$dir) = @_;
    return [$$self{file_list}];
}


sub qcall_provides
{
    my ($self,$dir) = @_;
    my @provides = ('qcall.done');
    return \@provides;
}

sub qcall
{
    my ($self,$dir,$lock_file) = @_;

    if ( !$$self{fa_ref} ) { $self->throw("Missing the option fa_ref.\n"); }
    if ( !$$self{fai_ref} ) { $self->throw("Missing the option fai_ref.\n"); }
    if ( !$$self{qcall_cmd} ) { $self->throw("Missing the option qcall_cmd.\n"); }
    if ( !$$self{mpileup_cmd} ) { $self->throw("Missing the option mpileup_cmd.\n"); }

    my $bams = [ $$self{file_list} ];
    my %opts =
    (
        dir        => $dir,
        work_dir   => "$dir/qcall",
        job_type   => 'qcall',
        bsub_opts  => {bsub_opts=>$$self{bsub_opts_qcall},dont_wait=>1,append=>0},
        bams       => $bams,
        split_size => $$self{split_size_qcall},

        split_chunks    => \&qcall_split_chunks,
        merge_chunks    => \&glue_vcf_chunks,
        merge_vcf_files => undef,
    );

    $self->run_in_parallel(\%opts);

    return $$self{Yes};
}

sub qcall_split_chunks
{
    my ($self,$bam,$chunk_name,$chunk) = @_;

    my $opts = $self->dump_opts(qw(file_list fa_ref fai_ref mpileup_cmd split_size qcall_cmd sort_cmd prefix));

    return qq[
use strict;
use warnings;
use VertRes::Pipelines::SNPs;

my $opts

my \$var = VertRes::Pipelines::SNPs->new(%\$opts);
\$var->run_qcall_chunk(q[$bam],q[$chunk_name],q[$chunk]);
    ];
}


sub run_qcall_chunk
{
    my ($self,$file_list,$chunk_name,$chunk) = @_;

    # QCall needs the list of samples in advance: open the BAM files,
    #   get the sample names and put them in a file.
    my $files  = $self->read_files($$self{file_list});
    my $groups = $self->bam_file_groups($files);
    open(my $fh,'>',"$$self{prefix}$chunk.names") or $self->throw("$$self{prefix}$chunk.names: $!");
    for my $grp (sort keys %$groups)
    {
        print $fh "$grp\n";
    }
    close($fh);

    if ( scalar keys %$groups < 2 ) { $self->throw("QCall is a population based caller and must be run with at least 2 samples.\n"); }

    # Execute mpileup with QCall
    Utils::CMD(qq[samtools mpileup -b $file_list -D -g -r $chunk -f $$self{fa_ref} | $$self{bcftools} view -Q - | $$self{qcall_cmd} -sn $$self{prefix}$chunk.names -co $chunk.vcf.part],{verbose=>1});

    # Compress the VCF file
    Utils::CMD("cat $chunk.vcf.part | bgzip -c > $chunk.vcf.gz.part",{verbose=>1});
    unlink("$chunk.vcf.part");
    Utils::CMD("tabix -p vcf -f $chunk.vcf.gz.part");

    # Clean
    unlink("$$self{prefix}$chunk.names");
    rename("$chunk.vcf.gz.part.tbi","$chunk_name.vcf.gz.tbi") or $self->throw("rename $chunk.vcf.gz.part.tbi $chunk_name.vcf.gz.tbi: $!");
    rename("$chunk.vcf.gz.part","$chunk_name.vcf.gz") or $self->throw("rename $chunk.vcf.gz.part $chunk_name.vcf.gz: $!");
    Utils::CMD("touch _$chunk_name.done");
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
    my @requires;

    # need a done file from each caller
    for my $task (keys %{$self->{task}}) {
        next if ($task eq "update_db" or $task eq "cleanup");

        if ($task eq "pseudo_genome") {
            #adding a "." prefix, otherwise "cleanup" would delete this file
            push @requires, ".$task.done";
            next;
        }

        push @requires, "$task.done";
    }   
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
    my @provides;

    if ($self->{task}{cleanup}){
        push @provides, '.snps_done';
    }
    else {
        @provides = @{$self->cleanup_requires($lane_path)};
        @provides or $self->throw("Something went wrong; cleanup action doesn't seem to require any files");
    }

    return \@provides;
}


=head2 cleanup

 Title   : cleanup
 Usage   : $obj->cleanup('/path/to/lane', 'lock_filename');
 Function: Cleans up files made by pipeline
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub cleanup {
    my ($self, $lane_path, $action_lock) = @_;
    my $fs = VertRes::Utils::FileSystem->new();

    # for each caller, we move the files we want to keep to the root directory
    # and the delete that caller's directory
    for my $task (keys %{$self->{task}}) {
        next if ($task eq "update_db" or $task eq "cleanup" or $task eq "pseudo_genome");
         
        my @suffixes = qw(vcf.gz vcf.gz.stats vcf.gz.tbi);

        # keep the unfiltered files from mpileup as well as filtered
        if ($task eq 'mpileup') {
            push @suffixes, qw(unfilt.vcf.gz unfilt.vcf.gz.stats unfilt.vcf.gz.tbi);
        }

        # move the files
        for my $suff (@suffixes) {
            my $old_file = File::Spec->catfile($lane_path, $task, "$task.$suff");
            my $new_file = File::Spec->catfile($lane_path, "$task.$suff");
            next unless (-s $old_file);
            $fs->copy($old_file, $new_file) or $self->throw("Error copying files:\nold: $old_file\nnew: $new_file");
        }
        
        if ($task eq 'gatk') {
            my $tranches = File::Spec->catfile($lane_path, $task, "file.tranches");
            if (-s $tranches) {
                $fs->copy($tranches, File::Spec->catfile($lane_path, "gatk.tranches")) or $self->throw("Error copying files:\nold: $tranches\n");
            }
        }

        # cleanup files in the root of the snps output directory
        Utils::CMD("rm " . File::Spec->catfile($lane_path, "$task.done"));

        # delete the caller directory
        Utils::CMD("rm -r " . File::Spec->catdir($lane_path, $task));
    } 

    unless ($self->{task}{update_db})
    {
        my $job_status =  File::Spec->catfile($lane_path, $self->{prefix} . 'job_status');
        Utils::CMD("rm $job_status") if (-e $job_status);
    }
    
    my $file_list = File::Spec->catfile($lane_path, 'file.list');
    Utils::CMD("rm $file_list") if (-e $file_list);
    Utils::CMD("touch " . File::Spec->catfile($lane_path, '.snps_done'));
    return $$self{'Yes'};
}


=head2 update_db_requires

 Title   : update_db_requires
 Usage   : my $required_files = $obj->update_db_requires('/path/to/lane');
 Function: Find out what files the update_db action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub update_db_requires {
    my ($self, $lane_path) = @_;
    my @requires = @{$self->cleanup_provides($lane_path)};
    @requires or $self->throw("Something went wrong; update_db action doesn't seem to require any files");
    return \@requires;
}

=head2 update_db_provides

 Title   : update_db_provides
 Usage   : my $provided_files = $obj->update_db_provides('/path/to/lane');
 Function: Find out what files the update_db action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub update_db_provides {
    return [];
}

=head2 update_db

 Title   : update_db
 Usage   : $obj->update_db('/path/to/lane', 'lock_filename');
 Function: Records in the database that the lane has been improved.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut



sub update_db {
    my ($self, $lane_path, $action_lock) = @_;

    # all the done files are there, so just need to update the processed
    # flag in the database (if not already updated)
    my $vrlane = $self->{vrlane};
    my $vrtrack = $vrlane->vrtrack;

    # Ignore snp_called status. Used for lanes with bam files from different mappers/assemblies.
    my $ignore_snp_called = (exists $$self{'ignore_snp_called_status'} && $$self{'ignore_snp_called_status'}) ? 1:0;

    return $$self{'Yes'} if $vrlane->is_processed('snp_called') && !$ignore_snp_called;

    unless($vrlane->is_processed('snp_called')){
	$vrtrack->transaction_start();
	$vrlane->is_processed('snp_called',1);
	$vrlane->update() || $self->throw("Unable to set improved status on lane $lane_path");
	$vrtrack->transaction_commit();
    }

    if ($self->{task}{cleanup}) {
        my $job_status =  File::Spec->catfile($lane_path, $self->{prefix} . 'job_status');
        Utils::CMD("rm $job_status") if (-e $job_status);
    }

    return $$self{'Yes'};
}


sub is_finished {
    my ($self, $lane_path, $action) = @_;

    if ($action->{name} eq 'update_db') {
        return $self->{No};
    }

    return $self->SUPER::is_finished($lane_path, $action);
}


#---------- Debugging and error reporting -----------------

sub format_msg
{
    my ($self,@msg) = @_;
    return '['. scalar gmtime() ."]\t". join('',@msg);
}

sub warn
{
    my ($self,@msg) = @_;
    my $msg = $self->format_msg(@msg);
    if ($self->verbose > 0) 
    {
        print STDERR $msg;
    }
    $self->log($msg);
}

sub debug
{
    # The granularity of verbose messaging does not make much sense
    #   now, because verbose cannot be bigger than 1 (made Base.pm
    #   throw on warn's).
    my ($self,@msg) = @_;
    if ($self->verbose > 0) 
    {
        my $msg = $self->format_msg(@msg);
        print STDERR $msg;
        $self->log($msg);
    }
}

sub throw
{
    my ($self,@msg) = @_;
    my $msg = $self->format_msg(@msg);
    Utils::error($msg);
}

sub log
{
    my ($self,@msg) = @_;

    my $msg = $self->format_msg(@msg);
    my $status  = open(my $fh,'>>',$self->log_file);
    if ( !$status ) 
    {
        print STDERR $msg;
    }
    else 
    { 
        print $fh $msg; 
    }
    if ( $fh ) { close($fh); }
}


1;

