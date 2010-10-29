package VertRes::Pipelines::SNPs;
use base qw(VertRes::Pipeline);

use strict;
use warnings;
use LSF;
use Utils;
use VertRes::Parser::sam;
use VertRes::Utils::FileSystem;
use VertRes::Wrapper::GATK;

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
);

our $options = 
{
    bcf_fix         => 'bcf-fix.pl',
    bcftools        => 'bcftools',
    bsub_opts_varfilter  => "-q normal -M3000000 -R 'select[type==X86_64 && mem>3000] rusage[mem=3000,thouio=1]'",
    bsub_opts_mpileup    => "-q normal -M3000000 -R 'select[type==X86_64 && mem>3000] rusage[mem=3000,thouio=1]'",
    bsub_opts_qcall      => "-q normal -M3000000 -R 'select[type==X86_64 && mem>3000] rusage[mem=3000,thouio=1]'",
    bsub_opts_gatk       => "-q normal -M3000000 -R 'select[type==X86_64 && mem>3000] rusage[mem=3000,thouio=1]'",
    fai_chr_regex   => '\d+|x|y',
    gatk_opts       => { all=>{verbose=>1} },
    dbSNP_rod       => undef,
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
    samtools_pileup_params => '-d 500',   # mouse: '-r 0.0001 -d 500'
    vcf_rmdup       => 'vcf-rmdup',
    vcf_stats       => 'vcf-stats',
    vcfutils        => 'vcfutils.pl',
};


# --------- OO stuff --------------

=head2 new

        Options    : See Pipeline.pm for general options and the code for default values.

                    bcf_fix         .. Script to fix invalid VCF
                    bcftools        .. The bcftools binary for mpileup task
                    dbSNP_rod       .. The dbSNP file for GATK
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
                    qcall_cmd       .. The qcall command.
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
                    vcfutils        .. vcfutils.pl from bcftools package (filter4vcf)

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
            $pos += $split_size - 250;
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

    my $out = qq[
use strict;
use warnings;
use Utils;
];

    $out .= qq[Utils::CMD("rm -f $name.vcf-tmp.gz.part");\n];
    $out .= qq[Utils::CMD("rm -f $name.columns");\n];
    for my $chunk (@$chunks)
    {
        # There are too many files, the argument list is often too long for the shell to accept.
        $out .= qq[Utils::CMD("zcat ${name}_$chunk.vcf.gz | gzip -c >> $name.vcf-tmp.gz.part");\n];

        # Check that the columns are ordered in the same order: collect the header line with column names
        #   The zcat followed by head will return SIGPIPE status if more than 200 lines are present. 
        $out .= qq[Utils::CMD("zcat ${name}_$chunk.vcf.gz | head -200 | grep ^#CHROM >> $name.columns",{ignore_errno=>36096});\n];
    }

    $out .= qq[
# Check that the columns are ordered in the same order
my \@out=Utils::CMD(qq[cat $name.columns | uniq | wc -l]);
if ( scalar \@out!=1 or !(\$out[0]=~/^1\$/) ) { Utils::error("FIXME, the column names do not agree: $name.columns\\n"); }

# Take the VCF header from one file and sort the rest. Empty files are allowed by ignoring the exit status 256 from grep.
Utils::CMD(qq[(zcat ${name}_$$chunks[0].vcf.gz | grep ^#; zcat $name.vcf-tmp.gz.part | grep -v ^# | $$self{sort_cmd} -k1,1 -k2,2n) | $$self{vcf_rmdup} | bgzip -c > $name.vcf.gz.part],{ignore_errno=>256});
Utils::CMD(qq[zcat $name.vcf.gz.part | $$self{vcf_stats} > $name.vcf.gz.stats]);
rename('$name.vcf.gz.part','$name.vcf.gz') or Utils::error("rename $name.vcf.gz.part $name.vcf.gz: \$!");
Utils::CMD(qq[tabix -f -p vcf $name.vcf.gz]);
unlink('$name.vcf-tmp.gz.part');
unlink('$name.columns');
];

    for my $chunk (@$chunks)
    {
        $out .= qq[unlink('${name}_$chunk.vcf.gz');\n]; 
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
    # Make sure xargs does not split the arguments to multiple commands by sort -c
    Utils::CMD("cat $name.chunks.list | xargs sort -k1,1 -k2,2n -m | sort -c -k1,1 -k2,2n | gzip -c > $name.pileup.gz.part");
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
    if ( !$$self{bcf_fix} ) { $self->throw("Missing the option bcf_fix.\n"); }

    my $bams = [ $$self{file_list} ];
    my %opts = 
    (
        dir        => $dir,
        work_dir   => "$dir/mpileup",
        job_type   => 'mpileup',
        bsub_opts  => {bsub_opts=>$$self{bsub_opts_mpileup},dont_wait=>1,append=>0},
        bams       => $bams,
        split_size => $$self{split_size_mpileup},

        split_chunks    => \&mpileup_split_chunks,    
        merge_chunks    => \&glue_vcf_chunks,
        merge_vcf_files => \&mpileup_postprocess,
    );

    $self->run_in_parallel(\%opts);

    return $$self{Yes};
}

sub run_mpileup
{
    my ($self,$file_list,$name,$chunk) = @_;

    Utils::CMD(qq[cat $file_list | xargs $$self{mpileup_cmd} -r $chunk -f $$self{fa_ref} | $$self{bcftools} view -gcv - | $$self{bcf_fix} | bgzip -c > $name.vcf.gz.part],{verbose=>1});
    rename("$name.vcf.gz.part","$name.vcf.gz") or $self->throw("rename $name.vcf.gz.part $name.vcf.gz: $!");

    Utils::CMD("touch _$name.done",{verbose=>1});
}


sub mpileup_split_chunks
{
    my ($self,$bam,$chunk_name,$chunk) = @_;

    my $opts = $self->dump_opts(qw(file_list fa_ref fai_ref mpileup_cmd bcftools bcf_fix vcfutils));

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

rename("$basename.vcf.gz","$name.unfilt.vcf.gz") or Utils::error("rename $basename.vcf.gz $name.unfilt.vcf.gz: \$!");
rename("$basename.vcf.gz.tbi","$name.unfilt.vcf.gz.tbi") or Utils::error("rename $basename.vcf.gz.tbi $name.unfilt.vcf.gz.tbi: \$!");
rename("$basename.vcf.gz.stats","$name.unfilt.vcf.gz.stats") or Utils::error("rename $basename.vcf.gz.stats $name.unfilt.vcf.gz.stats: \$!");

Utils::CMD("zcat $name.unfilt.vcf.gz | $$self{vcfutils} filter4vcf | bgzip -c > $basename.filt.vcf.gz");
Utils::CMD("zcat $basename.filt.vcf.gz | $$self{vcf_stats} > $name.vcf.gz.stats");
Utils::CMD("tabix -f -p vcf $basename.filt.vcf.gz");
rename("$basename.filt.vcf.gz","$name.vcf.gz") or Utils::error("rename $basename.filt.vcf.gz $name.vcf.gz: \$!");
rename("$basename.filt.vcf.gz.tbi","$name.vcf.gz.tbi") or Utils::error("rename $basename.filt.vcf.gz.tbi $name.vcf.gz.tbi: \$!");

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
    if ( !$$self{dbSNP_rod} ) { $self->throw("Missing the option dbSNP_rod.\n"); }
    if ( !$$self{gatk_opts} ) { $self->throw("Missing the option gatk_opts.\n"); }

    if ( !exists($$self{gatk_opts}{all}) ) { $$self{gatk_opts}{all}={}; }
    my $gopts = $$self{gatk_opts}{all};
    if ( exists($$gopts{dbsnp}) && $$gopts{dbsnp} ne $$self{dbSNP_rod} )
    {
        $self->throw("Conflicting options for dbSNP, which one to use: $$gopts{dbsnp} or $$self{dbSNP_rod}?\n");
    }
    if ( !exists($$gopts{dbsnp}) ) { $$gopts{dbsnp} = $$self{dbSNP_rod}; }

    if ( exists($$gopts{reference}) && $$gopts{reference} ne $$self{fa_ref} )
    {
        $self->throw("Conflicting options for reference, which one to use: $$gopts{reference} or $$self{fa_ref}?\n");
        $$gopts{reference} = $$self{fa_ref};
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

    my $opts = $self->dump_opts(qw(file_list fa_ref fai_ref sam2vcf gatk_opts dbSNP_rod indel_mask));

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
    $gatk->unified_genotyper($bam,"$name.vcf.gz",L=>$chunk);
    Utils::CMD("tabix -f -p vcf $name.vcf.gz");

    if ( !$$self{indel_mask} or !-e $$self{indel_mask} )
    {
        %opts = $self->get_gatk_opts(qw(all indel_genotyper));
        $gatk = VertRes::Wrapper::GATK->new(%opts);
        $gatk->indel_genotyper($bam, "$name.indels.raw.bed", "$name.indels.detailed.bed",L=>$chunk);
        Utils::CMD("make_indel_mask.pl $name.indels.raw.bed 10 $name.indels.mask.bed",{verbose=>1});

        %opts = $self->get_gatk_opts(qw(all variant_filtration));
        $gatk = VertRes::Wrapper::GATK->new(%opts);
        $gatk->set_b("variant,VCF,$name.vcf.gz", "mask,Bed,$name.indels.mask.bed");
        $gatk->variant_filtration("$name.filtered.vcf.gz");

        # unlink("$name.vcf.gz") or $self->throw("unlink $name.vcf.gz: $!");
        # unlink("$name.indels.detailed.bed");
        # unlink("$name.indels.mask.bed");
        # unlink("$name.indels.mask.bed.idx");
        # unlink("$name.indels.raw.bed");
        warn(qq[delete: $name.vcf.gz $name.indels.detailed.bed $name.indels.mask.bed $name.indels.mask.bed.idx $name.indels.raw.bed]);

        rename("$name.filtered.vcf.gz","$name.vcf.gz") or $self->throw("rename $name.filtered.vcf.gz $name.vcf.gz: $!");
    }

    Utils::CMD("touch _$name.done",{verbose=>1});
}


sub gatk_postprocess
{
    my ($self,$vcfs,$name) = @_;

    if ( scalar @$vcfs > 1 ) { $self->throw("FIXME: expected one file only\n"); }

    my $opts = $self->dump_opts(qw(file_list fa_ref fai_ref sam2vcf gatk_opts dbSNP_rod indel_mask));

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
        $gatk->variant_filtration("$basename.filtered.vcf.gz");
        rename("$basename.filtered.vcf.gz","$basename.tmp.vcf.gz") or $self->throw("rename $basename.filtered.vcf.gz $basename.tmp.vcf.gz: $!");
        Utils::CMD("tabix -f -p vcf $basename.tmp.vcf.gz");
    }
    else
    {
        Utils::CMD("cp $basename.vcf.gz $basename.tmp.vcf.gz");
        Utils::CMD("tabix -f -p vcf $basename.tmp.vcf.gz");
    }

    # Recalibrate
    %opts = $self->get_gatk_opts(qw(all generate_variant_clusters));
    $gatk = VertRes::Wrapper::GATK->new(%opts);
    $gatk->set_b("input,VCF,$basename.tmp.vcf.gz");
    $gatk->set_annotations('HaplotypeScore', 'SB', 'QD', 'HRun');
    $gatk->generate_variant_clusters("$basename.filtered.clusters");

    %opts = $self->get_gatk_opts(qw(all variant_recalibrator));
    $gatk = VertRes::Wrapper::GATK->new(%opts);
    $gatk->set_b("input,VCF,$basename.tmp.vcf.gz");
    $gatk->variant_recalibrator("$basename.filtered.clusters","$basename.recalibrated.vcf.gz");
    rename("$basename.recalibrated.vcf.gz","$basename.tmp.vcf.gz") or $self->throw("rename $basename.recalibrated.vcf.gz $basename.tmp.vcf.gz: $!");

    # Filter the calls down to an implied 10.0% novel false discovery rate
    %opts = $self->get_gatk_opts(qw(all apply_variant_cuts));
    $gatk = VertRes::Wrapper::GATK->new(%opts);
    $gatk->set_b("input,VCF,$basename.tmp.vcf.gz");
    Utils::CMD("tabix -f -p vcf $basename.tmp.vcf.gz");
    $gatk->apply_variant_cuts("$basename.recalibrated.vcf.gz.dat.tranches","$basename.cutted.vcf.gz");
    rename("$basename.cutted.vcf.gz","$basename.tmp.vcf.gz") or $self->throw("rename $basename.cutted.vcf.gz $basename.tmp.vcf.gz: $!");
    Utils::CMD("tabix -f -p vcf $basename.tmp.vcf.gz");

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

    # Execute mpileup with QCall
    Utils::CMD(qq[cat $file_list | xargs samtools mpileup -g -r $chunk -f $$self{fa_ref} | $$self{bcftools} view -Q - | $$self{qcall_cmd} -sn $$self{prefix}$chunk.names -co $chunk.vcf.part],{verbose=>1});

    # Compress the VCF file
    Utils::CMD("cat $chunk.vcf.part | gzip -c > $chunk.vcf.gz.part",{verbose=>1});
    unlink("$chunk.vcf.part");

    # Clean
    unlink("$$self{prefix}$chunk.names");
    rename("$chunk.vcf.gz.part","$chunk_name.vcf.gz") or $self->throw("rename $chunk.vcf.gz.part $chunk_name.vcf.gz: $!");
    Utils::CMD("touch _$chunk_name.done");
}



#---------- qcall_ori ---------------------


# Requires the bam files listed in the file
sub qcall_ori_requires
{
    my ($self,$dir) = @_;
    return [$$self{file_list}];
}


sub qcall_ori_provides
{
    my ($self,$dir) = @_;
    my @provides = ('qcall_ori.done');
    return \@provides;
}

sub qcall_ori
{
    my ($self,$dir,$lock_file) = @_;

    if ( !$$self{fa_ref} ) { $self->throw("Missing the option fa_ref.\n"); }
    if ( !$$self{fai_ref} ) { $self->throw("Missing the option fai_ref.\n"); }
    if ( !$$self{qcall_cmd} ) { $self->throw("Missing the option qcall_cmd.\n"); }
    if ( !$$self{tmp_dir} ) { $self->throw("Missing the option tmp_dir.\n"); }

    my ($count) = `cat $$self{file_list} | wc -l`;
    chomp($count);
    if ( $count<3 ) { $self->throw("QCall is population-based algorithm, requires at least three samples (BAM files) to run.\n"); }

    my $bams = [ $$self{file_list} ];
    my %opts =
    (
        dir        => $dir,
        work_dir   => "$dir/qcall_ori",
        job_type   => 'qcall_ori',
        bsub_opts  => {bsub_opts=>$$self{bsub_opts_qcall},dont_wait=>1,append=>0},
        bams       => $bams,
        split_size => $$self{split_size_qcall},
        tmp_dir    => $$self{tmp_dir},

        split_chunks    => \&qcall_ori_split_chunks,
        merge_chunks    => \&glue_vcf_chunks,
        merge_vcf_files => undef,
    );

    $self->run_in_parallel(\%opts);

    return $$self{Yes};
}

sub qcall_ori_split_chunks
{
    my ($self,$bam,$chunk_name,$chunk) = @_;

    my $opts = $self->dump_opts(qw(file_list fa_ref fai_ref split_size samtools_pileup_params qcall_cmd sort_cmd prefix tmp_dir));

    return qq[
use strict;
use warnings;
use VertRes::Pipelines::SNPs;

my $opts

my \$var = VertRes::Pipelines::SNPs->new(%\$opts);
\$var->run_qcall_ori_chunk(q[$bam],q[$chunk_name],q[$chunk]);
    ];
}


# Run pileups in a more sensible way which takes *MUCH MORE SPACE* on lustre
#   but does not fill tmps on nodes.
#
#   samtools view -b A-RG1.bam chr:from-to > tmp_dir/A-RG1.bam
#   samtools view -b B-RG1.bam chr:from-to > tmp_dir/B-RG1.bam
#   samtools merge tmp_dir/RG1.bam tmp_dir/A-RG1.bam tmp_dir/B-RG1.bam 
#   samtools index tmp_dir/RG1.bam
#   samtools pileup tmp_dir/RG1.bam | samtools glfview - > tmp_dir/RG1.glf.txt
#   rm -f tmp_dir/A-RG1.bam tmp_dir/B-RG1.bam tmp_dir/RG1.bam
#   ...
#   samtools view -b A-RGn.bam chr:from-to > tmp_dir/A-RGn.bam
#   samtools view -b B-RGn.bam chr:from-to > tmp_dir/B-RGn.bam
#   samtools merge tmp_dir/RGn.bam tmp_dir/A-RGn.bam tmp_dir/B-RGn.bam 
#   samtools index tmp_dir/RGn.bam
#   samtools pileup tmp_dir/RGn.bam | samtools glfview - > tmp_dir/RGn.glf.txt
#   rm -f tmp_dir/A-RGn.bam tmp_dir/B-RGn.bam tmp_dir/RGn.bam
#
#   ls tmp_dir/ | grep .glf.txt$ | xargs sort -m | QCall | ...
#   rm -f tmp_dir/*.glf.txt
#
sub run_qcall_ori_chunk
{
    my ($self,$file_list,$chunk_name,$chunk) = @_;

    my $files  = $self->read_files($$self{file_list});
    my $groups = $self->bam_file_groups($files);    # Cluster the BAM files into groups by RG
    my $prefix = $self->common_prefix($files);      # The common prefix will be stripped off the names
    my %names;  # Unique IDs, column names for the output VCF file

    my $tmp_dir = $$self{tmp_dir} .'/'. $chunk;
    Utils::CMD("mkdir -p $tmp_dir",{verbose=>1});

    # Prepare the QCall command line arguments
    for my $grp (keys %$groups)
    {
        my $id;
        if ( $grp eq $$groups{$grp}[0] ) 
        {
            # If we are here, the file name is considered to be the group name
            my ($dir,$name,$suffix) = Utils::basename($grp);
            $dir =~ s{^/*$prefix/*}{};
            $id  = $dir ? "$dir/$name" : $name;
        }
        else
        {
            $id = $grp;
            $id =~ s/\s/_/g;    # Not sure if this can happen, just in case
        }

        if ( exists($names{$id}) ) { $self->throw("FIXME: the names not unique [$grp] -> [$id]\n"); }
        $names{$id} = 1;

        # How many files are in the group? If there is only one, the merge step can be skipped.
        my $nfiles = @{$$groups{$grp}};
        my $bam_to_use = $$groups{$grp}[0];
        if ( $nfiles>1 )
        {
            my $to_be_merged;
            for (my $i=0; $i<@{$$groups{$grp}}; $i++)
            {
                my $file = $$groups{$grp}[$i];
                Utils::CMD(qq[samtools view -b $file $chunk > $tmp_dir/$id.$i.bam],{verbose=>1});
                $to_be_merged .= " $tmp_dir/$id.$i.bam";
            }
            Utils::CMD(qq[samtools merge $tmp_dir/$id.bam $to_be_merged],{verbose=>1});
            Utils::CMD(qq[rm -f $to_be_merged],{verbose=>1});
            Utils::CMD(qq[samtools index $tmp_dir/$id.bam],{verbose=>1});
            $bam_to_use = "$tmp_dir/$id.bam";
        }

        Utils::CMD(qq[samtools pileup $$self{samtools_pileup_params} -gs -f $$self{fa_ref} $bam_to_use | samtools glfview - | awk '{printf("%s\\t$id\\n",\$0);}' > $tmp_dir/$id.glf.txt],{verbose=>1});
        Utils::CMD(qq[rm -f $tmp_dir/$id.bam*]);
    }

    # Write the column names for QCall
    open(my $fh,'>',"$$self{prefix}$chunk.names") or $self->throw("$$self{prefix}$chunk.names: $!");
    for my $id (keys %names)
    {
        print $fh "$id\n";
    }
    close($fh);

    # Execute QCall
    Utils::CMD(qq[/bin/ls $tmp_dir | grep .glf.txt\$ | xargs -I '{}' sort -k1,1n -k2,2n -m '$tmp_dir/{}' | $$self{qcall_cmd} -sn $$self{prefix}$chunk.names -co $chunk.vcf.part],{verbose=>1});
    Utils::CMD(qq[rm -rf $tmp_dir/],{verbose=>1});

    # Compress the VCF file
    Utils::CMD("cat $chunk.vcf.part | gzip -c > $chunk.vcf.gz.part",{verbose=>1});
    unlink("$chunk.vcf.part");

    # Clean
    rename("$chunk.vcf.gz.part","$chunk_name.vcf.gz") or $self->throw("rename $chunk.vcf.gz.part $chunk_name.vcf.gz: $!");
    Utils::CMD("touch _$chunk_name.done");
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

