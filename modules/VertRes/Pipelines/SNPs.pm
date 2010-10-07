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
    bsub_opts       => "-q normal -M6000000 -R 'select[type==X86_64 && mem>6000] rusage[mem=6000,thouio=1]'",
    bsub_opts_long  => "-q normal -M7000000 -R 'select[type==X86_64 && mem>7000] rusage[mem=7000,thouio=1]'",
    fai_chr_regex   => '\d+|x|y',
    gatk_opts       => { verbose=>1, java_memory=>2800, U=>'ALLOW_UNSET_BAM_SORT_ORDER' },
    dbSNP_rod       => undef,
    max_jobs        => undef,
    merge_vcf       => 'merge-vcf -d',
    qcall_cmd       => 'QCALL -ct 0.01 -snpcan', # mouse: -pphet 0',
    sort_cmd        => 'sort',
    sam2vcf         => 'sam2vcf.pl',
    split_size      => 1_000_000,
    gatk_split_size => 1_000_000,
    varfilter       => 'samtools.pl varFilter -S 20 -i 20',
    pileup_rmdup    => 'pileup-rmdup',
    samtools_pileup_params => '-d 500',   # mouse: '-r 0.0001 -d 500'
    vcf_rmdup       => 'vcf-rmdup',
    vcf_stats       => 'vcf-stats',
};


# --------- OO stuff --------------

=head2 new

        Options    : See Pipeline.pm for general options and the code for default values.

                    dbSNP_rod       .. The dbSNP file for GATK
                    file_list       .. File name containing a list of bam files (e.g. the 17 mouse strains). 
                    fa_ref          .. The reference sequence in fasta format
                    fai_ref         .. The reference fai file to read the chromosomes and lengths.
                    fai_chr_regex   .. The chromosomes to be processed.
                    gatk_opts       .. Options to pass to the GATK wrapper
                    gatk_split_size .. GATK seems to require more memory
                    indel_mask      .. GATK precomputed indel mask for filtering SNPs [optional]
                    max_jobs        .. The maximum number of running jobs per task
                    merge_vcf       .. The merge-vcf script.
                    pileup_rmdup    .. The script to remove duplicate positions.
                    qcall_cmd       .. The qcall command.
                    sam2vcf         .. The convertor from samtools pileup format to VCF.
                    samtools_pileup_params .. The options to samtools.pl varFilter (Used by Qcall and varFilter.)
                    sort_cmd        .. Change e.g. to 'sort -T /big/space'
                    split_size      .. The size of the chunks (default is 1Mb).
                    varfilter       .. The samtools varFilter command (samtools.pl varFilter).
                    vcf_rmdup       .. The script to remove duplicate positions.

=cut

sub new 
{
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(%$options,'actions'=>\@actions,@args);
    $self->write_logs(1);

    if ( !$$self{file_list} ) { $self->throw("Missing the option file_list.\n"); }
    $$self{fsu} = VertRes::Utils::FileSystem->new();

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
        if ( !-e $line ) { $self->throw("The file does not exist: [$line]\n"); }
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

        if ( @samples != 1 ) 
        { 
            warn("Multiple samples present, coalescing into one named \"$file\"\n");
            $groups{$file} = $file; 
        }
        else { push @{$groups{$samples[0]}}, $file; }
    }
    return \%groups;
}


sub read_chr_lengths
{
    my ($self,$fai) = @_;

    # Determine the chromosomes and their lengths
    open(my $fh,'<',$fai) or $self->throw("$fai: $!"); 
    my %chr_lengths;
    while (my $line=<$fh>)
    {
        if ( !($line=~/^($$self{fai_chr_regex})\t(\d+)/i) ) { next; }
        $chr_lengths{$1} = $2;
    }
    close($fh);
    return \%chr_lengths;
}


# Create a list of chromosomes splitted into slightly overlapping chunks in the format chr:start-end.
sub chr_chunks
{
    my ($self,$fai_ref,$split_size) = @_;

    # It takes too long to call pileup on the big bam files - split them into chunks. Here, create
    #   the list of chunks.
    my @chunks;
    my $chr_lengths = $self->read_chr_lengths($fai_ref);
    while (my ($chr,$len)=each %$chr_lengths)
    {
        my $pos=1;
        while ($pos<$len)
        {
            my $from = $pos;
            my $to   = $from+$split_size-1;

            # GATK will fail if the region goes beyond the end of the chromosome
            if ( $to>$len ) { $to=$len; }

            push @chunks, "$chr:$from-$to";

            # Create some overlap between the chunks, duplicate records will be removed
            #   later. (There would be overlap SNPs anyway, at least for samtools pileup,
            #   and it would be hard to decide which of the SNPs should be kept. This way,
            #   all covering reads should be there and the SNPs should have similar/identical 
            #   depth.)
            $pos += $split_size - 250;
            if ( $pos<1 ) { $self->throw("The split size too small [$split_size]?\n"); }
        
            # Eearly exit for debugging: do two chunks for one chromosome only for testing.
            #   if ( scalar @chunks>1 ) { return \@chunks; }
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
#   work_dir .. the working directory, should be of the form dir/(varFilter|qcall|gatk)
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
    my $chunks = $self->chr_chunks($$self{fai_ref},$split_size);

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
        # There are too many files, the argument list is often too long for the shell to accept
        $out .= qq[Utils::CMD("zcat ${name}_$chunk.vcf.gz | gzip -c >> $name.vcf-tmp.gz.part");\n];

        # Check that the columns are ordered in the same order: collect the header line with column names
        #   The head will return non-zero status if less than 200 lines are present.
        $out .= qq[Utils::CMD("zcat ${name}_$chunk.vcf.gz | head -200 | grep ^#CHROM >> $name.columns",{exit_on_error=>0});\n];
    }

    $out .= qq[
# Check that the columns are ordered in the same order
my \@out=Utils::CMD(qq[cat $name.columns | uniq | wc -l]);
if ( scalar \@out!=1 or !(\$out[0]=~/^1\$/) ) { Utils::error("FIXME, the column names do not agree: $name.columns\\n"); }

# Take the VCF header from one file and sort the rest
Utils::CMD(qq[(zcat ${name}_$$chunks[0].vcf.gz | grep ^#; zcat $name.vcf-tmp.gz.part | grep -v ^# | $$self{sort_cmd} -k1,1 -k2,2n) | $$self{vcf_rmdup} | bgzip -c > $name.vcf.gz.part]);
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
        dir       => $dir,
        work_dir  => "$dir/varFilter",
        job_type  => 'varFilter',
        bsub_opts => {bsub_opts=>$$self{bsub_opts_long},dont_wait=>1,append=>0},
        bams      => $bams,

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
    if ( !-e "$name.pileup.gz" )
    {
        Utils::CMD(qq[samtools view -bh $bam $chunk | samtools pileup $$self{samtools_pileup_params} -c -f $$self{fa_ref} - | $$self{varfilter} | gzip -c > $name.pileup.gz.part],{verbose=>1});
        rename("$name.pileup.gz.part","$name.pileup.gz") or $self->throw("rename $name.pileup.gz.part $name.pileup.gz: $!");
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
    my ($self,$chunks,$name) = @_;

    my $concat  = qq[Utils::CMD("rm -f $name.pileup-tmp.gz.part");\n];
    my $rmfiles = '';
    for my $chunk (@$chunks)
    {
        my $chunk_name = $name.'_'.$chunk;

        # There are too many files, the argument list is often too long for the shell to accept
        $concat  .= qq[Utils::CMD("zcat $chunk_name.pileup.gz | gzip -c >> $name.pileup-tmp.gz.part");\n];
        $rmfiles .= qq[Utils::CMD("rm -f $chunk_name.pileup.gz");\n];
    }

    return qq[
use strict;
use warnings;
use Utils;
use VertRes::Utils::FileSystem;

if ( !-e "$name.pileup.gz" )
{
    $concat
    Utils::CMD("zcat $name.pileup-tmp.gz.part | sort -k1,1 -k2,2n | gzip -c > $name.pileup.gz.part");
    rename("$name.pileup.gz.part","$name.pileup.gz") or Utils::error("rename $name.pileup.gz.part $name.pileup.gz: \$!");
    $rmfiles
    Utils::CMD("rm -f $name.pileup-tmp.gz.part");
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

    my $bams = [ $$self{file_list} ];
    my %opts =
    (
        dir        => $dir,
        work_dir   => "$dir/gatk",
        job_type   => 'gatk',
        bsub_opts  => {bsub_opts=>$$self{bsub_opts_long},dont_wait=>1,append=>0},
        split_size => $$self{gatk_split_size},
        bams       => $bams,

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


sub run_gatk_chunk
{
    my ($self,$bam,$name,$chunk) = @_;

    my %opts = ( %{$$self{gatk_opts}} );
    if ( $$self{dbSNP_rod} ) 
    { 
        $opts{dbsnp} = $$self{dbSNP_rod}; 
    }

    my $gatk;

    $gatk = VertRes::Wrapper::GATK->new(%opts);
    $gatk->unified_genotyper($bam,"$name.vcf.gz",L=>$chunk);
    Utils::CMD("tabix -f -p vcf $name.vcf.gz");
    # $self->log("UnifiedGenotyper done: $name.vcf.gz\n");

    if ( !$$self{indel_mask} or !-e $$self{indel_mask} )
    {
        $gatk = VertRes::Wrapper::GATK->new(%{$$self{gatk_opts}});
        $gatk->indel_genotyper($bam, "$name.indels.raw.bed", "$name.indels.detailed.bed",L=>$chunk);
        Utils::CMD("make_indel_mask.pl $name.indels.raw.bed 10 $name.indels.mask.bed",{verbose=>1});
        #$self->log("IndelGenotyperV2 done: $name.indels.raw.bed\n");

        $gatk = VertRes::Wrapper::GATK->new(%{$$self{gatk_opts}});
        $gatk->set_b("variant,VCF,$name.vcf.gz", "mask,Bed,$name.indels.mask.bed");
        $gatk->variant_filtration("$name.filtered.vcf.gz");
        Utils::CMD("tabix -f -p vcf $name.filtered.vcf.gz");

        unlink("$name.vcf.gz") or $self->throw("unlink $name.vcf.gz: $!");
        rename("$name.filtered.vcf.gz","$name.vcf.gz") or $self->throw("rename $name.filtered.vcf.gz $name.vcf.gz: $!");
        #$self->log("renamed: $name.filtered.vcf.gz -> $name.vcf.gz\n");
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

    my $gatk;

    # Filter
    if ( $$self{indel_mask} && -e $$self{indel_mask} )
    {
        $gatk = VertRes::Wrapper::GATK->new(%{$$self{gatk_opts}});
        $gatk->set_b("variant,VCF,$vcf", "mask,Bed,$$self{indel_mask}");
        $gatk->variant_filtration("$basename.filtered.vcf.gz");
        rename("$basename.filtered.vcf.gz","$basename.tmp.vcf.gz") or $self->throw("rename $basename.filtered.vcf.gz $basename.tmp.vcf.gz: $!");
        Utils::CMD("tabix -f -p vcf $basename.tmp.vcf.gz");
        #$self->log("VariantFiltration: $basename.filtered.vcf.gz -> $basename.tmp.vcf.gz\n");
    }
    else
    {
        Utils::CMD("cp $basename.vcf.gz $basename.tmp.vcf.gz");
        Utils::CMD("tabix -f -p vcf $basename.tmp.vcf.gz");
    }

    # Recalibrate
    $gatk = VertRes::Wrapper::GATK->new(%{$$self{gatk_opts}});
    $gatk->set_b("input,VCF,$basename.tmp.vcf.gz");
    $gatk->set_annotations('HaplotypeScore', 'SB', 'QD', 'HRun');
    $gatk->generate_variant_clusters("$basename.filtered.clusters");
    #$self->log("GenerateVariantClusters done\n");

    $gatk = VertRes::Wrapper::GATK->new(%{$$self{gatk_opts}});
    $gatk->set_b("input,VCF,$basename.tmp.vcf.gz");
    $gatk->variant_recalibrator("$basename.filtered.clusters","$basename.recalibrated.vcf.gz");
    rename("$basename.recalibrated.vcf.gz","$basename.tmp.vcf.gz") or $self->throw("rename $basename.recalibrated.vcf.gz $basename.tmp.vcf.gz: $!");
    Utils::CMD("tabix -f -p vcf $basename.tmp.vcf.gz");
    #$self->log("VariantRecalibrator done: $basename.recalibrated.vcf.gz -> $basename.tmp.vcf.gz\n");

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

    my ($count) = `cat $$self{file_list} | wc -l`;
    chomp($count);
    if ( $count<3 ) { $self->throw("QCall is population-based algorithm, requires at least three samples (BAM files) to run.\n"); }

    my $bams = [ $$self{file_list} ];
    my %opts =
    (
         dir       => $dir,
         work_dir  => "$dir/qcall",
         job_type  => 'qcall',
         bsub_opts => {bsub_opts=>$$self{bsub_opts_long},dont_wait=>1,append=>0},
         bams      => $bams,

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

    my $opts = $self->dump_opts(qw(file_list fa_ref fai_ref split_size samtools_pileup_params qcall_cmd sort_cmd prefix));

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

    my $files  = $self->read_files($$self{file_list});
    my $groups = $self->bam_file_groups($files);
    my $prefix = $self->common_prefix($files);
    my %names;

    # Prepare the QCall command line arguments
    my $cmd = "(\n";
    for my $grp (keys %$groups)
    {
        my $id;
        if ( $grp eq $$groups{$grp}[0] ) 
        {
            # The group name is the file name
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

        my $grp_cmd = "\t(\n";
        for my $file (@{$$groups{$grp}})
        {
            $grp_cmd .= qq[samtools view $file $chunk; ];
        }
        $grp_cmd .= "\t) ";
        if ( @{$$groups{$grp}}>1 ) { $grp_cmd .= qq[ | $$self{sort_cmd} -k3,3n -k4,4n ]; }
        $grp_cmd .= qq[ | samtools pileup $$self{samtools_pileup_params} -gsS -f $$self{fa_ref} - | samtools glfview - | awk '{printf("%s\\t$id\\n",\$0);}';\n];
        $cmd .= $grp_cmd;
    }

    # Write the column names for QCall
    open(my $fh,'>',"$$self{prefix}$chunk.names") or $self->throw("$$self{prefix}$chunk.names: $!");
    for my $id (keys %names)
    {
        print $fh "$id\n";
    }
    close($fh);

    # Execute QCall
    $cmd .= ") | $$self{sort_cmd} -k1,1n -k2,2n | $$self{qcall_cmd} -sn $$self{prefix}$chunk.names -co $chunk.vcf.part";
    Utils::CMD($cmd,{verbose=>1});

    # used to go via qcall-to-vcf, but it doesn't like 1/1 genotypes, so skip
    # it for now
    Utils::CMD("cat $chunk.vcf.part | gzip -c > $chunk.vcf.gz.part",{verbose=>1});
    unlink("$chunk.vcf.part");

    rename("$chunk.vcf.gz.part","$chunk_name.vcf.gz") or $self->throw("rename $chunk.vcf.gz.part $chunk_name.vcf.gz: $!");
    Utils::CMD("touch _$chunk_name.done");
}



#---------- Debugging and error reporting -----------------

sub warn
{
    my ($self,@msg) = @_;
    my $msg = join('',@msg);
    if ($self->verbose > 0) 
    {
        print STDERR $msg;
    }
    $self->log($msg);
}

sub debug
{
    my ($self,@msg) = @_;
    if ($self->verbose > 0) 
    {
        my $msg = join('',@msg);
        print STDERR $msg;
        $self->log($msg);
    }
}

sub throw
{
    my ($self,@msg) = @_;
    Utils::error(@msg);
}

sub log
{
    my ($self,@msg) = @_;

    my $msg_str = join('',@msg);
    my $status  = open(my $fh,'>>',$self->log_file);
    if ( !$status ) 
    {
        print STDERR $msg_str;
    }
    else 
    { 
        print $fh $msg_str; 
    }
    if ( $fh ) { close($fh); }
}


1;

