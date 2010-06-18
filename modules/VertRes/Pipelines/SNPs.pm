package VertRes::Pipelines::SNPs;
use base qw(VertRes::Pipeline);

use strict;
use warnings;
use LSF;
use Utils;
use VertRes::Parser::sam;

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
    gatk_cmd        => 'java -Xmx6500m -jar /nfs/users/nfs_p/pd3/sandbox/call-snps/gatk/GenomeAnalysisTK/GenomeAnalysisTK.jar -T UnifiedGenotyper -hets 0.0001 -confidence 30 -mmq 25 -mc 1000 -mrl 10000000 --platform Solexa -U ALLOW_UNSET_BAM_SORT_ORDER',
    max_jobs        => undef,
    merge_vcf       => 'merge-vcf -d',
    qcall_cmd       => 'QCALL -ct 0.01 -snpcan -pphet 0',
    sort_cmd        => 'sort',
    sam2vcf         => 'sam2vcf.pl',
    split_size      => 1000000,
    gatk_split_size => 10_000_000,
    varfilter       => 'samtools.pl varFilter -S 20 -i 20',
    pileup_rmdup    => 'pileup-rmdup',
    samtools_pileup_params => '-r 0.0001 -d 500',
    vcf_rmdup       => 'vcf-rmdup',
    vcf_stats       => 'vcf-stats',
};


# --------- OO stuff --------------

=head2 new

        Options    : See Pipeline.pm for general options and the code for default values.

                    file_list       .. File name containing a list of bam files (e.g. the 17 mouse strains). 
                    fa_ref          .. The reference sequence in fasta format
                    fai_ref         .. The reference fai file to read the chromosomes and lengths.
                    fai_chr_regex   .. The chromosomes to be processed.
                    gatk_cmd        .. The command to launch GATK SNP caller
                    gatk_split_size .. GATK seems to require more memory
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
    if ( -e "$dir/_$name.o" ) { rename("$dir/_$name.o","$outdir/_$name.o"); }
    if ( -e "$dir/_$name.e" ) { rename("$dir/_$name.e","$outdir/_$name.e"); }
    if ( -e "$dir/_$name.jid" ) { unlink("$dir/_$name.jid"); }
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
Utils::CMD(qq[tabix -p vcf $name.vcf.gz]);
    ];
}


# Runs varFilter, gatk and QCall jobs in parallel. All have the same structure: 
#   1) For each bam file, get a list of chromosomes and split them into chunks ($$opts{split_chunks}).
#   2) Then merge the parts into one VCF file and delete the intermediate chunks ($$opts{merge_chunks}).
#   3) Finally, create one big VCF file with all the data, each column corresponds to one bam file ($$opts{merge_vcf_files}).
#
# The parameters:
#   work_dir .. the working directory, should be of the form dir/(varFilter|qcall|gatk)
#
# split_chunks function:
#   parameters .. $bam,$chunk_name,$chunk
#   must produce the file _$chunk_name.done, otherwise it will be run over and over again.
#
# merge_chunks function:
#   parameters .. $chunks,$name
#   must produce the file $name.vcf.gz out of the files from the $dir
#
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

    my $bams = $self->read_files($$self{file_list});
    my $prefix = $self->common_prefix($bams);
    my %names;
    my @to_be_merged;

    # Step 1, call the split_chunks command.
    for my $bam (@$bams)
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
        if ( -e "$work_dir/$outdir/$name.vcf.gz" ) 
        {
            $self->clean_files("$work_dir/$outdir",$name,"$work_dir/$outdir/_done"); 
            for my $chunk (@$chunks)
            {
                my $chunk_name = $name.'_'.$chunk;
                unlink("$work_dir/$outdir/_$chunk_name.done") unless !-e "$work_dir/$outdir/_$chunk_name.done";
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
            if ( -e "$work_dir/$outdir/_$chunk_name.done" || -e "$work_dir/$outdir/$name.vcf.gz" ) 
            { 
                $self->clean_files("$work_dir/$outdir",$chunk_name,"$work_dir/$outdir/_done");
                next; 
            }
            $chunks_finished = 0;
            
            Utils::CMD("mkdir -p $work_dir/$outdir") unless -e "$work_dir/$outdir";

            my $jids_file = "$work_dir/$outdir/_$chunk_name.jid";
            my $status = LSF::is_job_running($jids_file);
            if ( $status&$LSF::Running ) { next; }
            if ( $status&$LSF::Error ) { $self->warn("Some jobs failed: $jids_file\n"); }

            # Now get the splitting command
            my $cmd = &{$$opts{split_chunks}}($self,$bam,$chunk_name,$chunk);

            open(my $fh,'>',"$work_dir/$outdir/_$chunk_name.pl") or $self->throw("$work_dir/$outdir/_$chunk_name.pl: $!");
            print $fh $cmd;
            close($fh);
            LSF::run($jids_file,"$work_dir/$outdir","_$chunk_name",$$opts{bsub_opts},qq[perl -w _$chunk_name.pl]);
            $self->debug("Submitting $work_dir/$outdir .. _$chunk_name\n");
        }

        # Some of the chunks is not ready
        if ( !$chunks_finished ) { next; }

        # Now merge the chunks for one bam file into one VCF file.
        my $jids_file = "$work_dir/$outdir/_$name.jid";
        my $status = LSF::is_job_running($jids_file);
        if ( $status&$LSF::Running ) { next; }
        if ( $status&$LSF::Error ) { $self->warn("Some jobs failed: $jids_file\n"); }

        # Get the merging command
        my $cmd = &{$$opts{merge_chunks}}($self,$chunks,$name);

        open(my $fh,'>',"$work_dir/$outdir/_$name.pl") or $self->throw("$work_dir/$outdir/_$name.pl: $!");
        print $fh $cmd;
        close($fh);
        LSF::run($jids_file,"$work_dir/$outdir","_$name",$$opts{bsub_opts},qq[perl -w _$name.pl]);
        $self->debug("Submitting $work_dir/$outdir .. _$name\n");
    }

    # Some chunks still not done
    if ( !$is_finished ) { return; }

    # The big VCF file from the step 3. already exists.
    if ( -e "$work_dir/$$opts{job_type}.vcf.gz" ) 
    { 
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

    my %opts = 
    (
        dir       => $dir,
        work_dir  => "$dir/varFilter",
        job_type  => 'varFilter',
        bsub_opts => {bsub_opts=>$$self{bsub_opts_long},dont_wait=>1,append=>0},

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
    if ( ! -e "$name.pileup.gz" )
    {
        Utils::CMD(qq[samtools view -bh $bam $chunk | samtools pileup $$self{samtools_pileup_params} -c -f $$self{fa_ref} - | $$self{varfilter} | gzip -c > $name.pileup.gz.part],{verbose=>1});
        rename("$name.pileup.gz.part","$name.pileup.gz") or $self->throw("rename $name.pileup.gz.part $name.pileup.gz: $!");
        Utils::CMD("touch _$name.done",{verbose=>1});
    }
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

    my $files;
    for my $chunk (@$chunks)
    {
        my $chunk_name = $name.'_'.$chunk;
        $files .= " $chunk_name.pileup.gz";
    }

    return qq[
use strict;
use warnings;
use Utils;

if ( ! -e "$name.pileup.gz" )
{
    Utils::CMD("zcat $files | sort -k1,1 -k2,2n | gzip -c > $name.pileup.gz.part");
    rename("$name.pileup.gz.part","$name.pileup.gz") or Utils::error("rename $name.pileup.gz.part $name.pileup.gz: \$!");
    Utils::CMD("rm -f $files");
}
if ( ! -e "$name.vcf.gz" )
{
    Utils::CMD("zcat $name.pileup.gz | $$self{pileup_rmdup} | $$self{sam2vcf} -s -t $name | bgzip -c > $name.vcf.gz.part",{verbose=>1});
    Utils::CMD(qq[zcat $name.vcf.gz.part | $$self{vcf_stats} > $name.vcf.gz.stats]);
    rename("$name.vcf.gz.part","$name.vcf.gz") or Utils::error("rename $name.vcf.gz.part $name.vcf.gz: \$!");
    Utils::CMD(qq[tabix -p vcf $name.vcf.gz]);
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

    if ( !$$self{fai_ref} ) { $self->throw("Missing the option fai_ref.\n"); }
    if ( !$$self{gatk_cmd} ) { $self->throw("Missing the option gatk_cmd.\n"); }

    my %opts =
    (
        dir        => $dir,
        work_dir   => "$dir/gatk",
        job_type   => 'gatk',
        bsub_opts  => {bsub_opts=>$$self{bsub_opts_long},dont_wait=>1,append=>0},
        split_size => $$self{gatk_split_size},

        split_chunks    => \&gatk_split_chunks,
        merge_chunks    => \&gatk_merge_chunks,
        merge_vcf_files => \&merge_vcf_files,
    );

    $self->run_in_parallel(\%opts);

    return $$self{Yes};
}


sub gatk_split_chunks
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
    gatk_cmd   => q[$$self{gatk_cmd}],
};
my \$snps = VertRes::Pipelines::SNPs->new(%\$opts);
\$snps->run_gatk(q[$bam],q[$chunk_name],q[$chunk]);
    ];
}


sub gatk_merge_chunks
{
    my ($self,$chunks,$name) = @_;

    my $header;
    my $files;
    for my $chunk (@$chunks)
    {
        my $chunk_name = $name.'_'.$chunk;
        $files .= " $chunk_name.vcf.gz";
        if ( !$header ) { $header = "$chunk_name.vcf.gz"; }
    }

    return qq[
use strict;
use warnings;
use Utils;

if ( ! -e "$name.vcf.gz" )
{
    # Take the VCF header from one file and sort the rest
    Utils::CMD("(zcat $header | grep ^#; zcat $files | grep -v ^# | sort -k1,1 -k2,2n) | $$self{vcf_rmdup} -d DoC | bgzip -c > $name.vcf.gz.part");
    Utils::CMD(qq[zcat $name.vcf.gz.part | $$self{vcf_stats} > $name.vcf.gz.stats]);
    rename("$name.vcf.gz.part","$name.vcf.gz") or Utils::error("rename $name.vcf.gz.part $name.vcf.gz: \$!");
    Utils::CMD("tabix -p vcf $name.vcf.gz");
    Utils::CMD("rm -f $files");
}
    ];
}

sub run_gatk
{
    my ($self,$bam,$name,$chunk) = @_;
    if ( ! -e "$name.vcf.gz" )
    {
        Utils::CMD("$$self{gatk_cmd} -R $$self{fa_ref} -I $bam -L $chunk | gzip -c > $name.vcf.gz.part",{verbose=>1});
        rename("$name.vcf.gz.part","$name.vcf.gz") or $self->throw("rename $name.vcf.gz.part $name.vcf.gz: $!");
    }
    if ( -e "$name.vcf.gz" )
    {
        Utils::CMD("touch _$name.done",{verbose=>1});
    }
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

    my $files = $self->read_files($$self{file_list});
    if ( @$files<3 ) { $self->throw("QCall is population-based algorithm, requires at least three samples (BAM files) to run.\n"); }

    my $work_dir = "$dir/qcall";
    Utils::CMD("mkdir -p $work_dir") unless -e $work_dir;

    # The big VCF file already exists
    if ( -e "$work_dir/qcall.vcf.gz" )
    {
        Utils::CMD("touch $dir/qcall.done");
        return;
    }

    my $chunks = $self->chr_chunks($$self{fai_ref},$$self{split_size});

    my $max_jobs = $$self{max_jobs} ? $$self{max_jobs} : 0;
    my $done = $LSF::Done;
    my @to_be_merged;
    my $ntasks = 0;

    # For all chromosomes create the 1Mb chunks and schedule them to LSF queue
    for my $chunk (@$chunks)
    {
        # Add to the list of VCF files to be merged later
        push @to_be_merged, "$work_dir/$chunk.vcf.gz";

        if ( -e "$work_dir/$chunk.vcf.gz" )
        {
            $self->clean_files($work_dir,$chunk,"$work_dir/_done");
            next;
        }

        my $jids_file = "$work_dir/$$self{prefix}$chunk.jid";
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
            $self->warn("The command failed: $work_dir .. perl -w $$self{prefix}$chunk.pl\n"); 
        }

        open(my $fh,'>',"$work_dir/$$self{prefix}$chunk.pl") or $self->throw("$work_dir/$$self{prefix}$chunk.pl: $!");
        print $fh qq[
use strict;
use warnings;
use VertRes::Pipelines::SNPs;

my \$opts = {
    file_list  => q[$$self{file_list}],
    fa_ref     => q[$$self{fa_ref}],
    fai_ref    => q[$$self{fai_ref}],
    split_size => q[$$self{split_size}],
    samtools_pileup_params => q[$$self{samtools_pileup_params}],
    qcall_cmd  => q[$$self{qcall_cmd}],
    sort_cmd   => q[$$self{sort_cmd}],
    prefix     => q[$$self{prefix}],
};
my \$var = VertRes::Pipelines::SNPs->new(%\$opts);
\$var->run_qcall_chunk(q[$chunk]);
        \n];

        close($fh);
        LSF::run($jids_file,$work_dir,"$$self{prefix}$chunk",{%$self,bsub_opts=>$$self{bsub_opts_long},dont_wait=>1,append=>0},qq[perl -w $$self{prefix}$chunk.pl]);
        $self->debug("Submitting $work_dir .. perl -w $$self{prefix}$chunk.pl\n");

        $done |= $LSF::Running;
    }

    if ( $done&$LSF::Error ) { $self->throw("QCall failed for some of the files.\n"); }
    if ( $done&$LSF::Running ) { $self->debug("Some files not finished...\n"); return $$self{No}; }

    # Because this subroutine returns as if it has already finished, a custom jids_file must
    #   be used: Pipeline.pm will delete the $lock_file.
    my $jids_file = "$work_dir/$$self{prefix}qcall_merge.jid";
    my $status = LSF::is_job_running($jids_file);
    if ( $status&$LSF::Running ) { return; }
    if ( $status&$LSF::Error ) { $self->warn("Some jobs failed: $jids_file\n"); }

    # Get the VCF merging command
    my $cmd = $self->glue_vcf_chunks(\@to_be_merged,'qcall');

    open(my $fh,'>',"$work_dir/$$self{prefix}merge.pl") or $self->throw("$work_dir/$$self{prefix}merge.pl: $!");
    print $fh $cmd;
    close($fh);
    LSF::run($jids_file,$work_dir,"$$self{prefix}merge",{%$self,bsub_opts=>$$self{bsub_opts_long},dont_wait=>1},qq[perl -w $$self{prefix}merge.pl]);
    $self->debug("Submitting $work_dir $$self{prefix}merge\n");

    return $$self{Yes};
}

# Multiple BAM files can belong to the same sample and should be merged.
sub bam_file_groups
{
    my ($self,$files) = @_;
    my %groups;
    for my $file (@$files)
    {
        my $fh = Utils::CMD("samtools view -H $file",{rpipe=>1});
        my $pars = VertRes::Parser::sam->new(fh=>$fh);
        my @samples = $pars->samples();
        close($fh);

        unless ( @samples == 1 ) { $groups{$file} = $file; }
        else { push @{$groups{$samples[0]}}, $file; }
    }
    return \%groups;
}

sub run_qcall_chunk
{
    my ($self,$chunk) = @_;

    if ( -e "$chunk.vcf.gz" ) { return; }

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

    # Before it's fixed, convert to proper VCF
    Utils::CMD("cat $chunk.vcf.part | qcall-to-vcf | gzip -c > $chunk.vcf.gz.part",{verbose=>1});
    unlink("$chunk.vcf.part");

    rename("$chunk.vcf.gz.part","$chunk.vcf.gz") or $self->throw("rename $chunk.vcf.gz.part $chunk.vcf.gz: $!");
}

sub glue_vcf_chunks
{
    my ($self,$vcfs,$name) = @_;

    my $args = join(' ',@$vcfs);

    my $out = qq[
use strict;
use warnings;
use Utils;
# Take the VCF header from one file and sort the rest
Utils::CMD(qq[(zcat $$vcfs[0] | grep ^#; zcat $args | grep -v ^# | $$self{sort_cmd} -k1,1 -k2,2n) | $$self{vcf_rmdup} | bgzip -c > $name.vcf.gz.part]);
Utils::CMD(qq[zcat $name.vcf.gz.part | $$self{vcf_stats} > $name.vcf.gz.stats]);
rename('$name.vcf.gz.part','$name.vcf.gz') or Utils::error("rename $name.vcf.gz.part $name.vcf.gz: \$!");
Utils::CMD(qq[tabix -p vcf $name.vcf.gz]);

];

    for my $file (@$vcfs)
    {
        if ( !($file=~/.vcf.gz$/) ) { $self->throw("Could not parse the chunk file name: $file"); }
        my $chunk = $`;
        $out .= qq[unlink('$file');\n]; 
        $out .= qq[unlink('$$self{prefix}$chunk.names');\n]; 
    }

    return $out;
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

