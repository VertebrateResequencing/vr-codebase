package VertRes::TrackQC;
use base qw(VertRes::Pipeline);

use strict;
use warnings;
use LSF;
use HierarchyUtilities;
use VertRes::GTypeCheck;

our @actions =
(
    {
        'name'     => 'rename',
        'action'   => \&rename,
        'requires' => \&rename_requires, 
        'provides' => \&rename_provides,
    },
    {
        'name'     => 'subsample',
        'action'   => \&subsample,
        'requires' => \&subsample_requires, 
        'provides' => \&subsample_provides,
    },
    {
        'name'     => 'aln_fastqs',
        'action'   => \&aln_fastqs,
        'requires' => \&aln_fastqs_requires, 
        'provides' => \&aln_fastqs_provides,
    },
    {
        'name'     => 'map_sample',
        'action'   => \&map_sample,
        'requires' => \&map_sample_requires, 
        'provides' => \&map_sample_provides,
    },
    {
        'name'     => 'check_genotype',
        'action'   => \&check_genotype,
        'requires' => \&check_genotype_requires, 
        'provides' => \&check_genotype_provides,
    },
    {
        'name'     => 'create_graphs',
        'action'   => \&create_graphs,
        'requires' => \&create_graphs_requires, 
        'provides' => \&create_graphs_provides,
    },
);


our $options = 
{
    # Executables
    'bwa_exec'       => 'bwa-0.4.9',
    'gcdepth_R'      => '/nfs/users/nfs_p/pd3/cvs/maqtools/mapdepth/gcdepth/gcdepth.R',
    'glf'            => '/nfs/users/nfs_p/pd3/cvs/glftools/glfv3/glf',
    'mapviewdepth'   => 'mapviewdepth_sam',
    'samtools'       => 'samtools',

    'bsub_opts'      => "-q normal -M5000000 -R 'select[type==X86_64 && mem>5000] rusage[mem=5000]'",
    'gc_depth_bin'   => 20000,
    'min_glf_ratio'  => 5.0,
    'sample_dir'     => 'qc-sample',
    'sample_size'    => 50e6,
};


# --------- OO stuff --------------

=head2 new

        Example    : my $qc = VertRes::TrackDummy->new( 'sample_dir'=>'dir', 'sample_size'=>1e6 );
        Options    : See Pipeline.pm for general options.

                    # Executables
                    bwa_exec        .. bwa executable
                    gcdepth_R       .. gcdepth R script
                    glf             .. glf executable
                    mapviewdepth    .. mapviewdepth executable
                    samtools        .. samtools executable

                    # Options specific to TrackQC
                    bsub_opts       .. LSF bsub options for jobs
                    bwa_ref         .. overrides the value returned by HierarchyUtilities::lane_info
                    fa_ref          .. see HierarchyUtilities::lane_info
                    fai_ref         .. see HierarchyUtilities::lane_info
                    gc_depth_bin    .. the bin size for the gc-depth graph
                    min_glf_ratio   .. the minimum distinctive glf likelihood ratio
                    snps            .. see HierarchyUtilities::lane_info
                    sample_dir      .. where to put subsamples
                    sample_size     .. the size of the subsample

=cut

sub new 
{
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(%$options,'actions'=>\@actions,@args);
    $self->write_logs(1);
    return $self;
}


=head2 clean

        Description : If mrProper option is set, the entire QC directory will be deleted.
        Returntype  : None

=cut

sub clean
{
    my ($self) = @_;

    $self->SUPER::clean();

    if ( !$$self{'lane'} ) { $self->throw("Missing parameter: the lane to be cleaned.\n"); }
    if ( !$$self{'sample_dir'} ) { $self->throw("Missing parameter: the sample_dir to be cleaned.\n"); }

    if ( !$$self{'mrProper'} ) { return; }
    my $qc_dir = qq[$$self{'lane'}/$$self{'sample_dir'}];
    if ( ! -d $qc_dir ) { return; }

    $self->debug("rm -rf $qc_dir\n");
    Utils::CMD(qq[rm -rf $qc_dir]);
    return;
}

#---------- rename ---------------------

# Requires nothing
sub rename_requires
{
    my ($self,$lane) = @_;
    my @requires = ();
    return \@requires;
}

# It may provide also _2.fastq.gz, but we assume that if
#   there is _1.fastq.gz but _2 is missing, it is OK.
#
sub rename_provides
{
    my ($self,$lane) = @_;
    my $info = HierarchyUtilities::lane_info($lane);
    my @provides = ("$$info{'lane'}_1.fastq.gz");
    return \@provides;
}

# The naming convention is to name fastq files according
#   to the lane, e.g. 
#       project/sample/tech/libr/lane/lane_1.fastq.gz
#       project/sample/tech/libr/lane/lane_2.fastq.gz
#
sub rename
{
    my ($self,$lane_path,$lock_file) = @_;

    my $lane_info = HierarchyUtilities::lane_info($lane_path);
    my $name      = $$lane_info{'lane'};

    my $fastq_files = existing_fastq_files("$lane_path/$name");

    # They are named correctly.
    if ( scalar @$fastq_files  ) { return $$self{'Yes'}; }

    my @files = glob("$lane_path/*.fastq.gz");
    if ( scalar @files > 2 ) { Utils::error("FIXME: so far can handle up to 2 fastq files: $lane_path.\n") }

    my $i = 0;
    for my $file (sort cmp_last_number @files)
    {
        $i++;
        Utils::relative_symlink($file,"$lane_path/${name}_$i.fastq.gz");
    }
    return $$self{'Yes'};
}

sub cmp_last_number($$)
{
    my ($a,$b) = @_;
    if ( !($a=~/(\d+)\D*$/) ) { return 0 } # leave unsorted if there is no number in the name
    my $x = $1;
    if ( !($b=~/(\d+)\D*$/) ) { return 0 }
    my $y = $1;
    return $x<=>$y;
}

# How many fastq files there are of the name ${prefix}_[123].fastq.gz?
sub existing_fastq_files
{
    my ($prefix) = @_;

    my @files = ();
    my $i = 1;
    while ( -e "${prefix}_$i.fastq.gz" )
    {
        push @files, "${prefix}_$i.fastq.gz";
        $i++;
    }
    return \@files;
}


#---------- subsample ---------------------

# At least one fastq file is required. We assume that either all fastq files
#   are in place, or there is none.
#
sub subsample_requires
{
    my ($self,$lane) = @_;

    my $info = HierarchyUtilities::lane_info($lane);
    my @requires = ("$$info{'lane'}_1.fastq.gz");
    return \@requires;
}

# We must sample all fastq files but we do not know how many there are.
sub subsample_provides
{
    my ($self,$lane) = @_;

    my $sample_dir = $$self{'sample_dir'};
    my $info = HierarchyUtilities::lane_info($lane);
    my $name = $$info{'lane'};

    my $fastq_files = existing_fastq_files("$lane/$name");
    my $nfiles = scalar @$fastq_files;
    if ( !$nfiles ) { Utils::error("No fastq.gz files in $lane??\n") }

    my @provides = ();
    for (my $i=1; $i<=$nfiles; $i++)
    {
        push @provides, "$sample_dir/${name}_$i.fastq.gz";
    }
    return \@provides;
}

sub subsample
{
    my ($self,$lane_path,$lock_file) = @_;

    my $sample_dir = $$self{'sample_dir'};
    my $lane_info  = HierarchyUtilities::lane_info($lane_path);
    my $name       = $$lane_info{'lane'};

    Utils::create_dir("$lane_path/$sample_dir/");

    # Dynamic script to be run by LSF.
    open(my $fh, '>', "$lane_path/$sample_dir/_qc-sample.pl") or Utils::error("$lane_path/$sample_dir/_qc-sample.pl: $!");

    my $fastq_files = existing_fastq_files("$lane_path/$name");
    my $nfiles = scalar @$fastq_files;
    if ( !$nfiles ) { Utils::error("No fastq files in $lane_path??") }

    # The files will be created in reverse order, so that that _1.fastq.gz is created 
    #   last - the next action checks only for the first one. If there is only one file,
    #   the variable $seq_list is not used. If there are multiple, $seq_list is passed
    #   only to subsequent calls of FastQ::sample.
    #
    print $fh "use FastQ;\n";
    print $fh '$seq_list = ' unless $nfiles==1;
    for (my $i=$nfiles; $i>0; $i--)
    {
        my $print_seq_list = ($i==$nfiles) ? '' : ', $seq_list';
        print $fh "FastQ::sample(q{../${name}_$i.fastq.gz},q{${name}_$i.fastq.gz}, $$self{'sample_size'}$print_seq_list);\n";
    }
    close $fh;

    LSF::run($lock_file,"$lane_path/$sample_dir","_${name}_sample", $self, qq{perl -w _qc-sample.pl});

    return $$self{'No'};
}



#----------- aln_fastqs ---------------------

# If one sample is in place (_1.fastq.gz), we assume all are in place.
sub aln_fastqs_requires
{
    my ($self,$lane) = @_;

    my $sample_dir = $$self{'sample_dir'};
    my $info = HierarchyUtilities::lane_info($lane);

    my @requires = ("$sample_dir/$$info{'lane'}_1.fastq.gz");
    return \@requires;
}

sub aln_fastqs_provides
{
    my ($self,$lane) = @_;

    my @provides = ();
    my $sample_dir = $$self{'sample_dir'};
    my $info = HierarchyUtilities::lane_info($lane);
    my $name = $$info{'lane'};

    my $fastq_files = existing_fastq_files("$lane/$sample_dir/$name");
    my $nfiles = scalar @$fastq_files;
    if ( !$nfiles ) 
    {
        @provides = ("$sample_dir/${name}_1.sai");
        return \@provides;
    }
    
    for (my $i=1; $i<=$nfiles; $i++)
    {
        push @provides, "$sample_dir/${name}_$i.sai";
    }
    return \@provides;
}

sub aln_fastqs
{
    my ($self,$lane_path,$lock_file) = @_;

    my $lane_info = HierarchyUtilities::lane_info($lane_path);
    my $name      = $$lane_info{'lane'};
    my $work_dir  = "$lane_path/$$self{sample_dir}";
    my $prefix    = exists($$self{'prefix'}) ? $$self{'prefix'} : '_';

    my $bwa       = $$self{'bwa_exec'};
    my $bwa_ref   = exists($$self{'bwa_ref'}) ? $$self{'bwa_ref'} : $$lane_info{'bwa_ref'};
    my $fai_ref   = exists($$self{'fai_ref'}) ? $$self{'fai_ref'} : $$lane_info{'fai_ref'};
    my $samtools  = $$self{'samtools'};

    # How many files do we have?
    my $fastq_files = existing_fastq_files("$work_dir/$name");
    my $nfiles = scalar @$fastq_files;
    if ( $nfiles<1 || $nfiles>2 ) { Utils::error("FIXME: we can handle 1 or 2 fastq files in $work_dir, not $nfiles.\n") }

    for (my $i=1; $i<=$nfiles; $i++)
    {
        open(my $fh,'>', "$work_dir/${prefix}aln_fastq_$i.pl") or Utils::error("$work_dir/${prefix}aln_fastq_$i.pl: $!");
        print $fh
qq[
use Utils;

Utils::CMD("$bwa aln -l 32 $bwa_ref ${name}_$i.fastq.gz > ${name}_$i.saix");
if ( ! -s "${name}_$i.saix" ) { Utils::error("The command ended with an error:\n\t$bwa aln -l 32 $bwa_ref ${name}_$i.fastq.gz > ${name}_$i.saix\n") }
rename("${name}_$i.saix","${name}_$i.sai") or Utils::CMD("rename ${name}_$i.saix ${name}_$i.sai: \$!");

];
        LSF::run($lock_file,$work_dir,"_${name}_$i",$self,qq[perl -w ${prefix}aln_fastq_$i.pl]);
    }

    return $$self{'No'};
}



#----------- map_sample ---------------------


sub map_sample_requires
{
    my ($self,$lane) = @_;

    my @requires = ();
    my $sample_dir = $$self{'sample_dir'};
    my $info = HierarchyUtilities::lane_info($lane);
    my $name = $$info{'lane'};

    my $fastq_files = existing_fastq_files("$lane/$sample_dir/$name");
    my $nfiles = scalar @$fastq_files;
    if ( !$nfiles ) 
    {
        @requires = ("$sample_dir/${name}_1.sai");
        return \@requires;
    }
    
    for (my $i=1; $i<=$nfiles; $i++)
    {
        push @requires, "$sample_dir/${name}_$i.sai";
    }
    return \@requires;
}

sub map_sample_provides
{
    my ($self,$lane) = @_;
    my $sample_dir = $$self{'sample_dir'};
    my $info = HierarchyUtilities::lane_info($lane);
    my @provides = ("$sample_dir/$$info{'lane'}.bam");
    return \@provides;
}

sub map_sample
{
    my ($self,$lane_path,$lock_file) = @_;

    my $sample_dir = $$self{'sample_dir'};
    my $lane_info  = HierarchyUtilities::lane_info($lane_path);
    my $name       = $$lane_info{'lane'};
    my $work_dir   = "$lane_path/$$self{sample_dir}";

    my $bwa        = $$self{'bwa_exec'};
    my $bwa_ref    = exists($$self{'bwa_ref'}) ? $$self{'bwa_ref'} : $$lane_info{'bwa_ref'};
    my $fai_ref    = exists($$self{'fai_ref'}) ? $$self{'fai_ref'} : $$lane_info{'fai_ref'};
    my $samtools   = $$self{'samtools'};


    # How many files do we have?
    my $fastq_files = existing_fastq_files("$work_dir/$name");
    my $nfiles = scalar @$fastq_files;
    if ( $nfiles<1 || $nfiles>2 ) { Utils::error("FIXME: we can handle 1 or 2 fastq files in $work_dir, not $nfiles.\n") }


    my $bwa_cmd  = '';
    if ( $nfiles == 1 )
    {
        $bwa_cmd = "$bwa samse $bwa_ref ${name}_1.sai ${name}_1.fastq.gz";
    }
    else
    {
        $bwa_cmd = "$bwa sampe $bwa_ref ${name}_1.sai ${name}_2.sai ${name}_1.fastq.gz ${name}_2.fastq.gz";
    }

    # Dynamic script to be run by LSF. We must check that the bwa exists alright
    #   - samtools do not return proper status and create .bam file even if there was
    #   nothing read from the input.
    #
    open(my $fh,'>', "$work_dir/_map.pl") or Utils::error("$work_dir/_map.pl: $!");
    print $fh 
qq{
use Utils;

Utils::CMD("$bwa_cmd > ${name}.sam");
if ( ! -s "${name}.sam" ) { Utils::error("The command ended with an error:\n\t$bwa_cmd > ${name}.sam\n") }

Utils::CMD("$samtools import $fai_ref ${name}.sam ${name}.ubam");
Utils::CMD("$samtools sort ${name}.ubam ${name}x");     # Test - will this help from NFS problems?
Utils::CMD("rm -f ${name}.sam ${name}.ubam");
rename("${name}x.bam", "$name.bam") or Utils::CMD("rename ${name}x.bam $name.bam: \$!");
};

    LSF::run($lock_file,$work_dir,"_${name}_sampe",$self, q{perl -w _map.pl});
    return $$self{'No'};
}




#----------- check_genotype ---------------------

sub check_genotype_requires
{
    my ($self,$lane) = @_;
    my $sample_dir = $$self{'sample_dir'};
    my $info = HierarchyUtilities::lane_info($lane);
    my @requires = ("$sample_dir/$$info{'lane'}.bam");
    return \@requires;
}

sub check_genotype_provides
{
    my ($self,$lane) = @_;
    my $sample_dir = $$self{'sample_dir'};
    my $info = HierarchyUtilities::lane_info($lane);
    my @provides = ("$sample_dir/$$info{'lane'}.gtype");
    return \@provides;
}

sub check_genotype
{
    my ($self,$lane_path,$lock_file) = @_;

    my $lane_info  = HierarchyUtilities::lane_info($lane_path);
    my $name       = $$lane_info{'lane'};

    my $options = {};
    $$options{'bam'}           = "$lane_path/$$self{'sample_dir'}/$name.bam";
    $$options{'bsub_opts'}     = $$self{'bsub_opts'};
    $$options{'fa_ref'}        = exists($$self{'fa_ref'}) ? $$self{'fa_ref'} : $$lane_info{'fa_ref'};
    $$options{'glf'}           = $$self{'glf'};
    $$options{'snps'}          = exists($$self{'snps'}) ? $$self{'snps'} : $$lane_info{'snps'};
    $$options{'samtools'}      = $$self{'samtools'};
    $$options{'genotype'}      = exists($$lane_info{'genotype'}) ? $$lane_info{'genotype'} : '';
    $$options{'min_glf_ratio'} = exists($$lane_info{'gtype_confidence'}) ? $$lane_info{'gtype_confidence'} : $$self{'min_glf_ratio'};
    $$options{'prefix'}        = $$self{'prefix'};
    $$options{'lock_file'}     = $lock_file;

    my $gtc = VertRes::GTypeCheck->new(%$options);
    $gtc->check_genotype();

    return $$self{'No'};
}


#----------- create_graphs ---------------------

sub create_graphs_requires
{
    my ($self,$lane) = @_;
    my $sample_dir = $$self{'sample_dir'};
    my $info = HierarchyUtilities::lane_info($lane);
    my @requires = ("$sample_dir/$$info{'lane'}.bam");
    return \@requires;
}

sub create_graphs_provides
{
    my ($self,$lane) = @_;
    my $sample_dir = $$self{'sample_dir'};
    my $info = HierarchyUtilities::lane_info($lane);
    my @provides = ("$sample_dir/chrom-distrib.png","$sample_dir/gc-content.png","$sample_dir/insert-size.png","$sample_dir/gc-depth.png");
    return \@provides;
}

sub create_graphs
{
    my ($self,$lane_path,$lock_file) = @_;

    my $sample_dir = $$self{'sample_dir'};
    my $lane_info  = HierarchyUtilities::lane_info($lane_path);

    # Dynamic script to be run by LSF.
    open(my $fh, '>', "$lane_path/$sample_dir/_graphs.pl") or Utils::error("$lane_path/$sample_dir/_graphs.pl: $!");
    print $fh 
qq[
use VertRes::TrackQC;

my \$params = 
{
    'gc_depth_bin' => q[$$self{'gc_depth_bin'}],
    'mapviewdepth' => q[$$self{'mapviewdepth'}],
    'samtools'     => q[$$self{'samtools'}],
    'gcdepth_R'    => q[$$self{'gcdepth_R'}],
    'lane_path'    => q[$lane_path],
    'sample_dir'   => q[$$self{'sample_dir'}],
};

VertRes::TrackQC::run_graphs(\$params,\$\$params{lane_path});
];
    close $fh;

    LSF::run($lock_file,"$lane_path/$sample_dir","_$$lane_info{'lane'}_graphs", $self, qq{perl -w _graphs.pl});
    return $$self{'No'};
}


sub run_graphs
{
    my ($self,$lane_path) = @_;

    use Graphs;
    use SamTools;
    use Utils;
    use FastQ;

    my $sample_dir   = $$self{'sample_dir'};
    my $lane_info    = HierarchyUtilities::lane_info($lane_path);
    my $name         = $$lane_info{'lane'};
    my $outdir       = "$lane_path/$sample_dir/";
    my $bam_file     = "$outdir/$name.bam";

    my $samtools     = $$self{'samtools'};
    my $mapview      = $$self{'mapviewdepth'};
    my $refseq       = $$lane_info{'fa_ref'};
    my $gc_depth_bin = $$self{'gc_depth_bin'};
    my $bindepth    = "$outdir/gc-depth.bindepth";
    my $gcdepth_R    = $$self{'gcdepth_R'};

    my $stats_file  = "$outdir/_stats";
    my $other_stats = "$outdir/_detailed-stats.txt";
    my $dump_file   = "$outdir/_stats.dump";

    my $fastq_files = existing_fastq_files("$lane_path/$name");
    for (my $i=1; $i<=scalar @$fastq_files; $i++)
    {
        my $fastqcheck = "$lane_path/${name}_$i.fastq.gz.fastqcheck";
        if ( !-e $fastqcheck ) { next }

        my $data = FastQ::parse_fastqcheck($fastqcheck);
        $$data{'outfile'}    = "$outdir/fastqcheck_$i.png";
        $$data{'title'}      = "FastQ Check $i";
        $$data{'desc_xvals'} = 'Sequencing Quality';
        $$data{'desc_yvals'} = '1000 x Frequency / nBases';

        # Draw the 'Total' line as the last one and somewhat thicker
        my $total = shift(@{$$data{'data'}});
        $$total{'lines'} = ',lwd=3';
        push @{$$data{'data'}}, $total;

        Graphs::plot_stats($data);
    }

    if ( ! -e "$outdir/gc-depth.png" || Utils::file_newer($bam_file,$bindepth) )
    {
        Utils::CMD("$samtools view $bam_file | $mapview $refseq -b=$gc_depth_bin > $bindepth");
        Graphs::create_gc_depth_graph($bindepth,$gcdepth_R,qq[$outdir/gc-depth.png]);
    }

    my $all_stats = SamTools::collect_detailed_bam_stats($bam_file,$$lane_info{'fai_ref'});
    my $stats = $$all_stats{'total'};
    report_stats($stats,$lane_path,$stats_file);
    report_detailed_stats($stats,$lane_path,$other_stats);
    dump_detailed_stats($stats,$dump_file);

    my ($x,$y);
    $x = $$stats{'insert_size'}{'max'}{'x'};
    $y = $$stats{'insert_size'}{'max'}{'y'};
    Graphs::plot_stats({
            'outfile'    => qq[$outdir/insert-size.png],
            'title'      => 'Insert Size',
            'desc_yvals' => 'Frequency',
            'desc_xvals' => 'Insert Size',
            'data'       => [ $$stats{'insert_size'} ],
            'r_cmd'      => qq[text($x,$y,'$x',pos=4,col='darkgreen')\n],
            'r_plot'     => "xlim=c(0," . ($$lane_info{'insert_size'}*2.5) . ")",
            });

    $x = $$stats{'gc_content_forward'}{'max'}{'x'};
    $y = $$stats{'gc_content_forward'}{'max'}{'y'};
    my @gc_data = $$stats{'gc_content_forward'};
    if ( $$stats{'gc_content_reverse'} ) { push @gc_data, $$stats{'gc_content_reverse'}; }
    Graphs::plot_stats({
            'outfile'    => qq[$outdir/gc-content.png],
            'title'      => 'GC Content (both mapped and unmapped)',
            'desc_yvals' => 'Frequency',
            'desc_xvals' => 'GC Content [%]',
            'data'       => \@gc_data,
            'r_cmd'      => "text($x,$y,'$x',pos=4,col='darkgreen')\n",
            });

    Graphs::plot_stats({
            'barplot'    => 1,
            'outfile'    => qq[$outdir/chrom-distrib.png],
            'title'      => 'Chromosome Coverage',
            'desc_yvals' => 'Frequency/Length',
            'desc_xvals' => 'Chromosome',
            'data'       => [ $$stats{'reads_chrm_distrib'}, ],
            });
}


sub report_stats
{
    my ($stats,$lane_path,$outfile) = @_;

    my $info = HierarchyUtilities::lane_info($lane_path);

    my $avg_read_length = $$stats{'bases_total'}/$$stats{'reads_total'};

    open(my $fh,'>',$outfile) or Utils::error("$outfile: $!");
    print  $fh "MAPPED,,$$info{project},$$info{sample},$$info{technology},$$info{library},$$info{lane},";
    printf $fh ",0,$$info{lane}_1.fastq.gz,%.1f,$$info{lane}_2.fastq.gz,%.1f,", $avg_read_length,$avg_read_length;
    print  $fh "$$stats{reads_total},$$stats{bases_total},$$stats{reads_mapped},$$stats{bases_mapped_cigar},$$stats{reads_paired},";
    print  $fh "$$stats{rmdup_reads_total},0,$$stats{error_rate}\n";
    close $fh;
}



sub report_detailed_stats
{
    my ($stats,$lane_path,$outfile) = @_;

    open(my $fh,'>',$outfile) or Utils::error("$outfile: $!");

    printf $fh "reads total .. %d\n", $$stats{'reads_total'};
    printf $fh "     mapped .. %d (%.1f%%)\n", $$stats{'reads_mapped'}, 100*($$stats{'reads_mapped'}/$$stats{'reads_total'});
    printf $fh "     paired .. %d (%.1f%%)\n", $$stats{'reads_paired'}, 100*($$stats{'reads_paired'}/$$stats{'reads_total'});
    printf $fh "bases total .. %d\n", $$stats{'bases_total'};
    printf $fh "    mapped (read)  .. %d (%.1f%%)\n", $$stats{'bases_mapped_read'}, 100*($$stats{'bases_mapped_read'}/$$stats{'bases_total'});
    printf $fh "    mapped (cigar) .. %d (%.1f%%)\n", $$stats{'bases_mapped_cigar'}, 100*($$stats{'bases_mapped_cigar'}/$$stats{'bases_total'});
    printf $fh "duplication .. %f\n", $$stats{'duplication'};
    printf $fh "error rate  .. %f\n", $$stats{error_rate};
    printf $fh "\n";
    printf $fh "insert size        \n";
    printf $fh "    average .. %.1f\n", $$stats{insert_size}{average};
    printf $fh "    std dev .. %.1f\n", $$stats{insert_size}{std_dev};
    printf $fh "\n";
    printf $fh "chrm distrib dev .. %f\n", $$stats{'reads_chrm_distrib'}{'scaled_dev'};

    close $fh;
}


sub dump_detailed_stats
{
    my ($stats,$outfile) = @_;

    use Data::Dumper;
    open(my $fh,'>',$outfile) or Utils::error("$outfile: $!");
    print $fh Dumper($stats);
    close $fh;
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
    if ($self->verbose > 1) 
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
}


1;

