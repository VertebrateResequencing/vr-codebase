=head1 NAME

VertRes::Pipelines::TrackQC_Fastq - pipeline for QC of fastq files, inherits from VertRes::Pipelines::TrackQC_Bam.

=head1 SYNOPSIS

See /lustre/scratch102/conf/pipeline.conf and /lustre/scratch102/conf/qc-g1k.conf
for an example.

=cut

package VertRes::Pipelines::TrackQC_Fastq;
use base qw(VertRes::Pipelines::TrackQC_Bam);

use strict;
use warnings;
use LSF;
use VertRes::Parser::fastqcheck;
use VertRes::Wrapper::bwa;
use VertRes::Utils::Sam;
use VRTrack::Assembly;
use Pathogens::Parser::GenomeCoverage;

our @actions =
(
    # Takes care of naming convention, fastq.gz file names should match
    #   the lane name. 
    {
        'name'     => 'rename_files',
        'action'   => \&rename_files,
        'requires' => \&rename_files_requires, 
        'provides' => \&rename_files_provides,
    },

    # Creates a smaller subsample out of the fastq files to speed up
    #   the QC pipeline. 
    {
        'name'     => 'subsample',
        'action'   => \&subsample,
        'requires' => \&subsample_requires, 
        'provides' => \&subsample_provides,
    },

    # Runs bwa to create the .sai files and checks for the presence of 
    #   adapter sequences.
    {
        'name'     => 'process_fastqs',
        'action'   => \&process_fastqs,
        'requires' => \&process_fastqs_requires, 
        'provides' => \&process_fastqs_provides,
    },

    # Runs bwa to create the bam file.
    {
        'name'     => 'map_sample',
        'action'   => \&map_sample,
        'requires' => \&map_sample_requires, 
        'provides' => \&map_sample_provides,
    },

    # Runs glf to check the genotype. Inherited from TrackQC_Bam.
    {
        'name'     => 'check_genotype',
        'action'   => \&VertRes::Pipelines::TrackQC_Bam::check_genotype,
        'requires' => \&VertRes::Pipelines::TrackQC_Bam::check_genotype_requires, 
        'provides' => \&VertRes::Pipelines::TrackQC_Bam::check_genotype_provides,
    },
    
    # Finds percentage of transposons in reads if applicable
    {
      'name'     => 'transposon',
      'action'   => \&transposon,
      'requires' => \&transposon_requires, 
      'provides' => \&transposon_provides,
    },

    # Creates some QC graphs and generate some statistics.
    {
        'name'     => 'stats_and_graphs',
        'action'   => \&stats_and_graphs,
        'requires' => \&stats_and_graphs_requires, 
        'provides' => \&stats_and_graphs_provides,
    },

    # Checks the generated stats and attempts to auto pass or fail the lane. Inherited from TrackQC_Bam.
    {
        'name'     => 'auto_qc',
        'action'   => \&VertRes::Pipelines::TrackQC_Bam::auto_qc,
        'requires' => \&VertRes::Pipelines::TrackQC_Bam::auto_qc_requires, 
        'provides' => \&VertRes::Pipelines::TrackQC_Bam::auto_qc_provides,
    },

    # Writes the QC status to the tracking database. Inherited from TrackQC_Bam.
    {
        'name'     => 'update_db',
        'action'   => \&VertRes::Pipelines::TrackQC_Bam::update_db,
        'requires' => \&VertRes::Pipelines::TrackQC_Bam::update_db_requires, 
        'provides' => \&VertRes::Pipelines::TrackQC_Bam::update_db_provides,
    },
);

our $options = 
{
    # Executables
    'blat'            => '/software/pubseq/bin/blat',
    'bwa_exec'        => 'bwa-0.5.3',
    'gcdepth_R'       => '/software/vertres/bin-external/gcdepth.R',
    'glf'             => 'glf',
    'mapviewdepth'    => 'mapviewdepth_sam',
    'samtools'        => 'samtools',

    'adapters'        => '/software/pathogen/projects/protocols/ext/solexa-adapters.fasta',
    'bsub_opts'       => "-q normal -M5000 -R 'select[type==X86_64] select[mem>5000] rusage[mem=5000] rusage[thouio=1]'",
    'bsub_opts_stats_and_graphs'   => "-q normal -M1000 -R 'select[type==X86_64] select[mem>1000] rusage[mem=1000] rusage[thouio=1]'",
    'bsub_opts_map_sample'         => "-q normal -M5000 -R 'select[type==X86_64] select[mem>5000] rusage[mem=5000] rusage[thouio=1]'",
    'bsub_opts_process_fastqs'     => "-q normal -M3200 -R 'select[type==X86_64] select[mem>3200] rusage[mem=3200] rusage[thouio=1]'",
    'bsub_opts_subsample'          => "-q normal -M1000 -R 'select[type==X86_64] select[mem>1000] rusage[mem=1000] rusage[thouio=1]'",
    'bwa_clip'        => 20,
    'chr_regex'       => '^(?:\d+|X|Y)$',
    'gc_depth_bin'    => 20000,
    'gtype_confidence'=> 5.0,
    'sample_dir'      => 'qc-sample',
    'sample_size'     => 50e6,
    'stats'           => '_stats',
    'stats_detailed'  => '_detailed-stats.txt',
    'stats_dump'      => '_stats.dump',
};


# --------- OO stuff --------------

=head2 new

        Example    : my $qc = VertRes::Pipelines::TrackQC_Fastq->new( 'sample_dir'=>'dir', 'sample_size'=>1e6 );
        Options    : See Pipeline.pm for general options.

                    # Executables
                    blat            .. blat executable
                    bwa_exec        .. bwa executable
                    gcdepth_R       .. gcdepth R script
                    glf             .. glf executable
                    mapviewdepth    .. mapviewdepth executable
                    samtools        .. samtools executable

                    # Options specific to TrackQC
                    adapters        .. the location of .fa with adapter sequences
                    assembly        .. e.g. NCBI36
                    bsub_opts       .. LSF bsub options for jobs
                    bwa_clip        .. The value to the 'bwa aln -q ' command.
                    bwa_ref         .. the prefix to reference files, as required by bwa
                    chr_regex       .. For chromosome coverage graph (e.g. '^(?:\d+|X|Y)$')
                    fa_ref          .. the reference sequence in fasta format
                    fai_ref         .. the index to fa_ref generated by samtools faidx
                    gc_depth_bin    .. the bin size for the gc-depth graph
                    gtype_confidence.. the minimum expected glf likelihood ratio
                    insert_size     .. the maximum expected insert size (default is 1000)
                    paired          .. is the lane from paired-end sequencing?
                    snps            .. genotype file generated by hapmap2bin from glftools
                    sample_dir      .. where to put subsamples
                    sample_size     .. the size of the subsample (approx number of bases in one fastq file)
                    stats_ref       .. e.g. /path/to/NCBI36.stats

=cut

sub VertRes::Pipelines::TrackQC_Fastq::new 
{
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(%$options,'actions'=>\@actions,@args);
    $self->write_logs(1);

    if ( !$$self{bwa_exec} ) { $self->throw("Missing the option bwa_exec.\n"); }
    if ( !$$self{gcdepth_R} ) { $self->throw("Missing the option gcdepth_R.\n"); }
    if ( !$$self{glf} ) { $self->throw("Missing the option glf.\n"); }
    if ( !$$self{mapviewdepth} ) { $self->throw("Missing the option mapviewdepth.\n"); }
    if ( !$$self{samtools} ) { $self->throw("Missing the option samtools.\n"); }
    if ( !$$self{fa_ref} ) { $self->throw("Missing the option fa_ref.\n"); }
    if ( !$$self{fai_ref} ) { $self->throw("Missing the option fai_ref.\n"); }
    if ( !$$self{gc_depth_bin} ) { $self->throw("Missing the option gc_depth_bin.\n"); }
    if ( !$$self{gtype_confidence} ) { $self->throw("Missing the option gtype_confidence.\n"); }
    if ( !$$self{sample_dir} ) { $self->throw("Missing the option sample_dir.\n"); }
    if ( !$$self{sample_size} ) { $self->throw("Missing the option sample_size.\n"); }

    # Set mapper + version 
    my $mapper = VertRes::Wrapper::bwa->new(exe => $self->{'bwa_exec'});
    $self->{mapper} = 'bwa';
    $self->{mapper_version} = $mapper->version;

    return $self;
}



#---------- rename_files ---------------------

# Requires nothing
sub rename_files_requires
{
    my ($self) = @_;
    my @requires = ();
    return \@requires;
}

# It may provide also _2.fastq.gz and _3.fastq.gz, but we assume that if
#   there is _1.fastq.gz and the others are missing, it is OK.
#
sub rename_files_provides
{
    my ($self) = @_;
    my @provides = ("$$self{lane}_1.fastq.gz");
    return \@provides;
}

# The naming convention is to name fastq files according
#   to the lane, e.g. 
#       project/sample/tech/libr/lane/lane_1.fastq.gz    .. pair 1
#       project/sample/tech/libr/lane/lane_2.fastq.gz    .. pair 2
#       project/sample/tech/libr/lane/lane_3.fastq.gz    .. single
#
# The lanes named like "lane_s_1" complicate the algorithm: for 
#   paired reads, it should be ignored.  For unpaired, it should be used.
#   Luckily, these were treated already in Import.pm. Rename will have a
#   problem only with lanes which did not go through Import.pm and contain
#   these _s_ files. 
#
sub rename_files
{
    my ($self,$lane_path,$lock_file) = @_;

    my $name = $$self{lane};

    my $fastq_files = existing_fastq_files("$lane_path/$name");

    # They are named correctly.
    if ( scalar @$fastq_files  ) { return $$self{'Yes'}; }

    my @files = glob("$lane_path/*.fastq.gz");
    if ( scalar @files > 3 ) { Utils::error("FIXME: so far can handle up to 3 fastq files: $lane_path.\n") }

    # If there are three files, the one without a number or numbered as 3 is non-paired, the 1,2 are paired
    my $i = 1;
    for my $file (sort cmp_last_number @files)
    {
        if ( ($file=~/^\d+_s_\d+/) && scalar @files>1 )
        {
            # If there are more fastq files and some is named as _s_, it is not clear which are
            #   paired, which are non-paired and if they should be used at all.
            $self->throw(qq[Heuristic failed, did not expect "$file". Please symlink manually.]);
        }
        elsif ( $i<3 && !($file=~/_\d+\D*$/) && scalar @files>1 )
        {
            $self->throw(qq[Heuristic failed, which files are single-ended and which are paired? Please symlink manually.]); 
        }
        elsif ( $i==3 && ($file=~/_\d+\D*$/) ) 
        { 
            $self->throw(qq[Heuristic failed, this should be single-ended file "$file". Please symlink manually.]); 
        }

        if ( ! -e "$lane_path/${name}_$i.fastq.gz" )
        {
            Utils::relative_symlink($file,"$lane_path/${name}_$i.fastq.gz");
        }
        if ( -e "$file.fastqcheck" && ! -e "$lane_path/${name}_$i.fastq.gz.fastqcheck" )
        {
            Utils::relative_symlink("$file.fastqcheck","$lane_path/${name}_$i.fastq.gz.fastqcheck");
        }
        if ( -e "$file.md5" && ! -e "$lane_path/${name}_$i.fastq.gz.md5" )
        {
            Utils::relative_symlink("$file.md5","$lane_path/${name}_$i.fastq.gz.md5");
        }
        
        $i++;
    }
    return $$self{'Yes'};
}

sub cmp_last_number($$)
{
    my ($a,$b) = @_;
    if ( !($a=~/_(\d+)\D*$/) ) { return 1 } # is greater than the rest if there is no number in the name
    my $x = $1;
    if ( !($b=~/_(\d+)\D*$/) ) { return -1 }
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
    my ($self) = @_;

    my @requires = ("$$self{lane}_1.fastq.gz");
    return \@requires;
}

# We must sample all fastq files but we do not know how many there are.
sub subsample_provides
{
    my ($self,$lane_path) = @_;

    my $sample_dir = $$self{'sample_dir'};
    my $name = $$self{lane};

    my $fastq_files = existing_fastq_files("$lane_path/$name");
    my $nfiles = scalar @$fastq_files;
    if ( !$nfiles ) { $self->throw("No fastq.gz files in $lane_path??\n") }

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
    my $name       = $$self{'lane'};

    Utils::create_dir("$lane_path/$sample_dir/");

    # Dynamic script to be run by LSF.
    open(my $fh, '>', "$lane_path/$sample_dir/_qc-sample.pl") or $self->throw("$lane_path/$sample_dir/_qc-sample.pl: $!");

    my $fastq_files = existing_fastq_files("$lane_path/$name");
    my $nfiles = scalar @$fastq_files;
    if ( !$nfiles ) { $self->throw("No fastq files in $lane_path??") }

    # The files will be created in reverse order, so that that _1.fastq.gz is created 
    #   last - the next action checks only for the first one. If there is only one file,
    #   the variable $seq_list is not used. If there are multiple, $seq_list is passed
    #   only to subsequent calls of FastQ::sample.
    #
    print $fh "use FastQ;\n";
    for (my $i=$nfiles; $i>0; $i--)
    {
        my $print_seq_list = '';
        if ( $i==2 )
        {
            print $fh '$seq_list = ';
        }
        elsif ( $i==1 )
        {
            $print_seq_list = ', $seq_list';
        }
        print $fh "FastQ::sample(q{../${name}_$i.fastq.gz},q{${name}_$i.fastq.gz}, $$self{'sample_size'}$print_seq_list);\n";
    }
    close $fh;

    LSF::run($lock_file,"$lane_path/$sample_dir","_${name}_sample", {bsub_opts=>$$self{bsub_opts_subsample}}, qq{perl -w _qc-sample.pl});

    return $$self{'No'};
}



#----------- process_fastqs ---------------------

# If one sample is in place (_1.fastq.gz), we assume all are in place.
sub process_fastqs_requires
{
    my ($self) = @_;

    my $sample_dir = $$self{'sample_dir'};
    my @requires = ("$sample_dir/$$self{lane}_1.fastq.gz");
    return \@requires;
}

sub process_fastqs_provides
{
    my ($self,$lane_path) = @_;

    my @provides = ();
    my $sample_dir = $$self{'sample_dir'};
    my $name = $$self{lane};

    my $fastq_files = existing_fastq_files("$lane_path/$sample_dir/$name");
    my $nfiles = scalar @$fastq_files;

    # This should not happen, only when the import is broken.
    if ( !$nfiles ) { return 0; }
    
    for (my $i=1; $i<=$nfiles; $i++)
    {
        push @provides, "$sample_dir/${name}_$i.sai";
        push @provides, "$sample_dir/${name}_$i.nadapters";
    }
    return \@provides;
}

sub process_fastqs
{
    my ($self,$lane_path,$lock_file) = @_;

    if ( !$$self{bwa_ref} ) { $self->throw("Missing the option bwa_ref.\n"); }

    my $name      = $$self{lane};
    my $work_dir  = "$lane_path/$$self{sample_dir}";
    my $prefix    = exists($$self{'prefix'}) ? $$self{'prefix'} : '_';

    my $bwa       = $$self{'bwa_exec'};
    my $bwa_ref   = $$self{'bwa_ref'};
    my $fai_ref   = $$self{'fai_ref'};
    my $samtools  = $$self{'samtools'};

    # How many files do we have?
    my $fastq_files = existing_fastq_files("$work_dir/$name");
    my $nfiles = scalar @$fastq_files;
    if ( $nfiles<1 || $nfiles>3 ) { Utils::error("FIXME: we can handle 1-3 fastq files in $work_dir, not $nfiles.\n") }

    # Run bwa aln for each fastq file to create .sai files.
    for (my $i=1; $i<=$nfiles; $i++)
    {
        if ( -e qq[$work_dir/${name}_$i.sai] ) { next; }

        open(my $fh,'>', "$work_dir/${prefix}aln_fastq_$i.pl") or Utils::error("$work_dir/${prefix}aln_fastq_$i.pl: $!");
        print $fh
qq[
use strict;
use warnings;
use Utils;

Utils::CMD("$bwa aln -q $$self{bwa_clip} -l 32 $bwa_ref ${name}_$i.fastq.gz > ${name}_$i.saix");
if ( ! -s "${name}_$i.saix" ) 
{ 
    Utils::error("The command ended with an error:\n\t$bwa aln -q $$self{bwa_clip} -l 32 $bwa_ref ${name}_$i.fastq.gz > ${name}_$i.saix\n");
}
rename("${name}_$i.saix","${name}_$i.sai") or Utils::error("rename ${name}_$i.saix ${name}_$i.sai: \$!");

];
        close($fh);
        LSF::run($lock_file,$work_dir,"_${name}_$i",{bsub_opts=>$$self{bsub_opts_process_fastqs}},qq[perl -w ${prefix}aln_fastq_$i.pl]);
    }

    # Run blat for each fastq file to find out how many adapter sequences are in there.
    Utils::CMD(qq[cat $$self{adapters} | sed 's/>/>ADAPTER|/' > $work_dir/adapters.fa]);

    for (my $i=1; $i<=$nfiles; $i++)
    {
        if ( -e qq[$work_dir/${name}_$i.nadapters] ) { next; }

        open(my $fh,'>', "$work_dir/${prefix}blat_fastq_$i.pl") or Utils::error("$work_dir/${prefix}blat_fastq_$i.pl: $!");
        print $fh
qq[
use strict;
use warnings;
use Utils;

Utils::CMD(q[zcat ${name}_$i.fastq.gz | awk '{print ">"substr(\$1,2,length(\$1)); getline; print; getline; getline}' > ${name}_$i.fa ]);
Utils::CMD(q[$$self{blat} adapters.fa ${name}_$i.fa ${name}_$i.blat -out=blast8]);
Utils::CMD(q[cat ${name}_$i.blat | awk '{if (\$2 ~ /^ADAPTER/) print \$1}' | sort -u | wc -l > ${name}_$i.nadapters]);
unlink("${name}_$i.fa", "${name}_$i.blat");
];
        close($fh);
        LSF::run($lock_file,$work_dir,"_${name}_a$i",{bsub_opts=>$$self{bsub_opts_process_fastqs}},qq[perl -w ${prefix}blat_fastq_$i.pl]);
    }

    return $$self{'No'};
}



#----------- map_sample ---------------------


sub map_sample_requires
{
    my ($self,$lane_path) = @_;

    my @requires = ();
    my $sample_dir = $$self{'sample_dir'};
    my $name = $$self{lane};

    my $fastq_files = existing_fastq_files("$lane_path/$sample_dir/$name");
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
    my ($self) = @_;
    my $sample_dir = $$self{'sample_dir'};
    my @provides = ("$sample_dir/$$self{lane}.bam");
    return \@provides;
}

sub map_sample
{
    my ($self,$lane_path,$lock_file) = @_;

    if ( !$$self{bwa_ref} ) { $self->throw("Missing the option bwa_ref.\n"); }

    my $sample_dir = $$self{'sample_dir'};
    my $name       = $$self{lane};
    my $work_dir   = "$lane_path/$$self{sample_dir}";

    my $bwa        = $$self{'bwa_exec'};
    my $bwa_ref    = $$self{'bwa_ref'};
    my $fa_ref     = $$self{'fa_ref'};
    my $fai_ref    = $$self{'fai_ref'};
    my $samtools   = $$self{'samtools'};
    my $max_isize  = exists($$self{insert_size}) ? $$self{insert_size}*2 : 2000;


    # How many files do we have?
    my $fastq_files = existing_fastq_files("$work_dir/$name");
    my $nfiles = scalar @$fastq_files;
    if ( $nfiles<1 || $nfiles>3 ) { $self->throw("FIXME: we can handle 1-3 fastq files in $work_dir, not $nfiles.\n") }

    # If there are multiple fastqs, assume that the data come from the paired-end sequencing
    if ( !exists($$self{paired}) && $nfiles>1 ) { $$self{paired}=1; }

    # There can be a mixture of: 
    #   1) two paired-end fastqs, no single
    #   2) one single+two paired-end fastqs 
    #   3) arbitrary number of singles (i.e. 1-3)
    my @singles = ();
    my @paired;
    if ( $nfiles==1 || !$$self{paired} )
    {
        for my $file (@$fastq_files) { push @singles, $file; }
    }
    elsif ( $nfiles==3 ) 
    { 
        push @singles, $$fastq_files[2]; 
    }
    if ( $$self{paired} && $nfiles>1 )
    {
        @paired = ($$fastq_files[0],$$fastq_files[1]);  # paired-end must be named _1 and _2
    }

    # Dynamic script to be run by LSF. We must check that the bwa exists alright
    #   - samtools do not return proper status and create .bam file even when there was
    #   nothing read from the input.
    #
    open(my $fh,'>', "$work_dir/_map.pl") or Utils::error("$work_dir/_map.pl: $!");
    print $fh "use Utils;\n";
    my @single_bams;
    for my $single (@singles)
    {
        if ( !($single=~m{/?(${name}_\d+).fastq.gz$})) { Utils::error("Failed to parse [$single]"); }
        my $sname = $1;
        my $bwa_cmd = "$bwa samse $bwa_ref $sname.sai $sname.fastq.gz";
        push @single_bams, "$sname.bam";

        print $fh
qq[
Utils::CMD("$bwa_cmd > $sname.sam");
if ( ! -s "$sname.sam" ) { Utils::error("The command ended with an error:\n\t$bwa_cmd > $sname.sam\\n"); }
Utils::CMD("$samtools import $fai_ref $sname.sam $sname.ubam");
Utils::CMD("$samtools sort $sname.ubam ${sname}x");
Utils::CMD("$samtools calmd -b ${sname}x.bam $fa_ref >${sname}.bam.part 2>/dev/null");
Utils::CMD("rm -f $sname.sam $sname.ubam ${sname}x.bam");
rename("${sname}.bam.part", "$sname.bam") or Utils::error("rename ${sname}.bam.part $sname.bam: \$!");
];
    }
    my $paired_name;
    if ( scalar @paired )
    {
        my $bwa_cmd  = "$bwa sampe -a $max_isize $bwa_ref ${name}_1.sai ${name}_2.sai ${name}_1.fastq.gz ${name}_2.fastq.gz";
        $paired_name = scalar @singles ? "tmp_$name" : $name;

        print $fh
qq[
Utils::CMD("$bwa_cmd > $paired_name.sam");
if ( ! -s "$paired_name.sam" ) { Utils::error("The command ended with an error:\\n\\t$bwa_cmd > $paired_name.sam\\n"); }
Utils::CMD("$samtools import $fai_ref $paired_name.sam $paired_name.ubam");
Utils::CMD("$samtools sort $paired_name.ubam ${paired_name}x");
Utils::CMD("$samtools calmd -b ${paired_name}x.bam $fa_ref >${paired_name}.bam.part 2>/dev/null");
Utils::CMD("rm -f $paired_name.sam $paired_name.ubam ${paired_name}x.bam");
rename("${paired_name}.bam.part", "$paired_name.bam") or Utils::error("rename ${paired_name}.bam.part $paired_name.bam: \$!");
];
    }

    # Now merge the files if necessary
    if ( scalar @single_bams + scalar @paired > 1 )
    {
        my $bams = join(' ', @single_bams);
        if ( $paired_name ) { $bams .= " $paired_name.bam"; }

        print $fh 
qq[
Utils::CMD("$samtools merge x$name.bam $bams");
if ( ! -s "x$name.bam" ) { Utils::error("The command ended with an error:\\n\\t$samtools merge x$name.bam $bams\\n"); }
rename("x$name.bam","$name.bam") or Utils::error("rename x$name.bam $name.bam: \$!");
];
    }
    elsif ( scalar @single_bams == 1 && $single_bams[0] ne "$name.bam" )
    {
        print $fh qq[rename('$single_bams[0]',"$name.bam") or Utils::error("rename $single_bams[0] $name.bam: \$!");\n];
    }
    close($fh);

    LSF::run($lock_file,$work_dir,"_${name}_sampe",{bsub_opts=>$$self{bsub_opts_map_sample}}, q{perl -w _map.pl});
    return $$self{'No'};
}

#----------- transposon ---------------------

sub transposon_requires
{
  my ($self) = @_;

  my $sample_dir = $$self{'sample_dir'};
  my @requires = ("$sample_dir/$$self{lane}_1.fastq.gz");
  return \@requires;
}

sub transposon_provides
{
  my ($self) = @_;
  my @provides = ();
  
  if(( defined $$self{reads_contain_transposon}) && $$self{reads_contain_transposon})
  {
    my $sample_dir = $$self{'sample_dir'};
    @provides = ("$sample_dir/$$self{lane}.transposon");
  }
  
  return \@provides;
}

sub transposon
{
   my ($self,$lane_path,$lock_file) = @_;
   my $sample_dir = $$self{'sample_dir'};
   
   if(( defined $$self{reads_contain_transposon}) && $$self{reads_contain_transposon})
   {
     my $tag_length = $$self{transposon_length} || 7;
     my $output_file = "$lane_path/$sample_dir/".$$self{lane}.".transposon";
     
     my $transposon_sequence_str;
     if(defined $$self{transposon_sequence})
     { 
       $transposon_sequence_str = 'tag => "'.$$self{transposon_sequence}.'"';
     }
   
     # Dynamic script to be run by LSF.
     open(my $fh, '>', "$lane_path/$sample_dir/_transposon.pl") or Utils::error("$lane_path/$sample_dir/_transposon.pl: $!");
     print $fh 
     qq[ 
          use strict;
          use warnings;
          use Pathogens::Parser::Transposon;
          my \$transposon = Pathogens::Parser::Transposon->new(
            'filename'   => '$$self{lane}_1.fastq.gz',
            'tag_length' => $tag_length,
            $transposon_sequence_str
          );
          open(OUT,'+>','$output_file') or die "couldnt open transposon output file \$!";
          print OUT \$transposon->percentage_reads_with_tag;
          close(OUT);
     ];
     close $fh;

     LSF::run($lock_file,"$lane_path/$sample_dir","_$$self{lane}_transposon", $self, qq{perl -w _transposon.pl});
   }
   return $$self{'No'};
}


#----------- stats_and_graphs ---------------------

sub stats_and_graphs_requires
{
    my ($self) = @_;
    my $sample_dir = $$self{'sample_dir'};
    my @requires = ("$sample_dir/$$self{lane}.bam");
    return \@requires;
}

sub stats_and_graphs_provides
{
    my ($self) = @_;
    my $sample_dir = $$self{'sample_dir'};
    my @provides = ("$sample_dir/_graphs.done","$sample_dir/$$self{lane}.cover");
    return \@provides;
}

sub stats_and_graphs
{
    my ($self,$lane_path,$lock_file) = @_;

    my $sample_dir = $$self{'sample_dir'};
    my $lane  = $$self{lane};
    my $stats_ref = exists($$self{stats_ref}) ? $$self{stats_ref} : '';

    # Get size of assembly
    my $vrtrack  = VRTrack::VRTrack->new($$self{db}) or $self->throw("Could not connect to the database: ",join(',',%{$$self{db}}),"\n");
    my $assembly = VRTrack::Assembly->new_by_name($vrtrack, $$self{assembly});
    my $reference_size = $assembly->reference_size();
    unless($reference_size){ $self->throw("Failed to find reference genome size for lane $lane_path\n"); }

    my $name = $$self{lane};
    my $vrlane      = VRTrack::Lane->new_by_hierarchy_name($vrtrack,$name) or $self->throw("No such lane in the DB: [$name]\n");
    my $insert_size = (VRTrack::Library->new($vrtrack, $vrlane->library_id())->insert_size() )*3 || 8000;


    # Dynamic script to be run by LSF.
    open(my $fh, '>', "$lane_path/$sample_dir/_graphs.pl") or Utils::error("$lane_path/$sample_dir/_graphs.pl: $!");
    print $fh 
qq[
use VertRes::Pipelines::TrackQC_Fastq;

my \%params = 
(
    'gc_depth_bin' => q[$$self{'gc_depth_bin'}],
    'mapviewdepth' => q[$$self{'mapviewdepth'}],
    'samtools'     => q[$$self{'samtools'}],
    'gcdepth_R'    => q[$$self{'gcdepth_R'}],
    'lane_path'    => q[$lane_path],
    'lane'         => q[$$self{lane}],
    'sample_dir'   => q[$$self{'sample_dir'}],
    'fa_ref'       => q[$$self{fa_ref}],
    'fai_ref'      => q[$$self{fai_ref}],
    'stats_ref'    => q[$stats_ref],
    'bwa_clip'     => q[$$self{bwa_clip}],
    'chr_regex'    => q[$$self{chr_regex}],
    'bwa_exec'     => q[$$self{bwa_exec}],
    'do_samtools_rmdup' => q[$$self{do_samtools_rmdup}]
);

my \$qc = VertRes::Pipelines::TrackQC_Fastq->new(\%params);
\$qc->run_graphs(\$params{lane_path}, $reference_size, $insert_size);
];
    close $fh;

    LSF::run($lock_file,"$lane_path/$sample_dir","_${lane}_graphs", {bsub_opts=>$$self{bsub_opts_stats_and_graphs}}, qq{perl -w _graphs.pl});
    return $$self{'No'};
}


sub run_graphs
{
    my ($self,$lane_path, $reference_size,$insert_size) = @_;

    $self->SUPER::run_graphs($lane_path,$insert_size);

    # Get coverage, depth and sd.
    my $bamcheck_file = qq[$lane_path/$$self{'sample_dir'}/$$self{'lane'}.bam.bc];
    my $cover_file    = qq[$lane_path/$$self{'sample_dir'}/$$self{'lane'}.cover];

    my $genomecover = Pathogens::Parser::GenomeCoverage->new( bamcheck => $bamcheck_file,
                                                              ref_size => $reference_size );
    my($coverage, $depth, $depth_sd) = $genomecover->coverage_depth();

    # Output cover file.
    open(my $cov_fh, "> $cover_file") or $self->throw("Cannot open: $cover_file\n");
    printf $cov_fh "[%d, %.2f, %.2f]\n", $coverage, $depth, $depth_sd;
    close($cov_fh);
}

1;

