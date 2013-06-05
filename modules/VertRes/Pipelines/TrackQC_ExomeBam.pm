=head1 NAME

VertRes::Pipelines::TrackQC_ExomeBam - pipeline for QC of exome BAM files, inherits from VertRes::Pipelines::TrackQC_Bam.

=head1 SYNOPSIS

# Make a conf file with database information and mapping/snp info.

# A template conf file is:
root    => '/path/to/data/hirearchy/root',
module  => 'VertRes::Pipelines::TrackQC_ExomeBam',
prefix  => '_',
log     => '/path/to/file.log',
max_failures => 20,

db  => 
{
    database => 'dbname',
    host     => 'dbhost',
    port     => 3306,
    user     => 'dbuser',
},

data => 
{
    exit_on_errors=>0,
    #bamcheck => 'bamcheck -q 20',

    db  => 
    {
        database => 'dbname',
        host     => 'dbhost',
        port     => 3306,
        user     => 'dbuser',
        password => 'xxx',
    },

    bwa_ref => '/abs/path/to/bwa_ref.fasta',
    fa_ref  => '/abs/path/to/ref.fasta',
    fai_ref => '/abs/path/to/ref.fasta.fai',
    
    glf              => '/path/to/glf',
    snps             => '/foo/genotype.bin',
    snp_sites        => '/foo/genotype.snp_sites',
    gtype_confidence => 1.2,
    exome_design => 'name of exome design',
    exome_coords => '/abs/path/to/coords/file',
},

# Make a pipeline file:
echo '__VRTrack_QC__  exome_qc.conf' > exome_qc.pipeline

# run the pipeline:
run-pipeline -c exome_qc.pipeline

# (and make sure it keeps running by adding that last to a regular cron job)

=head1 AUTHOR

Martin Hunt: mh12@sanger.ac.uk

=cut

package VertRes::Pipelines::TrackQC_ExomeBam;
use base qw(VertRes::Pipelines::TrackQC_Bam);

use strict;
use warnings;
use VertRes::LSF;
use File::Spec;
use Data::Dumper;
use Utils;


our @actions =
(
    # Takes care of merging of the (possibly) multiple bam files.  Innherited from TrackQC_Bam
    {
        'name'     => 'rename_and_merge',
        'action'   => \&VertRes::Pipelines::TrackQC_Bam::rename_and_merge,
        'requires' => \&VertRes::Pipelines::TrackQC_Bam::rename_and_merge_requires, 
        'provides' => \&VertRes::Pipelines::TrackQC_Bam::rename_and_merge_provides,
    },

    # Check sanity, namely the reference sequence.  Inherited from TrackQC_Bam
    {
        'name'     => 'check_sanity',
        'action'   => \&VertRes::Pipelines::TrackQC_Bam::check_sanity,
        'requires' => \&VertRes::Pipelines::TrackQC_Bam::check_sanity_requires, 
        'provides' => \&VertRes::Pipelines::TrackQC_Bam::check_sanity_provides,
    },

    # Check the genotype. Inherited from TrackQC_Bam
    {
        'name'     => 'check_genotype',
        'action'   => \&VertRes::Pipelines::TrackQC_Bam::check_genotype,
        'requires' => \&VertRes::Pipelines::TrackQC_Bam::check_genotype_requires, 
        'provides' => \&VertRes::Pipelines::TrackQC_Bam::check_genotype_provides,
    },

    # Creates some QC graphs and generate some statistics.
    {
        'name'     => 'stats_and_graphs',
        'action'   => \&stats_and_graphs,
        'requires' => \&stats_and_graphs_requires, 
        'provides' => \&stats_and_graphs_provides,
    },

    # Checks the generated stats and attempts to auto pass or fail the lane.  NOT YET IMPLEMENTED:
    # we'll have to see what exome data look like
    {
        'name'     => 'auto_qc',
        'action'   => \&auto_qc,
        'requires' => \&auto_qc_requires, 
        'provides' => \&auto_qc_provides,
    },

    # Writes the QC status to the tracking database.
    {
        'name'     => 'update_db',
        'action'   => \&update_db,
        'requires' => \&update_db_requires, 
        'provides' => \&update_db_provides,
    },
);

our $options = {
    # Executables
    'bamcheck'        => 'bamcheck -q 20',
    'blat'            => '/software/pubseq/bin/blat',
    'samtools'        => 'samtools',
    'clean_fastqs'    => 0,

    'adapters'        => '/software/pathogen/projects/protocols/ext/solexa-adapters.fasta',
    'bsub_opts'       => "-q normal -M5000 -R 'select[type==X86_64] select[mem>5000] rusage[mem=5000] rusage[thouio=1]'",
    'bsub_opts_merge' => "-q normal -M5000 -R 'select[type==X86_64] select[mem>5000] rusage[mem=5000] rusage[thouio=5]'",
    'bsub_opts_stats' => "-q normal -M3500 -R 'select[type==X86_64] select[mem>3500] rusage[mem=3500]'",
    'gc_depth_bin'    => 20000,
    'gtype_confidence'=> 5.0,
    'mapstat_id'      => 'mapstat_id.txt',
    'sample_dir'      => 'qc-sample',
    'stats'           => '_stats',
    'stats_detailed'  => '_detailed-stats.txt',
    'chr_regex'       => '^(?:\d+|X|Y)$',
};


# --------- OO stuff --------------

=head2 new

        Example    : my $qc = VertRes::Pipelines::TrackQC_ExomeBam->new('sample_dir'=>'dir');
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
                    bsub_opts_stats .. LSF bsub options for the stats_and_graphs script
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
                    qc              .. where to put qc files
                    sample_size     .. the size of the subsample (approx number of bases in one fastq file)
                    skip_genotype   .. boolean (default false).  If true, skip genotype check, and
                                       snps and snp_sites are not required
                    stats_ref       .. e.g. /path/to/NCBI36.stats 


                    # exome options
                    exome_coords    .. REQUIRED. File containing info on targets and baits, made by
                                       VertRes::Utils::Sam::bam_exome_qc_stats() with the
                                       'dump_intervals' option
                    exome_design    .. REQUIRED. Name of exome design.  This is what will appear in the 
                                       exome_design table in the database
                    snps            .. REQUIRED. File made by hapmap2bin to be used for genotyping, containing
                                       genotypes of all the samples
                    snp_sites       .. REQUIRED. File containing locations of SNPs to be used for genotyping. 
                                       One site per line, tab separated:  chromosome<tab>position

=cut

sub VertRes::Pipelines::TrackQC_ExomeBam::new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(%$options,'actions'=>\@actions,@args);
    $self->write_logs(1);
    if ( !$$self{samtools} ) { $self->throw("Missing the option samtools.\n"); }
    if ( !$$self{fa_ref} ) { $self->throw("Missing the option fa_ref.\n"); }
    if ( !$$self{fai_ref} ) { $self->throw("Missing the option fai_ref.\n"); }
    if ( !$$self{sample_dir} ) { $self->throw("Missing the option sample_dir.\n"); }

    unless ( $self->{skip_genotype} ) {
        if ( !$self->{snps} ) { $self->throw("Missing the option snps\n"); }
        if ( !$self->{snp_sites} ) { $self->throw("Missing the option snp_sites\n"); }
        if ( !$$self{glf} ) { $self->throw("Missing the option glf.\n"); }
        if ( !$$self{gtype_confidence} ) { $self->throw("Missing the option gtype_confidence.\n"); }
    }

    if ( !$self->{exome_coords} ) { $self->throw("Missing the option exome_coords.\n"); }
    if ( !$self->{exome_design} ) { $self->throw("Missing the option exome_design.\n"); }

    return $self;
}


sub stats_and_graphs_requires
{
    my ($self) = @_;
    my $sample_dir = $$self{'sample_dir'};
    my @requires = (File::Spec->catfile($sample_dir, "$$self{lane}.bam"));
    return \@requires;
}


sub stats_and_graphs_provides
{
    my ($self) = @_;
    my $sample_dir = $$self{'sample_dir'};
    my @provides = (File::Spec->catfile($sample_dir, "_stats_and_graphs.done"));
    return \@provides;
}


sub stats_and_graphs
{
    my ($self,$lane_path,$lock_file) = @_;
    my $sample_dir = $self->{'sample_dir'};
    my $lane = $self->{lane};
    my $outdir = File::Spec->catdir($lane_path, $sample_dir);
    my $bam = File::Spec->catfile($outdir, $lane . '.bam');
    my $qc_files_prefix = 'exomeQC';
    
    #Get current bc file
	my $bcfile;
	if ( -e "$outdir/$$self{mapstat_id}" )
    {
        my $mapstats_id = `cat $outdir/$$self{mapstat_id}`;
        chomp $mapstats_id;
    	my @bcfiles = glob("$lane_path/$mapstats_id*.bam.bc");
    	for my $file ( @bcfiles ) {
    		if ( $file=~m{$mapstats_id\.[ps]e\.raw\.sorted\.bam\.bc} ) {
    			$bcfile = $file;
    		}
    	}
    }
	#what if $bcfile is uninitialised - can it ever be at this stage? Run bamcheck?
	my $bc = VertRes::Parser::bamcheck->new(file=>"$bcfile");
	print STDERR "CHECK: ", $bc->get('sequences'), ", ", $bc->get('total_length'), ", ", $bc->get('reads_mapped'), "\n";
    if ( $bc->get('sequences') > 0 && $bc->get('total_length') > 0 && $bc->get('reads_mapped') > 0 )
    {
    	#my $stats_ref = exists($$self{stats_ref}) ? $$self{stats_ref} : '';

    	# Dynamic script to be run by LSF.
    	my $script = File::Spec->catfile($lane_path, $sample_dir, "_stats_and_graphs.pl");
    	#my $fakefile = '/nfs/users/nfs_m/mh12/Random/fake_exome_qc_dump';
    	open my $fh, '>', $script or Utils::error("$script: $!");
    	print $fh 
	qq[
use strict;
use warnings;
use VertRes::Utils::Sam;
use Data::Dumper;

my \$o = VertRes::Utils::Sam->new();
my \%opts = (
    load_intervals_dump_file => q[$self->{exome_coords}],
    bam => q[$bam],
    verbose => 1,
);
my \$dump_file = q[$qc_files_prefix] . q[.dump];
my \$stats;

if (-e \$dump_file) {
    \$stats = do \$dump_file or die("Could not load stats from file \$dump_file");
}
else {
    \$stats = \$o->bam_exome_qc_stats(\%opts);
    open my \$fh, '>', \$dump_file or die "error opening \$dump_file";
    print \$fh  Dumper \$stats;
    close \$fh;
}
\$o->bam_exome_qc_make_plots(\$stats, q[$qc_files_prefix], q[png]);
die "error touching done file" if (system "touch _stats_and_graphs.done");
];

    	close $fh;

    	my $orig_bsub_opts = $self->{bsub_opts};
    	$self->{bsub_opts} = $self->{bsub_opts_stats};
    	VertRes::LSF::run($lock_file,$outdir,"_${lane}_stats_and_graphs", $self, qq{perl -w _stats_and_graphs.pl});
    	$self->{bsub_opts} = $orig_bsub_opts;
    }
    else {
    	`touch $outdir/_stats_and_graphs.done`;
    }
    return $$self{'No'};
}

sub auto_qc_requires
{
    my ($self) = @_;
    my $sample_dir = $$self{'sample_dir'};
    my $name = $$self{lane};
    my @requires = (File::Spec->catfile($sample_dir, "_stats_and_graphs.done"));
    return \@requires;
}

# See description of update_db.
#
sub auto_qc_provides {
    my ($self) = @_;

    if ( exists($$self{db}) ) { return 0; }

    my @provides = ();
    return \@provides;
}

sub auto_qc {
    my ($self,$lane_path,$lock_file) = @_;

    my $sample_dir = "$lane_path/$$self{sample_dir}";
    if ( !$$self{db} ) { $self->throw("Expected the db key.\n"); }

    my $vrtrack   = VRTrack::VRTrack->new($$self{db}) or $self->throw("Could not connect to the database: ",join(',',%{$$self{db}}),"\n");
    my $name      = $$self{lane};
    my $vrlane    = VRTrack::Lane->new_by_hierarchy_name($vrtrack,$name) or $self->throw("No such lane in the DB: [$name]\n");

    if ( !$vrlane->is_processed('import') ) { return $$self{Yes}; }

    # NOT YET IMPLEMENTED ...

    return $$self{'Yes'};
}



sub update_db_requires {
    my ($self) = @_;
    my $sample_dir = $$self{'sample_dir'};
    my $name = $$self{lane};
    my @requires = ("$sample_dir/_stats_and_graphs.done");
    return \@requires;
}

# This subroutine will check existence of the key 'db'. If present, it is assumed
#   that QC should write the stats and status into the VRTrack database. In this
#   case, 0 is returned, meaning that the task must be run. The task will change the
#   QC status from 'no_qc' to something else, therefore we will not be called again.
#
#   If the key 'db' is absent, the empty list is returned and the database will not
#   be written.
#
sub update_db_provides {
    my ($self) = @_;

    if ( exists($$self{db}) ) { return 0; }

    my @provides = ();
    return \@provides;
}

sub update_db {
    my ($self,$lane_path,$lock_file) = @_;
    my $sample_dir = File::Spec->catfile($lane_path, $self->{sample_dir});
    if ( !$self->{db} ) { $self->throw("Expected the db key.\n"); }

    # First check if the 'no_qc' status is still present. Another running pipeline
    #   could have queued the job a long time ago and the stats might have been
    #   already written.
    my $vrtrack   = VRTrack::VRTrack->new($$self{db}) or $self->throw("Could not connect to the database: ",join(',',%{$$self{db}}),"\n");
    my $name      = $$self{lane};
    my $vrlane    = VRTrack::Lane->new_by_hierarchy_name($vrtrack,$name) or $self->throw("No such lane in the DB: [$name]\n");

    if ( !$vrlane->is_processed('import') ) { return $$self{Yes}; }

    my $stats_file = File::Spec->catfile($sample_dir, 'exomeQC.dump');

    $vrtrack->transaction_start();

    # Now call the database API and fill the mapstats object with values
    my $mapping;
    my $has_mapstats = 0;
    my $mapstats_id;

    if ( -e "$sample_dir/$$self{mapstat_id}" )
    {
        # When run on bam files created by the mapping pipeline, reuse existing
        #   mapstats, so that the mapper and assembly information is not overwritten.
        ($mapstats_id) = `cat $sample_dir/$$self{mapstat_id}`;
        chomp($mapstats_id);
        $mapping = VRTrack::Mapstats->new($vrtrack, $mapstats_id);
        if ( $mapping )  { $has_mapstats=1; }
    }
    if ( !$mapping ) { $mapping = $vrlane->add_mapping(); }

    my $gtype;
    unless(exists $$self{'skip_genotype'} && $$self{'skip_genotype'}) {
        # Get genotype results
        $gtype = VertRes::Utils::GTypeCheck::get_status("$sample_dir/${name}.gtype");
    }
    else { # Skip genotype results
    $gtype = {
            status   => 'unchecked',
            expected => undef,
            found    => undef,
            ratio    => undef,
        };
    }

	if ( -e $stats_file ) {
		my $stats = do $stats_file or $self->throw("Could not load stats from file $stats_file");

    	# Write the mapstats values
    	$mapping->raw_reads($stats->{raw_reads});
    	$mapping->raw_bases($stats->{raw_bases});
    	$mapping->reads_mapped($stats->{reads_mapped});
    	$mapping->reads_paired($stats->{reads_paired});
    	$mapping->bases_mapped($stats->{bases_mapped});
    	$mapping->error_rate($stats->{pct_mismatches} / 100); # keep value in [0,1] to be consistent with standard QC
    	$mapping->rmdup_reads_mapped($stats->{rmdup_reads_mapped});
    	$mapping->rmdup_bases_mapped($stats->{rmdup_bases_mapped});
    	# TODO: need adapter search from BAM file. $mapping->adapter_reads($nadapters);
    	$mapping->clip_bases($stats->{raw_bases} - $stats->{clip_bases});
    	$mapping->mean_insert($stats->{mean_insert_size});
    	$mapping->sd_insert($stats->{insert_size_sd});
    	$mapping->bait_near_bases_mapped($stats->{bait_near_bases_mapped});
    	$mapping->target_near_bases_mapped($stats->{target_near_bases_mapped});
    	$mapping->bait_bases_mapped($stats->{bait_bases_mapped});
    	$mapping->mean_bait_coverage($stats->{mean_bait_coverage});
    	$mapping->bait_coverage_sd($stats->{bait_coverage_sd});
    	$mapping->off_bait_bases($stats->{off_bait_bases});
    	$mapping->reads_on_bait($stats->{reads_on_bait});
    	$mapping->reads_on_bait_near($stats->{reads_on_bait_near});
    	$mapping->reads_on_target($stats->{reads_on_target});
    	$mapping->reads_on_target_near($stats->{reads_on_target_near});
    	$mapping->target_bases_mapped($stats->{target_bases_mapped});
    	$mapping->mean_target_coverage($stats->{mean_target_coverage});
    	$mapping->target_coverage_sd($stats->{target_coverage_sd});
    	$mapping->target_bases_1X($stats->{target_bases_1X});
    	$mapping->target_bases_2X($stats->{target_bases_2X});
    	$mapping->target_bases_5X($stats->{target_bases_5X});
    	$mapping->target_bases_10X($stats->{target_bases_10X});
    	$mapping->target_bases_20X($stats->{target_bases_20X});
    	$mapping->target_bases_50X($stats->{target_bases_50X});
    	$mapping->target_bases_100X($stats->{target_bases_100X});
		# update the exome_design		
		my $exome_design = $mapping->exome_design($self->{exome_design});
    	unless ($exome_design) { $exome_design = $mapping->add_exome_design($self->{exome_design}) }
    	$exome_design->bait_bases($stats->{bait_bases});
    	$exome_design->target_bases($stats->{target_bases});
    	$exome_design->update();
    	$mapping->exome_design_id($exome_design->id());
    }
    # TODO: genotyping
    $mapping->genotype_expected($$gtype{expected});
    $mapping->genotype_found($$gtype{found});
    $mapping->genotype_ratio($$gtype{ratio});
    $vrlane->genotype_status($$gtype{status});

    # If there is no mapstats present, the mapper and assembly must be filled in.
    $self->_update_mapper_and_assembly($sample_dir, $mapping) unless $has_mapstats;
	$mapping->update;

    # Write the images
    my %images;
    my $prefix = 'exomeQC';
    if (-e File::Spec->catfile($sample_dir, $prefix . '.bait_gc_vs_cvg.png'))      {$images{$prefix . '.bait_gc_vs_cvg.png'} = 'Bait GC vs coverage'};
    if (-e File::Spec->catfile($sample_dir, $prefix . '.coverage_per_base.png'))   {$images{$prefix . '.coverage_per_base.png'} = 'Bait and target coverage per base'};
    if (-e File::Spec->catfile($sample_dir, $prefix . '.cumulative_coverage.png')) {$images{$prefix . '.cumulative_coverage.png'} = 'Cumulative coverage distribution'};
    if (-e File::Spec->catfile($sample_dir, $prefix . '.gc_mapped.png'))           {$images{$prefix . '.gc_mapped.png'} = 'GC of mapped reads'};
    if (-e File::Spec->catfile($sample_dir, $prefix . '.gc_unmapped.png'))         {$images{$prefix . '.gc_unmapped.png'} = 'GC of unmapped reads'};
    if (-e File::Spec->catfile($sample_dir, $prefix . '.insert_size.png'))         {$images{$prefix . '.insert_size.png'} = 'Insert size distribution'};
    if (-e File::Spec->catfile($sample_dir, $prefix . '.mean_coverage.png'))       {$images{$prefix . '.mean_coverage.png'} = 'Mean bait and target coverage distribution'};
    if (-e File::Spec->catfile($sample_dir, $prefix . '.normalised_coverage.png')) {$images{$prefix . '.normalised_coverage.png'} = 'Normalised coverage distribution'};
    if (-e File::Spec->catfile($sample_dir, $prefix . '.quality_scores_1.png'))    {$images{$prefix . '.quality_scores_1.png'} = 'Quality scores 1'};
    if (-e File::Spec->catfile($sample_dir, $prefix . '.quality_scores_2.png'))    {$images{$prefix . '.quality_scores_2.png'} = 'Quality scores 2'};
    if (-e File::Spec->catfile($sample_dir, $prefix . '.target_gc_vs_cvg.png'))    {$images{$prefix . '.target_gc_vs_cvg.png'} = 'Target GC vs coverage'};
    if (-e File::Spec->catfile($sample_dir, $prefix . '.target_gc_vs_cvg.scaled.png'))    {$images{$prefix . '.target_gc_vs_cvg.scaled.png'} = 'Target scaled GC vs coverage'};
    if (-e File::Spec->catfile($sample_dir, $prefix . '.bait_gc_vs_cvg.scaled.png'))    {$images{$prefix . '.bait_gc_vs_cvg.scaled.png'} = 'Bait scaled GC vs coverage'};

    while (my ($imgname,$caption) = each %images) {
        my $img = $mapping->add_image_by_filename(File::Spec->catfile($sample_dir, $imgname));
        $img->caption($caption);
        $img->update;
    }

    # Write the QC status. Never overwrite a QC status set previously by human. Only NULL or no_qc can be overwritten.
    $self->_write_QC_status($vrtrack, $vrlane, $name);
    $vrtrack->transaction_commit();

    # Clean the big files
    for my $file ('gc-depth.bindepth',"$$self{lane}.bam.bai","$$self{lane}*.sai","$$self{lane}*.fastq.gz","$$self{lane}.bam","$$self{lane}.glf")
    {
        Utils::CMD("rm -f $sample_dir/$file");
    }

    if ( $$self{clean_fastqs} )
    {
        Utils::CMD("rm -f $lane_path/$$self{lane}*.fastq.gz");
    }

    # Rename the sample dir by mapstats ID, cleaning existing one
    if ( $has_mapstats && $mapstats_id )
    {
        # Check if the sample dir is a symlink to reused $sample_dir.$mapstats_id directory.
        #   In that case keep the .id directory and remove $sample_dir. Otherwise rename
        #   $sample_dir to the .id dir.
        my $link;
        if ( -l $sample_dir && ($link = readlink($sample_dir)) && $link=~/\.$mapstats_id$/ )
        {
            unlink($sample_dir);
        }
        else
        {
            Utils::CMD("rm -rf $sample_dir.$mapstats_id");
            rename($sample_dir,"$sample_dir.$mapstats_id");
        }
    }

    return $$self{'Yes'};
}


1;
