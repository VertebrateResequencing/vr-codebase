=head1 NAME

VertRes::Wrapper::GATK - wrapper for Broad's GenomeAnalysisToolKit

=head1 SYNOPSIS

use VertRes::Wrapper::GATK;

my $wrapper = VertRes::Wrapper::GATK->new();

# run something with the toolkit
$wrapper->count_covariates('in.bam', 'in.bam');
$wrapper->table_recalibration('in.bam', 'in.bam.recal_data.csv', 'out.bam');

# or for your convienience:
$wrapper->recalibrate('in.bam', 'out.bam');

# check the status
my $status = $wrapper->run_status;
if ($status == -1) {
    # try and run it again?...
}

=head1 DESCRIPTION

A wrapper for Broad's GenomeAnalysisToolKit, focusing on all the steps needed
to call SNPs (from recalibration, through realignment to genotyping). See:
http://www.broadinstitute.org/gsa/wiki/index.php/Whole_genome,_low-pass

For default "exe" path assumes you have the env variable GATK pointing to the
directory containing the GATK .jar files etc.
For default reference and dbsnp file paths assumes you have the env variable
GATK_RESOURCES pointing to the directory containing the .rod and .fa files etc.
We also expect a vcfs subdirectory to be added there containing standard vcf
files for use.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Wrapper::GATK;

use strict;
use warnings;
use Cwd qw (abs_path);
use File::Basename;
use File::Spec;
use File::Copy;

use base qw(VertRes::Wrapper::WrapperI);
use VertRes::Wrapper::samtools;
use VertRes::Utils::FileSystem;

our $DEFAULT_GATK_JAR = File::Spec->catfile($ENV{GATK}, 'GenomeAnalysisTK.jar');
our $DEFAULT_LOGLEVEL = 'ERROR';
our $DEFAULT_PLATFORM = 'ILLUMINA';
our $DEFAULT_RG = 'RG';
my $fsu = VertRes::Utils::FileSystem->new();

=head2 new

 Title   : new
 Usage   : my $wrapper = VertRes::Wrapper::GATK->new();
 Function: Create a VertRes::Wrapper::GATK object.
 Returns : VertRes::Wrapper::GATK object
 Args    : exe   => string (full path to GenomeAnalysisTK.jar; a TEAM145 default
                            exists)
           java_memory => int (the amount of memory in MB to give java; default
                               2800)
           reference => ref.fa (path to reference fasta; can be overriden in
                                individual methods with the R option)
           dbsnp     => snp.rod (path to dbsnp rod file; can be overriden in
                                 individual methods with the DBSNP option)
           covs      => [] (as per set_covs())
           bs        => [] (as per set_b())
           build     => NCBI36|NCBI37|NCBIM37 (default NCBI37: sets defaults for
                        reference, dbsnp, covs and bs as appropriate for the
                        build; overriden by the above 4 options if they are set
                        manually)
           log_level => DEBUG|INFO|WARN|ERROR|FATAL|OFF (set the log level;
                                can be overriden in individual methods with the
                                l option; default 'ERROR')
           default_platform => ILLUMINA (when not set in the RG header, this
                                         platform will be used)
           default_readgroup => string (when not set in the bam records, this
                                        readgroup will be used)
           tmp_dir => /tmp/dir (VertRes::Utils::FileSystem->tempdir by default;
                                Any supplied directory is used as the root for a
                                new directory that will be auto-deleted)

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(exe => $DEFAULT_GATK_JAR, @args);
    
    my $temp_dir = delete $self->{tmp_dir};
    $temp_dir = $fsu->tempdir($temp_dir ? (DIR => $temp_dir) : ());
    
    my $java_mem = delete $self->{java_memory} || 2800;
    my $xss = 280; # int($java_mem / 10);
    if ($java_mem > 1000) {
        $xss = "-Xss${xss}m";
    }
    else {
        # login node with small memory limit doesn't like Xss option at all
        $xss = '';
    }
    $self->{java_memory} = $java_mem;
    $self->exe("java -Xmx${java_mem}m -Xms${java_mem}m $xss -Djava.io.tmpdir=$temp_dir -server -XX:+UseParallelGC -XX:ParallelGCThreads=2 -jar ".$self->exe);
    
    # our bsub jobs will get killed if we don't select high-mem machines
    $self->bsub_options(M => ($java_mem * 1000), R => "'select[mem>$java_mem] rusage[mem=$java_mem]'");
    
    # default settings
    my $build = delete $self->{build} || '';
    my ($default_ref, $default_dbsnp, $default_covs, $default_bs);
    if ($build eq 'NCBI36') {
        $default_ref = File::Spec->catfile($ENV{GATK_RESOURCES}, 'human_b36_both.fasta');
        $default_dbsnp = File::Spec->catfile($ENV{GATK_RESOURCES}, 'dbsnp_129_b36.rod');
    }
    elsif ($build eq 'NCBI37') {
        $default_ref = File::Spec->catfile($ENV{GATK_RESOURCES}, 'human_g1k_v37.fasta');
        $default_dbsnp = File::Spec->catfile($ENV{GATK_RESOURCES}, 'dbsnp_130_b37.rod');
        # no default bs for now, since issues with not being in NCBI37 coords?
        # ... and even now we have them, they're population specific...
        # ... actually, we can and should use all the populations regardless of
        #     which pop we're analysing
        #$default_bs = ['pilot1_CEU,VCF,'.File::Spec->catfile($ENV{GATK_RESOURCES}, 'vcfs', 'CEU.2and3_way.vcf'),
        #               'pilot1_YRI,VCF,'.File::Spec->catfile($ENV{GATK_RESOURCES}, 'vcfs', 'YRI.2and3_way.vcf')];
    }
    elsif ($build eq 'NCBIM37') {
        # no defaults yet, could set them up if desired...
    }
    elsif ($build) {
        $self->throw("bad build option '$build'");
    }
    
    $self->{_default_R} = delete $self->{reference} || $default_ref;
    # Allow overriding the default DBSNP by undef
    #   $self->{_default_DBSNP} = delete $self->{dbsnp} || $default_dbsnp;
    $self->{_default_DBSNP} = exists($$self{dbsnp}) ? $$self{dbsnp} : $default_dbsnp; delete($$self{dbsnp});
    $self->set_covs(@{delete $self->{covs} || $default_covs || []});
    $self->set_b(@{delete $self->{bs} || $default_bs || []});
    $self->{_default_loglevel} = delete $self->{log_level} || $DEFAULT_LOGLEVEL;
    $self->{_default_platform} = delete $self->{default_platform} || $DEFAULT_PLATFORM;
    $self->{_default_read_group} = delete $self->{default_read_group} || $DEFAULT_RG;
    
    return $self;
}

sub _handle_common_params {
    my ($self, $params) = @_;
    
    unless (defined $params->{R}) {
        $params->{R} = $self->{_default_R};
    }
    unless (defined $params->{DBSNP}) {
        $params->{DBSNP} = $self->{_default_DBSNP};
    }
    unless (defined $params->{l}) {
        $params->{l} = $self->{_default_loglevel};
    }
    unless (defined $params->{default_platform}) {
        $params->{default_platform} = $self->{_default_platform};
    }
    unless (defined $params->{default_read_group}) {
        $params->{default_read_group} = $self->{_default_read_group};
    }
}

=head2 version

 Title   : version
 Usage   : $wrapper->version();
 Function: Get the version of GATK.
 Returns : string of GATK version
 Args    : n/a

=cut

sub version {
    my $self = shift;
    
    open(my $fh, "java -jar $DEFAULT_GATK_JAR -h 2>&1 |") || $self->throw("Could not start $DEFAULT_GATK_JAR");
    my $version = 0;
    while (<$fh>) {
        if (/v([\d\.]+)/) {
            $version = $1;
            last;
        }
    }
    close $fh ;
    
    return $version;
}

=head2 count_covariates

 Title   : count_covariates
 Usage   : $wrapper->count_covariates('input.bam', 'output_prefix');
 Function: Generates a file necessary for recalibration. Will create first
           input.bam.bai using samtools index if it doesn't already exist.
 Returns : n/a
 Args    : path to input .bam file, path to output file (which will have its
           name suffixed with '.recal_data.csv' unless you suffix it yourself
           with .csv). Optionally, supply R, DBSNP, useOriginalQualities or l
           options (as a hash), as understood by GATK. -B and -cov should be set
           with the set_b() and set_covs() methods beforehand.

=cut

sub count_covariates {
    my ($self, $in_bam, $out_csv, @params) = @_;
    
    # java -Xmx2048m -jar GenomeAnalysisTK.jar \
    #   -R resources/Homo_sapiens_assembly18.fasta \
    #   --DBSNP resources/dbsnp_129_hg18.rod \
    #   -B mask,VCF,sitesToMask.vcf \
    #   -l INFO \
    #   -T CountCovariates \ 
    #   -I my_reads.bam \
    #   -cov ReadGroupCovariate \
    #   -cov QualityScoreCovariate \
    #   -cov CycleCovariate \
    #   -cov DinucCovariate \
    #   -recalFile my_reads.recal_data.csv
    
    # CycleCovariate: The machine cycle for this base (different definition for the various technologies and therefore platform [@PL tag] is pulled out of the read's read group).
    # DinucCovariate: The combination of this base and the previous base.
    # HomopolymerCovariate: The number of consecutive previous bases that match the current base.
    # MappingQualityCovariate: The mapping quality assigned to this read by the aligner.
    # MinimumNQSCovariate: The minimum base quality score in a small window in the read around this base.
    # PositionCovariate: The position along the length of the read. For Illumina this is the same as machine cycle but that is not the case for the other platforms.
    # PrimerRoundCovariate: The primer round for this base (only meaningful for SOLiD reads).
    # QualityScoreCovariate: The reported base quality score for this base.
    # ReadGroupCovariate: The read group this read is a member of.
    
    $self->switches([qw(useOriginalQualities)]);
    $self->params([qw(R DBSNP l T default_platform default_read_group)]);
    
    # used to take a fileroot, but now takes an output file
    my $recal_file = $out_csv;
    unless ($recal_file =~ /\.csv$/) {
        $recal_file .= '.recal_data.csv';
    }
    my $covs = $self->get_covs();
    my $bs = $self->get_b();
    my @file_args = (" -I $in_bam", " $covs $bs -recalFile $recal_file");
    
    my %params = @params;
    $params{T} = 'CountCovariates';
    unless (defined $params{useOriginalQualities}) {
        $params{useOriginalQualities} = 1;
    }
    $self->_handle_common_params(\%params);
    
    $self->register_output_file_to_check($recal_file);
    $self->_set_params_and_switches_from_args(%params);
    
    return $self->run(@file_args);
}

=head2 set_covs

 Title   : set_covs
 Usage   : $wrapper->set_covs('CycleCovariate', 'PositionCovariate');
 Function: Set which covariates to calculate. See http://www.broadinstitute.org/gsa/wiki/index.php/Base_quality_score_recalibration#Available_covariates
           for a list of valid ones.
 Returns : list currently set
 Args    : list of covariate name strings

=cut

sub set_covs {
    my $self = shift;
    if (@_) {
        $self->{covs} = [@_];
    }
    return @{$self->{covs} || []};
}

=head2 get_covs

 Title   : get_covs
 Usage   : my $covs_string = $wrapper->get_covs();
 Function: Get the command line options defining which covarates will be used
           (one or more -cov options). If none have been set with set_cov,
           returns the special --standard_covs option.
 Returns : string
 Args    : n/a

=cut

sub get_covs {
    my $self = shift;
    my @covs = $self->set_covs;
    
    if (@covs) {
        my $args = '';
        foreach my $cov (@covs) {
            $args .= "-cov $cov ";
        }
        return $args;
    }
    else {
        return '--standard_covs';
    }
}

=head2 set_annotations

 Title   : set_annotations
 Usage   : $wrapper->set_annotations('HaplotypeScore', 'SB', 'QD');
 Function: Set which annotations to use during generate_variant_clusters().
 Returns : list currently set
 Args    : list of covariate name strings

=cut

sub set_annotations {
    my $self = shift;
    if (@_) {
        $self->{ans} = [@_];
    }
    return @{$self->{ans} || []};
}

=head2 get_annotations

 Title   : get_annotations
 Usage   : my $ans_string = $wrapper->get_annotations();
 Function: Get the command line options defining which annotations will be used
           (one or more -an options). 
 Returns : string
 Args    : n/a

=cut

sub get_annotations {
    my $self = shift;
    my @ans = $self->set_annotations;

    # Let the user decide what annotations are needed. GATK sets its own defaults anyway.
    # unless (@ans) {
    #     @ans = ('HaplotypeScore', 'SB', 'QD', 'HRun');
    # }
    
    my $args = '';
    foreach my $an (@ans) {
        $args .= "-an $an ";
    }
    
    return $args;
}

=head2 set_b

 Title   : set_b
 Usage   : $wrapper->set_b('pilot1,VCF,file1.vcf', 'pilot2,VCF,file2.vcf');
 Function: Set which variant files (VCF, BED, dnSNP) to use (for setting B
           option).
 Returns : list currently set
 Args    : list of variant strings (name,type,filename)

=cut

sub set_b {
    my $self = shift;
    if (@_) {
        $self->{bs} = [@_];
    }
    return @{$self->{bs} || []};
}

=head2 add_b

 Title   : add_b
 Usage   : $wrapper->add_b('pilot1,VCF,file1.vcf', 'pilot2,VCF,file2.vcf');
 Function: Add more variant files to use (having previously used set_b or the bs
           option to new()).
 Returns : n/a
 Args    : list of variant strings (name,type,filename)

=cut

sub add_b {
    my $self = shift;
    if (@_) {
        my @current_bs = $self->set_b;
        $self->set_b(@current_bs, @_);
    }
}

=head2 get_b

 Title   : get_b
 Usage   : my $covs_string = $wrapper->get_b();
 Function: Get the command line options defining which vcf or bed etc. files
           will be used (one or more -B options). If none have been set with
           set_b, returns empty string.
 Returns : string
 Args    : n/a

=cut

sub get_b {
    my $self = shift;
    my @bs = $self->set_b;
    
    my $args = '';
    foreach my $b (@bs) {
        # -B:variant,VCF snps.raw.vcf
        # -B:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.b37.sites.vcf
        # it used to be -B variant,VCF,snps.raw.vcf, so convert latter to former
        my @b = split(',', $b);
        my $last = pop @b;
        if ($last =~ /\s/) {
            ($b[++$#b], $last) = split(/\s/, $last);
        }
        $args .= "-B:".(join ',', @b)." $last";
    }
    return $args;
}

=head2 set_filters

 Title   : set_filters
 Usage   : $wrapper->set_filters(
                     filters => {'filter_name' => 'filter expression', ... },
                     g_filters => {'filter_name' => 'filter expression', ... });
 Function: Set which filters to use during variant_filtration().
 Returns : hash (keys are filters and g_filters) of hash refs of those filters
           currently set
 Args    : hash with keys 'filters' and 'g_filters', with values as hash refs.
           Those hash refs should have keys as filter names (setting
           --filterName or --genotypeFilterName) and values as filter
           expressions (setting --filterExpression and
           --genotypeFilterExpression)

=cut

sub set_filters {
    my $self = shift;
    if (@_) {
        my %hash = @_;
        $self->{the_filters} = {filters => $hash{filters} || {},
                                g_filters => $hash{g_filters} || {}};
    }
    return %{$self->{the_filters} || {}};
}

=head2 get_filters

 Title   : get_filters
 Usage   : my $filters_string = $wrapper->get_filters();
 Function: Get the command line options defining which filters to use for
           variant_filtration().
 Returns : string
 Args    : n/a

=cut

sub get_filters {
    my $self = shift;
    my %both_filters = $self->set_filters;
    
    my $args = '';
    foreach my $type ('filters', 'g_filters') {
        my $filters = $both_filters{$type} || next;
        my $name_arg = $type eq 'filters' ? '--filterName' : '--genotypeFilterName';
        my $exp_arg = $type eq 'filters' ? '--filterExpression' : '--genotypeFilterExpression';
        
        while (my ($name, $exp) = each %{$filters}) {
            $args .= " $exp_arg \"$exp\" $name_arg \"$name\"";
        }
    }
    return $args;
}

=head2 table_recalibration

 Title   : table_recalibration
 Usage   : $wrapper->table_recalibration('in.bam', 'recal_data.csv', 'out.bam');
 Function: Recalibrates a bam using the csv file made with count_covariates().
 Returns : n/a
 Args    : path to input .bam file, path to csv file made by count_covariates(),
           path to output file. Optionally, supply R or l or useOriginalQualities
           options (as a hash), as understood by GATK. useOriginalQualities is on
           by default.

=cut

sub table_recalibration {
    my ($self, $in_bam, $csv, $out_bam, @params) = @_;
    
    # java -Xmx2048m -jar GenomeAnalysisTK.jar \
    #   -l INFO \ 
    #   -R resources/Homo_sapiens_assembly18.fasta \ 
    #   -T TableRecalibration \
    #   -I my_reads.bam \
    #   -out my_reads.recal.bam \
    #   -recalFile my_reads.recal_data.csv
    
    $self->switches([qw(useOriginalQualities fail_with_no_eof_marker
                        doNotWriteOriginalQuals disable_bam_indexing
                        generate_md5 simplifyBAM)]);
    $self->params([qw(R l T default_platform default_read_group)]);
    
    my @file_args = (" -I $in_bam", " -recalFile $csv", " --out $out_bam");
    
    my %params = (fail_with_no_eof_marker => 1, disable_bam_indexing => 1, @params);
    $params{T} = 'TableRecalibration';
    unless (defined $params{useOriginalQualities}) {
        $params{useOriginalQualities} = 1;
    }
    $self->_handle_common_params(\%params);
    
    $self->register_output_file_to_check($out_bam);
    $self->_set_params_and_switches_from_args(%params);
    
    return $self->run(@file_args);
}

=head2 realignment_targets

 Title   : realignment_targets
 Usage   : $wrapper->set_b('indels,VCF,indel_calls.vcf');
           $wrapper->realignment_targets('in.bam', 'out.intervals');
 Function: Finds target intervals in a quality-recalibrated bam that could be
           realigned with indel_realigner().
 Returns : n/a
 Args    : path to input .bam file, path to output intervals. If doing
           indel_realigner on a constant set of known indels only, input bam
           is not required here, and should be set to undef.
           Optionally, supply R or DBSNP options etc (as a hash), as understood
           by GATK. Use set_b() if you have known snps/indels.

=cut

sub realignment_targets {
    my ($self, $in_bam, $out_intervals, @params) = @_;
    
    # java -jar GenomeAnalysisTK.jar \
    #   -T RealignerTargetCreator \
    #   -I recalibrated.bam \
    #   -R resources/Homo_sapiens_assembly18.fasta \
    #   -o forRealigner.intervals \
    #   -B:indels,VCF /path/to/indel_calls.vcf \
    #   -D resources/dbsnp_129_hg18.rod
    
    $self->switches([]);
    $self->params([qw(R DBSNP T minReadsAtLocus maxIntervalSize
                      mismatchFraction windowSize)]);
    
    my $bs = $self->get_b();
    my $I = $in_bam ? "-I $in_bam" : '';
    my @file_args = (" $bs $I", " -o $out_intervals");
    
    my %params = @params;
    $params{T} = 'RealignerTargetCreator';
    $self->_handle_common_params(\%params);
    
    $self->register_output_file_to_check($out_intervals);
    $self->_set_params_and_switches_from_args(%params);
    
    return $self->run(@file_args);
}

=head2 indel_realigner

 Title   : indel_realigner
 Usage   : $wrapper->set_b('indels,VCF,indel_calls.vcf');
           $wrapper->indel_realigner('in.bam', 'intervals', 'out.bam');
 Function: Does local realignment around indels in intervals as determined by
           realignment_targets(), generating a "cleaned" bam suitable for
           calling SNPs on.
 Returns : n/a
 Args    : path to input .bam file, path to output of realignment_targets(),
           path to output bam.
           Optionally, supply R or DBSNP options etc (as a hash), as understood
           by GATK. LODThresholdForCleaning and useOnlyKnownIndels are set by
           default. Use set_b() if you have known snps/indels. maxReadsInMemory is
           set to 100x java_memory, minimum 500000 by default.

=cut

sub indel_realigner {
    my ($self, $in_bam, $intervals_file, $out_bam, @params) = @_;
    
    #java -Xmx4g -Djava.io.tmpdir=/path/to/tmpdir \
    #    -jar /path/to/GenomeAnalysisTK.jar \
    #    -I <lane-level.bam> \
    #    -R <ref.fasta> \
    #    -T IndelRealigner \
    #    -targetIntervals <intervalListFromStep1Above.intervals> \
    #    -o <realignedBam.bam> \
    #    -B:indels,VCF /path/to/indel_calls.vcf \
    #    -D /path/to/dbsnp.rod \
    #    -LOD 0.4
    #    -compress 0
    
    $self->switches([qw(noOriginalAlignmentTags
                        disable_bam_indexing generate_md5 simplifyBAM
                        noPGTag targetIntervalsAreNotSorted)]);
    $self->params([qw(R DBSNP T maxReadsForConsensuses maxConsensuses
                      entropyThreshold bam_compression 
                      consensusDeterminationModel
                      maxIsizeForMovement maxPositionalMoveAllowed
                      maxReadsForRealignment
                      LODThresholdForCleaning maxReadsInMemory)]);
    
    my $bs = $self->get_b();
    my @file_args = (" $bs -I $in_bam", " -targetIntervals $intervals_file", " -o $out_bam");
    
    my %params = (LODThresholdForCleaning => 0.4, bam_compression => 0, disable_bam_indexing => 1, 
                    consensusDeterminationModel => 'KNOWNS_ONLY', @params);
    if (! defined $params{maxReadsInMemory}) {
        $params{maxReadsInMemory} = 100 * ($self->{java_memory} >= 5000 ? $self->{java_memory} : 5000);
    }
    $params{T} = 'IndelRealigner';
    $self->_handle_common_params(\%params);
    
    $self->register_output_file_to_check($out_bam);
    $self->_set_params_and_switches_from_args(%params);
    
    return $self->run(@file_args);
}

=head2 indel_genotyper

 Title   : indel_genotyper
 Usage   : $wrapper->indel_genotyper('in.bam', 'out.raw.bed', 'out.detailed.bed');
 Function: Call indels on a bam file (preferably one that has been output by
           indel_realigner()).
 Returns : n/a
 Args    : path to input .bam file, paths to two bed output files.
           Optionally, supply R or DBSNP options (as a hash), as understood by
           GATK, along with the other options like minIndelCount etc (1000
           genomes defaults exist).

=cut

sub indel_genotyper {
    my ($self, $in_bam, $out_raw_bed, $out_detailed_bed, @params) = @_;
    
    # java -jar GenomeAnalysisTK.jar \
    #   -T IndelGenotyperV2 \
    #   -R resources/Homo_sapiens_assembly18.fasta \
    #   -I cleaned.bam \
    #   -O indels.raw.bed \
    #   -o detailed.output.bed \
    #   --verbose \
    #   -minCnt 2 \   (minIndelCount)
    #   -minFraction 0.03 \
    #   -minConsensusFraction 0.6 \
    #   -mnr 1000000  (maxNumberOfReads)
    
    $self->switches([qw(verbose)]);
    $self->params([qw(R DBSNP T L 1kg_format minCoverage minNormalCoverage
                      minFraction minConsensusFraction minIndelCount refseq
                      blacklistedLanes window_size maxNumberOfReads)]);
    
    my @file_args = (" -I $in_bam", " -bed $out_raw_bed -o $out_detailed_bed");
    
    my %params = (minIndelCount => 2, minFraction => 0.03,
                  minConsensusFraction => 0.6,
                  window_size => 324, @params);
    $params{T} = 'SomaticIndelDetector';
    $self->_handle_common_params(\%params);
    
    $self->register_output_file_to_check($out_raw_bed);
    $self->register_output_file_to_check($out_detailed_bed);
    $self->_set_params_and_switches_from_args(%params);
    
    return $self->run(@file_args);
}

=head2 unified_genotyper

 Title   : unified_genotyper
 Usage   : $wrapper->unified_genotyper('in.bam', 'out.vcf');
 Function: Call SNPs on a bam file (preferably one that has been output by
           indel_realigner()).
 Returns : n/a
 Args :    path to input .bam file (or bam fofn), path to output vcf file.
           Optionally, supply R or DBSNP options (as a hash), as understood by
           GATK, along with the other options like confidence etc (1000
           genomes defaults exist).

=cut

sub unified_genotyper {
    my ($self, $in_bam, $out_vcf, @params) = @_;
    
    # java -jar GenomeAnalysisTK.jar \
    #   -R resources/Homo_sapiens_assembly18.fasta \
    #   -T UnifiedGenotyper \
    #   -I cleaned.bam \
    #   -D resources/dbsnp_129_hg18.rod \
    #   -o snps.raw.vcf \
    #   --standard_min_confidence_threshold_for_calling 10.0
    
    $self->switches([qw(noSLOD)]);
    $self->params([qw(R DBSNP T L confidence genotype_likelihoods_model heterozygosity
                      standard_min_confidence_threshold_for_calling
                      standard_min_confidence_threshold_for_emitting
                      assume_single_sample_reads platform metrics_file
                      genotyping_mode output_mode pcr_error_rate
                      min_indel_count_for_genotyping indel_heterozygosity
                      min_base_quality_score min_mapping_quality_score
                      max_deletion_fraction debug_file annotation group)]);
    
    my @file_args = (" -I $in_bam", " -o $out_vcf ");
    
    my %params = (standard_min_confidence_threshold_for_calling => 10.0, @params);
    $params{T} = 'UnifiedGenotyper';
    $self->_handle_common_params(\%params);
    
    $self->register_output_file_to_check($out_vcf);
    $self->_set_params_and_switches_from_args(%params);
    
    # *** should check the metrics file output for the last line which says how
    # many lines of output there should be excluding header:
    #Visited bases                                3001
    #Callable bases                               2996
    #Confidently called bases                     74
    #% callable bases of all loci                 99.833
    #% confidently called bases of all loci       2.466
    #% confidently called bases of callable loci  2.466
    #Actual calls made                            15
    
    return $self->run(@file_args);
}

=head2 variant_annotator

 Title   : variant_annotator
 Usage   : $wrapper->set_b('variant,VCF,variants.vcf');
           $wrapper->set_annotations('HaplotypeScore', 'SB', 'QD');
           $wrapper->variant_annotator('the.bam.list', 'out.vcf');
 Function: Annotates VCF calls.
 Returns : n/a
 Args    : path to input bam or list of bams (undef if not needed), path to
           output vcf file.
           Optionally, supply R or DBSNP options (as a hash), as understood by
           GATK, along with the other options like useAllAnnotations etc (1000
           genomes defaults exist).
           Before calling this, you should use set_b() to set the input vcf
           to the VCF you want to annotate

=cut

sub variant_annotator {
    my ($self, $input, $out_vcf, @params) = @_;
    
    # java -jar GenomeAnalysisTK.jar \
    #    -T VariantAnnotator \
    #    -R ref.fasta \
    #    -o output.vcf \
    #    -B variant,VCF,calls.vcf \
    #    -G Standard \      [use all standard annotations; don't run with '-all']
    #    -BTI variant \      [speeds up the runtime]
    #    -D dbsnp.rod 

    $self->switches([qw(useAllAnnotations)]);
    $self->params([qw(R DBSNP T L G group
                      rodToIntervalTrackName)]);
    
    my $bs = $self->get_b();
    my $ans = $self->get_annotations();
    my $I = $input ? " -I $input" : '';
    my @file_args = (" $bs $ans$I -o $out_vcf");
    
    my %params = (G => 'Standard', rodToIntervalTrackName => 'variant', @params);
    $params{T} = 'VariantAnnotator';
    $self->_handle_common_params(\%params);
    
    $self->register_output_file_to_check($out_vcf);
    $self->_set_params_and_switches_from_args(%params);
    
    return $self->run(@file_args);
}

=head2 variant_filtration

 Title   : variant_filtration
 Usage   : $wrapper->set_b('variant,VCF,snps.vcf',
                           'mask,Bed,indels.mask.bed');
           $wrapper->set_filters(
                     filters => {'filter_name' => 'filter expression', ... },
                     g_filters => {'filter_name' => 'filter expression', ... });
           $wrapper->variant_filtration('out.vcf');
 Function: Filters SNPs generated by unified_genotyper() to remove dodgy calls.
 Returns : n/a
 Args    : path to output vcf file.
           Optionally, supply R or DBSNP options (as a hash), as understood by
           GATK, along with the other options like clusterWindowSize etc (1000
           genomes defaults exist).
           Before calling this, you should use set_b() to set the input
           snps.vcf as produced by unified_genotyper(), and an indel mask as
           produced by the output of makeIndelMask.py on the raw bed output of
           indel_genotyper().
           You can also call set_filters() to set the filter expressions you
           want to use. 1000 genomes defaults exist.

=cut

sub variant_filtration {
    my ($self, $out_vcf, @params) = @_;
    
    #java -jar GenomeAnalysisTK.jar \
    #  -T VariantFiltration \
    #  -R resources/Homo_sapiens_assembly18.fasta \
    #  -o snps.filtered.vcf \
    #  -B variant,VCF,snps.raw.vcf \
    #  -B mask,Bed,indels.mask.bed \
    #  --maskName InDel \
    #  --clusterWindowSize 10
    #  --filterExpression "MQ0 >= 4 && (MQ0 / (1.0 * DP)) > 0.1" \
    #  --filterName "HARD_TO_VALIDATE"
    
    $self->switches([qw(missingValuesInExpressionsShouldEvaluateAsFailing)]);
    $self->params([qw(R DBSNP T clusterSize clusterWindowSize maskExtension maskName)]);
    
    my $bs = $self->get_b();
    my $filters = $self->get_filters();
    unless ($filters) {
        $self->set_filters(filters => { HARD_TO_VALIDATE => "MQ0 >= 4 && (MQ0 / (1.0 * DP)) > 0.1" });
        $filters = $self->get_filters();
    }
    my @file_args = (" $bs $filters -o $out_vcf");
    
    my %params = (maskName => 'InDel', clusterWindowSize => 10, @params);
    $params{T} = 'VariantFiltration';
    $self->_handle_common_params(\%params);
    
    $self->register_output_file_to_check($out_vcf);
    $self->_set_params_and_switches_from_args(%params);
    
    return $self->run(@file_args);
}

=head2 variant_recalibrator

 Title   : variant_recalibrator
 Usage   : $wrapper->set_b('input,VCF,snps.filtered.vcf');
           $wrapper->variant_recalibrator('in.cluster', 'out.vcf');
 Function: Recalibrates variant calls.
 Returns : n/a
 Args    : path to output of generate_variant_clusters(), path to output vcf
           file (which is also a prefix of other files generated, including
           .dat.tranches file needed by apply_variant_cuts()).
           Optionally, supply R or DBSNP options (as a hash), as understood by
           GATK, along with the other options like numGaussians etc (1000
           genomes defaults exist).
           Before calling this, you should use set_b() to set the input
           snps.filtered.vcf as produced by variant_filtration(), and the
           annotations to look at with set_annotations().

=cut

sub variant_recalibrator {
    my ($self, $out_recal_file, $out_tranches_file, @params) = @_;
    
    #java -Xmx4g -jar GenomeAnalysisTK.jar \
    #   -R /broad/1KG/reference/human_b36_both.fasta \
    #   -B input,VCF,snps.filtered.vcf \
    #   -B input2,VCF,another.snps.filtered.vcf \
    #   --DBSNP resources/dbsnp_129_hg18.rod \
    #   -l INFO \
    #   -clusterFile output.cluster \
    #   -o optimizer_output \
    #   -tranchesFile path/to/output.dat.tranches \
    #   -reportDatFile path/to/output.dat \
    #   --target_titv 2.07 \
    #   -resources R/ \
    #   --ignore_filter HARD_TO_VALIDATE \
    #   -T VariantRecalibrator
    
    my $tranches = ' -tranche 10 -tranche 5 -tranche 1 -tranche 0.1';

    # if it's there, check the tranches file for how many present and adjust -tranches
    # options accordingly.  This is a workaround for the SNP pipeline to work: the 
    # tranches file only has a line per tranche when SNPs fall into that tranche, so
    # pipeline breaks on empty tranches.  It's a known bug in GATK, but for now
    # this gets around it, since the pipeline will resubmit the job, pick up the tranches
    # file (with some lines missing) and only use the tranches that contain SNPs
    # if (-e "$out_vcf.dat.tranches") {
    #     open my $fh, "$out_vcf.dat.tranches" or $self->throw("$out_vcf.dat.tranches: $!");
    #     my @lines = <$fh>;
    #     close $fh;
    #     $tranches = "";
    #     foreach (@lines) {
    #         chomp;
    #         if (m/,FDRtranche.*to(.*)$/) {
    #             $tranches .= " -tranche $1";
    #         }
    #     }
    #     if ($tranches eq "") {
    #          $tranches = '-tranche 10 -tranche 5 -tranche 1 -tranche 0.1';
    #     }
    # }

    $self->switches([]);
    $self->params([qw(R DBSNP T l mode maxGaussians maxIterations numKMeans stdThreshold
                        qualThreshold shrinkage dirichlet priorCounts percentBadVariants
                        minNumBadVariants target_titv ignore_filter path_to_resources
                        ts_filter_level)]);
    
    my $bs = $self->get_b();
    my $ans = $self->get_annotations();
    my @file_args = (" $bs -tranchesFile $out_tranches_file -recal_file $out_recal_file",
                     '-resources '.File::Spec->catdir($ENV{GATK}, 'resources'),
                     '-Rscript Rscript', " -rscriptFile XXXXXX",
                     $tranches);
    
    my %params = (target_titv => 2.07, ignore_filter => 'HARD_TO_VALIDATE', l => 'INFO', @params);
    $params{T} = 'VariantRecalibrator';
    $self->_handle_common_params(\%params);
    
    $self->register_output_file_to_check($out_recal_file);
    $self->register_output_file_to_check($out_tranches_file);
    $self->_set_params_and_switches_from_args(%params);
    
    return $self->run(@file_args);
}

=head2 apply_recalibration

 Title   : apply_recalibration
 Usage   : $wrapper->set_b('input,VCF,recalibrator.output.vcf');
           $wrapper->apply_recalibration('recal.dat.tranches', 'out.vcf');
 Function: Filteres recalibrated calls.
 Returns : n/a
 Args    : path to the .tranches output of variant_recalibrator() (the second
           arg you supplied to that suffixed with .dat.traches), path to output
           vcf file. Optionally, supply R or DBSNP options (as a hash), as
           understood by GATK, along with the other options like
           fdr_filter_level etc (1000 genomes defaults exist). Before calling
           this, you should use set_b() to set the input recalibrator.output.vcf
           as produced by variant_recalibrator().

=cut

sub apply_recalibration {
    my ($self, $in_recal_file, $in_tranches, $out_vcf, @params) = @_;
    
    #java -Xmx6g -jar GenomeAnalysisTK.jar \
    #   -R /broad/1KG/reference/human_b36_both.fasta \
    #   -B input,VCF,recalibrator_output.vcf \
    #   --DBSNP resources/dbsnp_129_b36.rod \
    #   -l INFO \
    #   --ts_filter_level 10.0 \
    #   -tranchesFile output.cluster.dat.tranches \
    #   -o recalibrator_output.filtered.vcf \
    #   -T ApplyRecalibration
    
    $self->switches([]);
    $self->params([qw(R DBSNP T l ts_filter_level mode ignore_filter)]);
    
    my $bs = $self->get_b();
    my @file_args = (" $bs -recalFile $in_recal_file -tranchesFile $in_tranches -o $out_vcf");
    
    my %params = (ts_filter_level => '0.1', ignore_filter => 'HARD_TO_VALIDATE', l => 'INFO', @params);
    $params{T} = 'ApplyRecalibration';
    $self->_handle_common_params(\%params);
    
    $self->register_output_file_to_check($out_vcf);
    $self->_set_params_and_switches_from_args(%params);
    
    return $self->run(@file_args);
}

=head2 analyze_covariates

 Title   : analyze_covariates
 Usage   : $wrapper->analyze_covariates('recal_data.csv', 'output_dir');
 Function: Generates pdf graphs in an output dir.
 Returns : n/a
 Args    : path to input .csv file as generated by count_covariates(), output
           directory

=cut

sub analyze_covariates {
    my ($self, $csv, $out_dir) = @_;
    
    # java -Xmx4g -jar AnalyzeCovariates.jar \
    #   -recalFile /path/to/recal.table.csv  \
    #   -outputDir /path/to/output_dir/  \
    #   -resources resources/  \
    #   -ignoreQ 5
    
    # -Rscript   <path> Path to your implementation of Rscript. Default value: /broad/tools/apps/R-2.6.0/bin/Rscript
    # -resources <path>	Path to resources folder holding the Sting R scripts. Default value: R/
    
    $self->switches([]);
    $self->params([]);
    
    my $outdir = abs_path($out_dir);
    mkdir($outdir);
    
    my @file_args = (" -recalFile $csv",
                     " -outputDir $outdir",
                     " -ignoreQ 5",
                     " -resources ".File::Spec->catdir($ENV{GATK}, 'resources'),
                     " -Rscript Rscript");
    
    my $orig_exe = $self->exe;
    my $a_exe = $orig_exe;
    $a_exe =~ s/GenomeAnalysisTK\.jar/AnalyzeCovariates.jar/;
    $self->exe($a_exe);
    
    my @return = $self->run(@file_args);
    
    $self->exe($orig_exe);
    return @return;
}

=head2 combine_variants

 Title   : variant_eval
 Usage   : $wrapper->combine_variants('merged.vcf',
                                      ['in1.vcf', 'in2.vcf'],
                                      'UNION');
 Function: Merges VCF files together.
 Returns : n/a
 Args    : path to output file, array ref of input files (in priority order),
           the type of merge (the string 'UNION' or 'INTERSECT').
           Optionally, supply the R options (as a hash), as understood by
           GATK, along with the other options like genotypemergeoption (1000
           genomes defaults exist).

=cut

sub combine_variants {
    my ($self, $out_file, $in_files, $type, @params) = @_;
    
    #java -jar GenomeAnalysisTK.jar \
    #    -R ref.fasta \
    #    -T CombineVariants \
    #    -variantMergeOptions UNION \
    #    -B foo,VCF,foo.vcf \
    #    -B bar,VCF,bar.vcf \
    #    -priority foo,bar \
    #    -o merged.vcf
    
    $self->switches([qw(printComplexMerges filteredAreUncalled)]);
    $self->params([qw(R T genotypemergeoption)]);
    
    my @priority;
    my @bs;
    my %prev_names;
    my $ui = 0;
    foreach my $in (@{$in_files}) {
        my $base = basename($in);
        my $name = $base;
        $name =~ s/\.vcf.*//g;
        if ($name =~ /qcall/i) {
            $name = 'QCALL';
        }
        elsif ($name =~ /gatk/i) {
            $name = 'UG';
        }
        elsif ($name =~ /mpileup/i) {
            $name = 'MPILEUP';
        }
        else {
            $self->throw("Only qcall, mpileup and gatk are currently supported, by having those strings in the input vcf filenames");
        }
        
        if (exists $prev_names{$name}) {
            $name .= ++$ui;
        }
        $prev_names{$name} = 1;
        
        push(@bs, $name.',VCF,'.$in);
        push(@priority, $name);
    }
    my $priority = join(',', @priority);
    $self->set_b(@bs);
    my $bs = $self->get_b();
    
    my @file_args = (" $bs -variantMergeOptions $type -priority $priority -o $out_file");
    
    my %params = (@params);
    $params{T} = 'CombineVariants';
    $self->_handle_common_params(\%params);
    
    $self->register_output_file_to_check($out_file);
    $self->_set_params_and_switches_from_args(%params);
    
    return $self->run(@file_args);
}

=head2 variant_eval

 Title   : variant_eval
 Usage   : $wrapper->set_b('eval,VCF,snps.recal.vcf', 'comp,VCF,other.vcf');
           $wrapper->set_samples('NAxxxx', 'NAyyyyy');
           $wrapper->set_selects('name' => 'exp');
           $wrapper->set_eval_modules('', '');
           $wrapper->variant_eval('out.R');
 Function: Get lots of summary stats about a vcf in comparison to another.
 Returns : n/a
 Args    : path to output file.
           Optionally, supply R or DBSNP options (as a hash), as understood by
           GATK, along with the other options like reportType etc (1000
           genomes defaults exist).
           Before calling this, you should use set_b() to set the input
           vcf as produced by variant_filtration() or similar, and the
           selects to use with set_selects().

=cut

sub variant_eval {
    my ($self, $out_file, @params) = @_;
    
    #java -Xmx2048m -jar GenomeAnalysisTK.jar \
    #   -T VariantEval -R human_b36_both.fasta \
    #   -l INFO \
    #   -B eval,VCF,NA12878.vcf \
    #   -B comp,VCF,NA12891.vcf \
    #   -D /humgen/gsa-scr1/GATK_Data/dbsnp_129_b36.rod \
    #   -E DbSNPPercentage
    #   -o out.txt
    
    $self->switches([qw(indelCalls useNoModules)]);
    $self->params([qw(R DBSNP T l family_structure
                      MendelianViolationQualThreshold InterestingSitesVCF
                      minPhredConfidenceScore minPhredConfidenceScoreForComp
                      rsID maxRsIDBuild reportType reportLocation nSamples)]);
    
    # evalModule (default '') known_names
    
    my $bs = $self->get_b();
    my $selects = $self->get_selects();
    my $samples = $self->get_samples();
    my @file_args = (" $bs $selects -reportLocation $out_file");
    
    my %params = (reportType => 'R', l => 'INFO', @params);
    $params{T} = 'VariantEval';
    $self->_handle_common_params(\%params);
    
    $self->register_output_file_to_check($out_file);
    $self->_set_params_and_switches_from_args(%params);
    
    return $self->run(@file_args);
}

=head2 produce_beagle_input

 Title   : produce_beagle_input
 Usage   : $wrapper->set_b('variant,VCF,snps.recal.vcf');
           $wrapper->produce_beagle_input($file_for_beagle)
 Function: Produce beagle input based on a GATK VCF file (that may have been
           filtered, recalibrated etc.).
 Returns : n/a
 Args    : output file for giving to beagle
           Optionally, supply R option (as a hash).
           Before calling this, you should use set_b() to set the input
           snps.filtered.vcf as produced by variant_filtration() or similar.

=cut

sub produce_beagle_input {
    my ($self, $out_file, @params) = @_;
    
    #java -Xmx4000m -jar dist/GenomeAnalysisTK.jar -L 20 \
    #   -R reffile.fasta -T ProduceBeagleInput \
    #   -B:variant,VCF path_to_input_vcf/inputvcf.vcf
    #   -beagle path_to_beagle_output/beagle_output
    
    $self->switches([]);
    $self->params([qw(R T genotypes_file inserted_nocall_rate)]);
    
    my $bs = $self->get_b();
    my @file_args = (" $bs -o $out_file");
    
    my %params = (@params);
    $params{T} = 'ProduceBeagleInput';
    $self->_handle_common_params(\%params);
    
    $self->register_output_file_to_check($out_file);
    $self->_set_params_and_switches_from_args(%params);
    
    return $self->run(@file_args);
}

=head2 beagle_output_to_vcf

 Title   : beagle_output_to_vcf
 Usage   : $wrapper->set_b('variant,VCF,snps.recal.vcf',
                           'beagleR2,BEAGLE,beagle_out.r2',
                           'beaglePhased,BEAGLE,beagle_out.phased',
                           'beagleProbs,BEAGLE,beagle_out.gprobs');
           $wrapper->beagle_output_to_vcf($out_vcf)
 Function: Update an original GATK VCF with the beagle data from a beagle run
           on the output of produce_beagle_input() on that same VCF.
 Returns : n/a
 Args    : output vcf
           Optionally, supply R option (as a hash).
           Before calling this, you should use set_b() to set the input
           snps.filtered.vcf as produced by variant_filtration() or similar, and
           to provide the 3 beagle output files.

=cut

sub beagle_output_to_vcf {
    my ($self, $out_file, @params) = @_;
    
    #java -Xmx4000m -jar dist/GenomeAnalysisTK.jar \
    #   -R reffile.fasta -T BeagleOutputToVCF \
    #   -B:variant,VCF input_vcf.vcf \
    #   -B:beagleR2,BEAGLE /myrun.beagle_output.r2 \
    #   -B:beaglePhased,BEAGLE /myrun.beagle_output.phased \
    #   -B:beagleProbs,BEAGLE /myrun.beagle_output.gprobs \ 
    #   -o output_vcf.vcf 
    
    $self->switches([]);
    $self->params([qw(R T nocall_threshold)]);
    
    my $bs = $self->get_b();
    my @file_args = (" $bs -o $out_file");
    
    my %params = (@params);
    $params{T} = 'BeagleOutputToVCF';
    $self->_handle_common_params(\%params);
    
    $self->register_output_file_to_check($out_file);
    $self->_set_params_and_switches_from_args(%params);
    
    return $self->run(@file_args);
}

=head2 set_selects

 Title   : set_selects
 Usage   : $wrapper->set_selects('selectName' => 'select', ...);
 Function: Set which selects to use during variant_eval().
 Returns : hash (keys are select_names, values are select_exps)
 Args    : hash with keys as select_names (setting -selectName) and values as
           selects (setting -select)

=cut

sub set_selects {
    my $self = shift;
    if (@_) {
        my %hash = @_;
        $self->{the_selects} = \%hash;
    }
    return %{$self->{the_selects} || {}};
}

=head2 get_selects

 Title   : get_selects
 Usage   : my $selects_string = $wrapper->get_selects();
 Function: Get the command line options defining which selects to use for
           variant_eval().
 Returns : string
 Args    : n/a

=cut

sub get_selects {
    my $self = shift;
    my %selects = $self->set_selects;
    
    my $args = '';
    while (my ($name, $exp) = each %selects) {
        $args .= qq[ -select '$exp' -selectName '$name'];
    }
    
    return $args;
}

=head2 set_samples

 Title   : set_samples
 Usage   : $wrapper->set_samples('NAyyyy', 'NAxxxx', ...);
 Function: Set which samples to consuder during variant_eval().
 Returns : list currently set
 Args    : list of sample name strings

=cut

sub set_samples {
    my $self = shift;
    if (@_) {
        $self->{samples} = [@_];
    }
    return @{$self->{samples} || []};
}

=head2 get_samples

 Title   : get_samples
 Usage   : my $samples_string = $wrapper->get_samples();
 Function: Get the command line options defining which samples to use during
           variant_eval(). If none have been set with set_samples() returns
           empty string.
 Returns : string
 Args    : n/a

=cut

sub get_samples {
    my $self = shift;
    my @samples = $self->set_samples;
    
    my $args = '';
    foreach my $s (@samples) {
        $args .= "--samples $s ";
    }
    return $args;
}

=head2 recalibrate

 Title   : recalibrate
 Usage   : $wrapper->recalibrate('input.bam', 'out.bam');
 Function: Easy-to-use recalibration. Just runs count_covariates() followed by
           table_recalibration(). Also ensures the output bam isn't truncated.
           Won't attempt to recalibrate if out.bam already exists.
 Returns : n/a
 Args    : path to input .bam file, path to output file. Optionally, supply R,
           DBSNP, useOriginalQualities or l options (as a hash), as understood by
           GATK. useOriginalQualities is on by default. -B and -cov should be set
           with the set_b() and set_covs() methods beforehand.

=cut

sub recalibrate {
    my ($self, $in_bam, $out_bam, @params) = @_;
    
    my $orig_run_method = $self->run_method;
    $self->run_method('system');
    
    # count_covariates
    my $csv = $in_bam.".recal_data.csv";
    my $csv_tmp = $in_bam.".recal_data.tmp.csv";
    unless (-s $csv) {
        $self->count_covariates($in_bam, $csv_tmp, @params);
        $self->throw("failed during the count_covariates step, giving up for now") unless $self->run_status >= 1;
        # *** need to be able to test csv_tmp for truncation... but just assume
        #     that if we're still alive, all is well
        move($csv_tmp, $csv) || $self->throw("Could not move $csv_tmp to $csv");
    }
    
    # table_recalibration
    unless (-s $out_bam) {
        my $tmp_out = $out_bam;
        $tmp_out =~ s/\.bam$/.tmp.bam/;
        $self->table_recalibration($in_bam, $csv, $tmp_out, @params);
        $self->throw("failed during the table_recalibration step, giving up for now") unless $self->run_status >= 1;
        
        # find out how many lines are in the input bam file
        my $st = VertRes::Wrapper::samtools->new(quiet => 1);
        $st->run_method('open');
        my $bam_count = 0;
        my $fh = $st->view($in_bam);
        while (<$fh>) {
            $bam_count++;
        }
        close($fh);
        
        # find out how many lines are in the recalibrated bam file
        my $recal_count = 0;
        $fh = $st->view($tmp_out);
        while (<$fh>) {
            $recal_count++;
        }
        close($fh);
        
        # check for truncation
        if ($recal_count >= $bam_count) {
            move($tmp_out, $out_bam) || $self->throw("Failed to move $tmp_out to $out_bam: $!");
            $self->_set_run_status(2);
        }
        else {
            $self->warn("$tmp_out is bad (only $recal_count lines vs $bam_count), will unlink it");
            $self->_set_run_status(-1);
            unlink($tmp_out);
        }
    }
    
    $self->run_method($orig_run_method);
    return;
}

sub _pre_run {
    my $self = shift;
    $self->_set_params_string(mixed_dash => 1);
    
    my $input_bam = $_[0];
    if ($input_bam =~ /\.bam$/) {
        $input_bam =~ s/.* -I\s+//;
        my $bai_file = $input_bam.'.bai';
        
        # check that the bai is older than the bam
        if (-e $bai_file && $self->_is_older($input_bam, $bai_file)) {
            unlink($bai_file);
        }
        
        # create a bai file if necessary
        unless (-e $bai_file) {
            my $sam_wrapper = VertRes::Wrapper::samtools->new(verbose => $self->verbose,
                                                              run_method => 'system');
            $sam_wrapper->index($input_bam, $bai_file);
            
            unless ($sam_wrapper->run_status >= 1) {
                if ($sam_wrapper->run_status == -1) {
                    $self->warn("Failed to create index file '$bai_file', will try one more time and then give up...");
                    sleep(5);
                    $sam_wrapper->index($input_bam, $bai_file);
                }
                unless ($sam_wrapper->run_status >= 1) {
                    $self->throw("Failed to create index file '$bai_file', giving up!");
                }
            }
        }
    }
    
    $self->register_for_unlinking('GATK_Error.log');
    
    return @_;
}

1;
