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

our $DEFAULT_GATK_JAR = File::Spec->catfile($ENV{GATK}, 'GenomeAnalysisTK.jar');
our $DEFAULT_LOGLEVEL = 'ERROR';
our $DEFAULT_PLATFORM = 'ILLUMINA';

=head2 new

 Title   : new
 Usage   : my $wrapper = VertRes::Wrapper::GATK->new();
 Function: Create a VertRes::Wrapper::GATK object.
 Returns : VertRes::Wrapper::GATK object
 Args    : quiet => boolean
           exe   => string (full path to GenomeAnalysisTK.jar; a TEAM145 default
                            exists)
           java_memory => int (the amount of memory in MB to give java; default
                               6000)
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

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(exe => $DEFAULT_GATK_JAR, @args);
    
    my $java_mem = delete $self->{java_memory} || 6000;
    $self->exe("java -Xmx${java_mem}m -jar ".$self->exe);
    
    # our bsub jobs will get killed if we don't select high-mem machines
    $self->bsub_options(M => ($java_mem * 1000), R => "'select[mem>$java_mem] rusage[mem=$java_mem]'");
    
    # default settings
    my $build = delete $self->{build} || 'NCBI37';
    my ($default_ref, $default_dbsnp, $default_covs, $default_bs);
    if ($build eq 'NCBI36') {
        $default_ref = File::Spec->catfile($ENV{GATK_RESOURCES}, 'human_b36_both.fasta');
        $default_dbsnp = File::Spec->catfile($ENV{GATK_RESOURCES}, 'dbsnp_129_b36.rod');
    }
    elsif ($build eq 'NCBI37') {
        $default_ref = File::Spec->catfile($ENV{GATK_RESOURCES}, 'human_g1k_v37.fasta');
        $default_dbsnp = File::Spec->catfile($ENV{GATK_RESOURCES}, 'dbsnp_130_b37.rod');
        # no default bs for now, since issues with not being in NCBI37 coords?
        #$default_bs = ['pilot1_CEU,VCF,'.File::Spec->catfile($ENV{GATK_RESOURCES}, 'vcfs', 'CEU.2and3_way.vcf'),
        #               'pilot1_YRI,VCF,'.File::Spec->catfile($ENV{GATK_RESOURCES}, 'vcfs', 'YRI.2and3_way.vcf')];
    }
    elsif ($build eq 'NCBIM37') {
        $self->throw("mouse .rod not available; suggest not attempting recalibration on mouse at the moment...");
    }
    else {
        $self->throw("bad build option '$build'");
    }
    
    $self->{_default_R} = delete $self->{reference} || $default_ref;
    $self->{_default_DBSNP} = delete $self->{dbsnp} || $default_dbsnp;
    $self->set_covs(@{delete $self->{covs} || $default_covs || []});
    $self->set_b(@{delete $self->{bs} || $default_bs || []});
    $self->{_default_loglevel} = delete $self->{log_level} || $DEFAULT_LOGLEVEL;
    $self->{_default_platform} = delete $self->{default_platform} || $DEFAULT_PLATFORM;
    
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

    $self->switches([qw(quiet_output_mode useOriginalQualities)]);
    $self->params([qw(R DBSNP l T max_reads_at_locus default_platform)]);
    
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
    $params{quiet_output_mode} = $self->quiet();
    unless (defined $params{useOriginalQualities}) {
        $params{useOriginalQualities} = 1;
    }
    $self->_handle_common_params(\%params);
    $params{max_reads_at_locus} ||= 50000; # stop it using tons of memory in repeat regions
    
    $self->throw("Non-existant rod file '$params{DBSNP}'") unless -s $params{DBSNP};
    
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
           (one or more -an options). If none have been set with
           set_annotations(), defaults to 'HaplotypeScore', 'SB', 'QD'.
 Returns : string
 Args    : n/a

=cut

sub get_annotations {
    my $self = shift;
    my @ans = $self->set_annotations;
    unless (@ans) {
        @ans = ('HaplotypeScore', 'SB', 'QD');
    }
    
    my $args = '';
    foreach my $an (@ans) {
        $args .= "-an $an ";
    }
    
    return $args;
}

=head2 set_b

 Title   : set_b
 Usage   : $wrapper->set_b('pilot1,VCF,file1.vcf', 'pilot2,VCF,file2.vcf');
 Function: Set which vcf files to use (for setting B option).
 Returns : list currently set
 Args    : list of vcf strings (name,type,filename)

=cut

sub set_b {
    my $self = shift;
    if (@_) {
        $self->{bs} = [@_];
    }
    return @{$self->{bs} || []};
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
        $args .= "-B $b ";
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
    #   -outputBAM my_reads.recal.bam \
    #   -recalFile my_reads.recal_data.csv
    
    $self->switches([qw(quiet_output_mode useOriginalQualities)]);
    $self->params([qw(R l T default_platform)]);
    
    my @file_args = (" -I $in_bam", " -recalFile $csv", " --output_bam $out_bam");
    
    my %params = @params;
    $params{T} = 'TableRecalibration';
    $params{quiet_output_mode} = $self->quiet();
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
 Usage   : $wrapper->realignment_targets('in.bam', 'out.intervals');
 Function: Finds target intervals in a quality-recalibrated bam that could be
           realigned with indel_realigner().
 Returns : n/a
 Args    : path to input .bam file, path to output intervals.
           Optionally, supply R or DBSNP options (as a hash), as understood by
           GATK.

=cut

sub realignment_targets {
    my ($self, $in_bam, $out_intervals, @params) = @_;
    
    # java -jar GenomeAnalysisTK.jar \
    #   -T RealignerTargetCreator \
    #   -I recalibrated.bam \
    #   -R resources/Homo_sapiens_assembly18.fasta \
    #   -o forRealigner.intervals \
    #   -D resources/dbsnp_129_hg18.rod  (-D == DBSNP ??)
    
    $self->switches([qw(quiet_output_mode)]);
    $self->params([qw(R DBSNP T)]);
    
    my @file_args = (" -I $in_bam", " -o $out_intervals");
    
    my %params = @params;
    $params{T} = 'RealignerTargetCreator';
    $params{quiet_output_mode} = $self->quiet();
    $self->_handle_common_params(\%params);
    
    $self->register_output_file_to_check($out_intervals);
    $self->_set_params_and_switches_from_args(%params);
    
    return $self->run(@file_args);
}

=head2 indel_realigner

 Title   : indel_realigner
 Usage   : $wrapper->indel_realigner('in.bam', 'intervals', 'out.bam');
 Function: Does local realignment around indels in intervals as determined by
           realignment_targets(), generating a "cleaned" bam suitable for
           calling SNPs on.
 Returns : n/a
 Args    : path to input .bam file, path to output of realignment_targets(),
           path to output bam.
           Optionally, supply R or DBSNP options (as a hash), as understood by
           GATK.

=cut

sub indel_realigner {
    my ($self, $in_bam, $intervals_file, $out_bam, @params) = @_;
    
    # java -Djava.io.tmpdir=/path/to/tmpdir -jar GenomeAnalysisTK.jar \
    #   -I recalibrated.bam \
    #   -R resources/Homo_sapiens_assembly18.fasta \
    #   -T IndelRealigner \
    #   -targetIntervals forRealigner.intervals \
    #   --output cleaned.bam \
    #   -D resources/dbsnp_129_hg18.rod
    
    $self->switches([qw(quiet_output_mode)]);
    $self->params([qw(R DBSNP T)]);
    
    my @file_args = (" -I $in_bam", " -targetIntervals $intervals_file", " --output $out_bam");
    
    my %params = @params;
    $params{T} = 'IndelRealigner';
    $params{quiet_output_mode} = $self->quiet();
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
    
    $self->switches([qw(quiet_output_mode verbose)]);
    $self->params([qw(R DBSNP T 1kg_format minCoverage minNormalCoverage
                      minFraction minConsensusFraction minIndelCount refseq
                      blacklistedLanes window_size maxNumberOfReads)]);
    
    my @file_args = (" -I $in_bam", " -O $out_raw_bed -o $out_detailed_bed");
    
    my %params = (verbose => 1, minIndelCount => 2, minFraction => 0.03,
                  minConsensusFraction => 0.6, maxNumberOfReads => 1000000,
                  @params);
    $params{T} = 'IndelGenotyperV2';
    $params{quiet_output_mode} = $self->quiet();
    $self->_handle_common_params(\%params);
    
    $self->register_output_file_to_check($out_raw_bed);
    $self->register_output_file_to_check($out_detailed_bed);
    $self->_set_params_and_switches_from_args(%params);
    
    return $self->run(@file_args);
}

=head2 unified_genotyper

 Title   : unified_genotyper
 Usage   : $wrapper->unified_genotyper('in.bam', 'out.vcf', 'out.beagle');
 Function: Call SNPs on a bam file (preferably one that has been output by
           indel_realigner()).
 Returns : n/a
 Args :    path to input .bam file (or bam fofn), paths to output files (vcf,
           and file suitable for input into Beagle.
           Optionally, supply R or DBSNP options (as a hash), as understood by
           GATK, along with the other options like confidence etc (1000
           genomes defaults exist).

=cut

sub unified_genotyper {
    my ($self, $in_bam, $out_vcf, $out_beagle, @params) = @_;
    
    # java -jar GenomeAnalysisTK.jar \
    #   -R resources/Homo_sapiens_assembly18.fasta \
    #   -T UnifiedGenotyper \
    #   -I cleaned.bam \
    #   -D resources/dbsnp_129_hg18.rod \
    #   -varout snps.raw.vcf \
    #   --standard_min_confidence_threshold_for_calling 10.0 \
    #   -beagle snps.beagle
    
    $self->switches([qw(quiet_output_mode genotype output_all_callable_bases
                        noSLOD)]);
    $self->params([qw(R DBSNP T confidence genotype_model base_model
                      heterozygosity max_reads_at_locus
                      standard_min_confidence_threshold_for_calling
                      standard_min_confidence_threshold_for_emitting
                      trigger_min_confidence_threshold_for_calling
                      trigger_min_confidence_threshold_for_emitting
                      assume_single_sample_reads platform
                      min_base_quality_score min_mapping_quality_score
                      max_mismatches_in_40bp_window use_reads_with_bad_mates
                      max_deletion_fraction cap_base_quality_by_mapping_quality
                      variant_output_format verbose_mode annotation group)]);
    
    my @file_args = (" -I $in_bam", " -varout $out_vcf -beagle $out_beagle");
    
    my %params = (standard_min_confidence_threshold_for_calling => 10.0, @params);
    $params{T} = 'UnifiedGenotyper';
    $params{quiet_output_mode} = $self->quiet();
    $params{max_reads_at_locus} ||= 50000; # stop it using tons of memory in repeat regions
    $self->_handle_common_params(\%params);
    
    $self->register_output_file_to_check($out_vcf);
    $self->register_output_file_to_check($out_beagle);
    $self->_set_params_and_switches_from_args(%params);
    
    return $self->run(@file_args);
}

=head2 variant_annotator

 Title   : variant_annotator
 Usage   : $wrapper->set_b('variant,VCF,variants.vcf');
           $wrapper->variant_annotator('the.bam.list', 'out.vcf');
 Function: Annotates VCF calls.
 Returns : n/a
 Args    : path to input bam or list of bams, path to output vcf file.
           Optionally, supply R or DBSNP options (as a hash), as understood by
           GATK, along with the other options like clusterWindowSize etc (1000
           genomes defaults exist).
           Before calling this, you should use set_b() to set the input vcf
           to the VCF you want to annotate

=cut

sub variant_annotator {
    my ($self, $input, $out_vcf, @params) = @_;
    
    # java -jar /path/to/dist/GenomeAnalysisTK.jar \
    #   -T VariantAnnotator \
    #   -l INFO
    #   --DBSNP resources/dbsnp_129_hg18.rod \
    #   -R /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta \
    #   -I /path/to/bam/file.bam.list \
    #   -o /path/to/output.vcf \
    #   -B variant,VCF,/path/to/input/variants.vcf \

    $self->switches([qw(quiet_output_mode)]);
    $self->params([qw(R DBSNP T L)]);
    
    my $bs = $self->get_b();
    my @file_args = (" $bs -I $input -o $out_vcf");
    
    my %params = (@params);
    $params{T} = 'VariantAnnotator';
    $params{quiet_output_mode} = $self->quiet();
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

    $self->switches([qw(quiet_output_mode)]);
    $self->params([qw(R DBSNP T clusterSize clusterWindowSize maskName)]);
    
    my $bs = $self->get_b();
    my $filters = $self->get_filters();
    # used for maq; variant recalibration will handle bwa alignments without
    # this
    #unless ($filters) {
    #    $self->set_filters(filters => { HARD_TO_VALIDATE => "(MQ0 / (1.0 * DP)) > 0.1" });
    #    $filters = $self->get_filters();
    #}
    my @file_args = (" $bs $filters -o $out_vcf");
    
    my %params = (maskName => 'InDel', clusterWindowSize => 10, @params);
    $params{T} = 'VariantFiltration';
    $params{quiet_output_mode} = $self->quiet();
    $self->_handle_common_params(\%params);
    
    $self->register_output_file_to_check($out_vcf);
    $self->_set_params_and_switches_from_args(%params);
    
    return $self->run(@file_args);
}

=head2 generate_variant_clusters

 Title   : generate_variant_clusters
 Usage   : $wrapper->set_b('input,VCF,snps.filtered.vcf');
           $wrapper->set_annotations('HaplotypeScore', 'SB', 'QD');
           $wrapper->generate_variant_clusters('cluster.out');
 Function: Clusters variants for later variant recalibration.
 Returns : n/a
 Args    : path to output cluster file.
           Optionally, supply R or DBSNP options (as a hash), as understood by
           GATK, along with the other options like numGaussians etc (1000
           genomes defaults exist).
           Before calling this, you should use set_b() to set the input
           snps.filtered.vcf as produced by variant_filtration(), and the
           annotations to look at with set_annotations().

=cut

sub generate_variant_clusters {
    my ($self, $out_cluster, @params) = @_;
    
    #java -Xmx4g -jar GenomeAnalysisTK.jar \
    #   -R /broad/1KG/reference/human_b36_both.fasta \
    #   -B input,VCF,snps.filtered.vcf \
    #   -B input2,VCF,another.snps.filtered.vcf \
    #   --DBSNP resources/dbsnp_129_hg18.rod \
    #   -l INFO \
    #   -nG 6 \  [--numGaussians]
    #   -nI 10 \ [--numIterations]
    #   -an HaplotypeScore -an SB -an QD \
    #   -clusterFile output.cluster \
    #   -resources R/ \
    #   -T GenerateVariantClusters
    
    $self->switches([qw(quiet_output_mode ignore_all_input_filters)]);
    $self->params([qw(R DBSNP T numGaussians numIterations ignore_filter
                      minVarInCluster weightKnowns weightHapMap
                      weight1000Genomes weightMQ1)]);
    
    my $bs = $self->get_b();
    my $ans = $self->get_annotations();
    my @file_args = (" $bs $ans -clusterFile $out_cluster",
                     '-resources '.File::Spec->catdir($ENV{GATK}, 'resources').'/',
                     '-Rscript Rscript');
    
    my %params = (numGaussians => 6, numIterations => 10, @params);
    $params{T} = 'GenerateVariantClusters';
    $params{quiet_output_mode} = $self->quiet();
    $self->_handle_common_params(\%params);
    
    $self->register_output_file_to_check($out_cluster);
    $self->_set_params_and_switches_from_args(%params);
    
    #*** after completion, remove $out_cluster.1 .. $out_cluster.$numIterations
    
    return $self->run(@file_args);
}

=head2 variant_recalibrator

 Title   : variant_recalibrator
 Usage   : $wrapper->set_b('input,VCF,snps.filtered.vcf');
           $wrapper->variant_recalibrator('in.cluster', 'out.vcf');
 Function: Recalibrates variant calls.
 Returns : n/a
 Args    : path to output of generate_variant_clusters(), path to output vcf
           file.
           Optionally, supply R or DBSNP options (as a hash), as understood by
           GATK, along with the other options like numGaussians etc (1000
           genomes defaults exist).
           Before calling this, you should use set_b() to set the input
           snps.filtered.vcf as produced by variant_filtration(), and the
           annotations to look at with set_annotations().

=cut

sub variant_recalibrator {
    my ($self, $in_cluster, $out_vcf, @params) = @_;
    
    #java -Xmx4g -jar GenomeAnalysisTK.jar \
    #   -R /broad/1KG/reference/human_b36_both.fasta \
    #   -B input,VCF,snps.filtered.vcf \
    #   -B input2,VCF,another.snps.filtered.vcf \
    #   --DBSNP resources/dbsnp_129_hg18.rod \
    #   -l INFO \
    #   -clusterFile output.cluster \
    #   -output optimizer_output \
    #   --target_titv 2.1 \
    #   -resources R/ \
    #   -T VariantRecalibrator
    
    $self->switches([qw(quiet_output_mode ignore_all_input_filters)]);
    $self->params([qw(R DBSNP T target_titv backOff desired_num_variants
                      ignore_filter known_prior novel_prior
                      quality_scale_factor)]);
    
    my $bs = $self->get_b();
    $out_vcf =~ s/\.vcf//; # it adds .vcf
    my @file_args = (" $bs -clusterFile $in_cluster -output $out_vcf",
                     '-resources '.File::Spec->catdir($ENV{GATK}, 'resources'),
                     '-Rscript Rscript');
    
    my %params = (target_titv => 2.1, @params);
    $params{T} = 'VariantRecalibrator';
    $params{quiet_output_mode} = $self->quiet();
    $self->_handle_common_params(\%params);
    
    $self->register_output_file_to_check($out_vcf.'.vcf');
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

=head2 variant_eval

 Title   : variant_eval
 Usage   : $wrapper->set_b('eval,VCF,snps.recal.vcf', 'comp,VCF,other.vcf');
           $wrapper->set_samples('NAxxxx', 'NAyyyyy');
           $wrapper->set_selects('name' => 'exp');
           $wrapper->set_eval_modules('', '');
           $wrapper->variant_eval('out.grepable');
 Function: Get lots of summary stats about a vcf in comparison to another.
 Returns : n/a
 Args    : path to output file.
           Optionally, supply R or DBSNP options (as a hash), as understood by
           GATK, along with the other options like numGaussians etc (1000
           genomes defaults exist).
           Before calling this, you should use set_b() to set the input
           snps.filtered.vcf as produced by variant_filtration(), and the
           annotations to look at with set_annotations().

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
    
    $self->switches([qw(quiet_output_mode indelCalls useNoModules)]);
    $self->params([qw(R DBSNP T family_structure
                      MendelianViolationQualThreshold InterestingSitesVCF
                      minPhredConfidenceScore minPhredConfidenceScoreForComp
                      rsID maxRsIDBuild reportType reportLocation nSamples)]);
    
    # evalModule (default '') known_names
    
    my $bs = $self->get_b();
    my $selects = $self->get_selects();
    my $samples = $self->get_samples();
    my @file_args = (" $bs $selects -o $out_file");
    
    my %params = (reportType => 'grep', @params);
    $params{T} = 'VariantEval';
    $params{quiet_output_mode} = $self->quiet();
    $self->_handle_common_params(\%params);
    
    $self->register_output_file_to_check($out_file);
    $self->_set_params_and_switches_from_args(%params);
    
    return $self->run(@file_args);
}

=head2 set_selects

 Title   : set_selects
 Usage   : $wrapper->set_selects('select_name' => 'select_exp', ...);
 Function: Set which selects to use during variant_eval().
 Returns : hash (keys are select_names, values are select_exps)
 Args    : hash with keys as select_names (setting --select_names) and values as
           select_exps (setting --select_exps)

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
        $args .= " --select_exps \"$exp\" --select_names \"$name\"";
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
        $input_bam =~ s/^ -I //;
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
