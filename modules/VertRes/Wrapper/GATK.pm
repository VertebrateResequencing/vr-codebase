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

Ostensibly a wrapper for Broad's GenomeAnalysisToolKit, this is primarily
focused on using it to recalibrate the quality values in bam files. See:
http://www.broadinstitute.org/gsa/wiki/index.php/Quality_scores_recalibration

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
           vcfs      => [] (as per set_vcfs())
           build     => NCBI36|NCBI37|NCBIM37 (default NCBI36: sets defaults for
                        reference, dbsnp, covs and vcfs as appropriate for the
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
    my $build = delete $self->{build} || 'NCBI36';
    my ($default_ref, $default_dbsnp, $default_covs, $default_vcfs);
    if ($build eq 'NCBI36') {
        $default_ref = File::Spec->catfile($ENV{GATK_RESOURCES}, 'human_b36_both.fasta');
        $default_dbsnp = File::Spec->catfile($ENV{GATK_RESOURCES}, 'dbsnp_129_b36.rod');
    }
    elsif ($build eq 'NCBI37') {
        $default_ref = File::Spec->catfile($ENV{GATK_RESOURCES}, 'human_g1k_v37.fasta');
        $default_dbsnp = File::Spec->catfile($ENV{GATK_RESOURCES}, 'dbsnp_130_b37.rod');
        # no default vcfs for now, since issues with not being in NCBI37 coords?
        #$default_vcfs = ['pilot1_CEU,VCF,'.File::Spec->catfile($ENV{GATK_RESOURCES}, 'vcfs', 'CEU.2and3_way.vcf'),
        #                 'pilot1_YRI,VCF,'.File::Spec->catfile($ENV{GATK_RESOURCES}, 'vcfs', 'YRI.2and3_way.vcf')];
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
    $self->set_vcfs(@{delete $self->{vcfs} || $default_vcfs || []});
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
           with the set_vcfs() and set_covs() methods beforehand.

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
    my $vcfs = $self->get_vcfs();
    my @file_args = (" -I $in_bam", " $covs $vcfs -recalFile $recal_file");
    
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

=head2 set_vcfs

 Title   : set_vcfs
 Usage   : $wrapper->set_vcfs('pilot1,VCF,file1.vcf', 'pilot2,VCF,file2.vcf');
 Function: Set which vcf files to use (for setting B option).
 Returns : list currently set
 Args    : list of vcf strings (name,type,filename)

=cut

sub set_vcfs {
    my $self = shift;
    if (@_) {
        $self->{vcfs} = [@_];
    }
    return @{$self->{vcfs} || []};
}

=head2 get_vcfs

 Title   : get_vcfs
 Usage   : my $covs_string = $wrapper->get_vcfs();
 Function: Get the command line options defining which vcf files will be used
           (one or more -B options). If none have been set with set_vcfs,
           returns empty string.
 Returns : string
 Args    : n/a

=cut

sub get_vcfs {
    my $self = shift;
    my @vcfs = $self->set_vcfs;
    
    my $args = '';
    foreach my $vcf (@vcfs) {
        $args .= "-B $vcf ";
    }
    return $args;
}

=head2  table_recalibration

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

=head2   analyze_covariates

 Title   : analyze_covariates
 Usage   : $wrapper->analyze_covariates('recal_data.csv', 'output.pdf');
 Function: Generates pdf graphs .
 Returns : boolean (true on success)
 Args    : path to input .csv file as generated by count_covariates()

=cut

sub analyze_covariates {
    my ($self, $csv, $out_pdf) = @_;
    
    # java -Xmx4g -jar AnalyzeCovariates.jar \
    #   -recalFile /path/to/recal.table.csv  \
    #   -outputDir /path/to/output_dir/  \
    #   -resources resources/  \
    #   -ignoreQ 5
    
    # -Rscript   <path> Path to your implementation of Rscript. Default value: /broad/tools/apps/R-2.6.0/bin/Rscript
    # -resources <path>	Path to resources folder holding the Sting R scripts. Default value: R/
    
    $self->switches([]);
    $self->params([]);
    
    my $outdir = abs_path("${out_pdf}_working");
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
           with the set_vcfs() and set_covs() methods beforehand.

=cut

sub recalibrate {
    my ($self, $in_bam, $out_bam, @params) = @_;
    
    my $orig_run_method = $self->run_method;
    $self->run_method('system');
    
    # count_covariates
    my $csv = $in_bam.".recal_data.csv";
    unless (-s $csv) {
        $self->count_covariates($in_bam, $in_bam, @params);
        $self->throw("failed during the count_covariates step, giving up for now") unless $self->run_status >= 1;
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
