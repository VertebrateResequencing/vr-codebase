package Pathogens::QC::HetSNPCalculator;

use Moose;
use Utils;
use File::Spec;
use Pathogens::Exception;

#executables
has 'samtools' => ( is => 'ro', isa => 'Str', required => 1 );
has 'bcftools' => ( is => 'ro', isa => 'Str', required => 1 );

#bcftools filters
has 'min_rawReadDepth'  => ( is => 'ro', isa => 'Str', required => 1 );
has 'min_hqNonRefBases' => ( is => 'ro', isa => 'Str', required => 1 );
has 'rawReadDepth_hqNonRefBases_ratio' =>
  ( is => 'ro', isa => 'Str', required => 1 );
has 'min_qual' => ( is => 'ro', isa => 'Str', required => 1 );
has 'hqRefReads_hqAltReads_ratio' =>
  ( is => 'ro', isa => 'Str', required => 1 );

#Misc required info
has 'fa_ref'         => ( is => 'ro', isa => 'Str', required => 1 );
has 'lane_path'      => ( is => 'ro', isa => 'Str', required => 1 );
has 'lane'           => ( is => 'ro', isa => 'Str', required => 1 );
has 'sample_dir'     => ( is => 'ro', isa => 'Str', required => 1 );
has 'het_report'     => ( is => 'ro', isa => 'Str', required => 1 );

#File path builders
has 'full_path' =>
  ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_full_path' );
has 'bam_file' =>
  ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_bam_file' );
has 'temp_vcf' =>
  ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_temp_vcf' );
has 'snp_called_vcf' =>
  ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_snp_called_vcf' );
has 'all_snps_csv' =>
  ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_all_snps_csv' );
has 'filtered_snp_called_csv' => (
    is      => 'rw',
    isa     => 'Str',
    lazy    => 1,
    builder => 'build_filtered_snp_called_csv'
);
has 'het_report_path' =>
  ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_het_report_path' );

#Command string builders
has 'mpileup_command' =>
  ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_mpileup_command' );
has 'all_snps_command' =>
  ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_all_snps_command' );
has 'bcf_query_command' =>
  ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_bcf_query_command' );


sub build_full_path {

    my ($self) = @_;
    my $full_path = File::Spec->catfile( $self->lane_path, $self->sample_dir );
    return ($full_path);
}

sub build_bam_file {

    my ($self) = @_;
    my $path = File::Spec->catfile( $self->full_path, $self->{lane} . q(.bam) );
    return ($path);
}

sub build_temp_vcf {

    my ($self) = @_;
    my $path = File::Spec->catfile( $self->full_path,
        $self->{lane} . q(_temp_vcf.vcf) );
    return ($path);
}

sub build_snp_called_vcf {

    my ($self) = @_;
    my $path = File::Spec->catfile( $self->full_path,
        $self->{lane} . q(_snp_called.vcf.gz) );
    return ($path);
}

sub build_all_snps_csv {

    my ($self) = @_;
    my $path = File::Spec->catfile( $self->full_path,
        $self->{lane} . q(_all_snps_list.csv) );
    return ($path);
}

sub build_filtered_snp_called_csv {

    my ($self) = @_;
    my $path = File::Spec->catfile( $self->full_path,
        $self->{lane} . q(_filtered_snp_called_list.csv) );
    return ($path);
}


sub build_het_report_path {

    my ($self) = @_;
    my $het_report_path = File::Spec->catfile( $self->lane_path,
        $self->lane . q(_) . $self->het_report );
    return ($het_report_path);
}

sub build_mpileup_command {

    my ($self) = @_;

    return $self->samtools
       . ' mpileup --skip-indels -d 500 -t INFO/AD,INFO/ADF,INFO/ADR -C50 -u'
       . " -f " . $self->fa_ref
       . ' ' . $self->bam_file
       . ' > ' . $self->temp_vcf;
}


sub build_snp_call_command {

    my ($self) = @_;

    my $cmd = $self->bcftools;
    $cmd .= q( call -vm -O z );
    $cmd .= $self->temp_vcf;
    $cmd .= q( > );
    $cmd .= $self->snp_called_vcf;

    return ($cmd);
}

sub build_all_snps_command {

    my ($self) = @_;

    my $cmd = $self->{bcftools} . q( query -f);
    $cmd .= q{ "%CHROM %POS\n" -i};
    $cmd .= q{ "MIN(DP) >= } . $self->{min_rawReadDepth};
    $cmd .= q{ & MIN(DV) >= } . $self->{min_hqNonRefBases};
    $cmd .= q{ & MIN(DV/DP)>= } . $self->{rawReadDepth_hqNonRefBases_ratio};
    $cmd .= q{ & QUAL >= } . $self->{min_qual};
    $cmd .= q{ & (GT='0/0' | GT='1/1' | GT='0/1' | GT='1/2')" };
    $cmd .= $self->snp_called_vcf;
    $cmd .= q{ > } . $self->all_snps_csv;

    return ($cmd);
}

sub build_bcf_query_command {

    my ($self) = @_;

    my $cmd = $self->{bcftools} . q( query -f);
    $cmd .= q{ "%CHROM %POS\n" -i};
    $cmd .= q{ "MIN(DP) >= } . $self->{min_rawReadDepth};
    $cmd .= q{ & MIN(DV) >= } . $self->{min_hqNonRefBases};
    $cmd .= q{ & MIN(DV/DP)>= } . $self->{rawReadDepth_hqNonRefBases_ratio};
    $cmd .= q{ & QUAL >= } . $self->{min_qual};
    $cmd .= q{ & (GT='0/1' | GT='1/2')};
    $cmd .= q{ & ((DP4[0]+DP4[1])/(DP4[2]+DP4[3]) > }
      . $self->{hqRefReads_hqAltReads_ratio} . q{)" };
    $cmd .= $self->snp_called_vcf;
    $cmd .= q{ > } . $self->filtered_snp_called_csv;

    return ($cmd);
}


sub _run_mpileup {
    my $self = shift;
    my $vcf_out = shift;

    my $mpileup_cmd = $self->samtools
       . ' mpileup --skip-indels -d 500 -t INFO/AD,INFO/ADF,INFO/ADR -C50 -u'
       . " -f " . $self->fa_ref
       . ' ' . $self->bam_file
       . ' > ' . $vcf_out;

    if (system($mpileup_cmd)) {
        Pathogens::Exception::SystemCall->throw(error => "Error running mpileup: $mpileup_cmd\n");
    }
}


sub _list_sum_and_max {
    my $self = shift;
    my $list = shift;
    my $sum = 0;
    my $max = 0;
    for my $num (@$list) {
        $sum += $num;
        $max = $num > $max ? $num : $max;
    }
    return ($sum, $max);
}


sub _parse_vcf_line {
    my $self = shift;
    my $line = shift;
    my $min_total_depth = 4;
    my $min_second_depth = 2;
    my $max_allele_freq = 0.9;

    my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format) = split(/\t/, $$line);
    return ($chrom, 0, 0) unless $alt =~/,/;

    if ($info =~ /ADF=([0-9,]+);.*ADR=([0-9,]+);/) {
        my @adf = split(/,/, $1);
        my @adr = split(/,/, $2);
        ($#adf == $#adr) or Pathogens::Exception::VcfParse->throw(error => "Mismatch in ADF and ADR lengths for this VCF line: $$line\n");

        my ($adf_sum, $adf_max) = $self->_list_sum_and_max(\@adf);
        return ($chrom, 0, 0) if $adf_sum < $min_total_depth;

        my ($adr_sum, $adr_max) = $self->_list_sum_and_max(\@adr);
        return ($chrom, 0, 0) if $adr_sum < $min_total_depth;

        my $good_indexes = 0;
        for my $i (0 .. $#adf) {
            $good_indexes++ if ($adf[$i] >= $min_second_depth
                                && $adf[$i] / $adf_sum <= $max_allele_freq
                                && $adr[$i] >= $min_second_depth
                                && $adr[$i] / $adr_sum <= $max_allele_freq);

            return ($chrom, 1, 1) if $good_indexes > 1;
        }

        if ($adf[0] < $adf_max && $adr[0] < $adr_max) {
            return ($chrom, 1, 0);
        }
        else {
            return ($chrom, 0, 0);
        }
    }
    else {
        Pathogens::Exception::VcfParse->throw(error => "Error getting ADF and ADR info from this line: $$line\n");
    }
}


sub _filter_vcf_and_count_snps {
    my $self = shift;
    my $infile = shift;
    my $outfile = shift;
    my $position_count = 0;
    my $het_count = 0;
    my $snp_count = 0;
    my %results;

    open FIN, $infile or Pathogens::Exception::FileIO->throw(error => "Error opening file $infile\n");
    open FOUT, ">$outfile" or Pathogens::Exception::FileIO->throw(error => "Error opening file $outfile\n");

    while (my $line = <FIN>) {
        if ($line =~/^#/) {
            print FOUT $line;
        }
        else {
            my ($chrom, $is_snp, $is_het) = $self->_parse_vcf_line(\$line);
            unless (exists $results{$chrom}) {
                $results{$chrom} = {};
            }
            $results{$chrom}{'positions'}++;
            $results{$chrom}{'snps'}++ if $is_snp;
            if ($is_het) {
                $results{$chrom}{'hets'}++;
                print FOUT $line;
            }
        }
    }

    close FIN or Pathogens::Exception::FileIO->throw(error => "Error closing file $infile\n");
    close FOUT or Pathogens::Exception::FileIO->throw(error => "Error closing file $outfile\n");
    return \%results;
}


sub _lengths_from_fai {
    my $self = shift;
    my $fai = shift;
    my %lengths;
    open F, $fai or Pathogens::Exception::FileIO->throw(error => "Error opening file $fai\n");
    while (<F>) {
        chomp;
        my ($contig, $length, undef, undef, undef) = split /\t/;
        $lengths{$contig} = $length;
    }
    return \%lengths;
}


sub _write_summary_report {
    my $self = shift;
    my $outfile = shift;
    my $totals = shift;

    open F, ">$outfile" or Pathogens::Exception::FileIO->throw(error => "Error opening file $outfile\n");
    print F "No. Het SNPs\t\% Het SNPs (Total Genome)\t\% Het SNPs (Genome Covered)\t\% Het SNPs (Total No. of SNPs)\n";
    my $het_total_genome = sprintf("%.2f", $totals->{length} == 0 ? 0 : 100 * $totals->{hets} / $totals->{length});
    my $het_genome_cov = sprintf("%.2f", $totals->{positions} == 0 ? 0 : 100 * $totals->{hets} / $totals->{positions});
    my $het_total_snp = sprintf("%.2f", $totals->{snps} == 0 ? 0 : 100 * $totals->{hets} / $totals->{snps});
    print F "$totals->{hets}\t$het_total_genome\t$het_genome_cov\t$het_total_snp\n";
    close F or Pathogens::Exception::FileIO->throw(error => "Error closing file $outfile\n");
}


sub _write_reports {
    my $self = shift;
    my $outfile_per_contig = shift;
    my $outfile_summary = shift;
    my $stats = shift;
    my $ctg_lengths = shift;
    my %totals = (length => 0, positions => 0, snps => 0, hets => 0);

    open F, ">$outfile_per_contig" or Pathogens::Exception::FileIO->throw(error => "Error opening file $outfile_per_contig\n");
    print F "Contig\tLength\tPositions_mapped\tSNPs\tHet_SNPs\t\% Het SNPs (Total Contig)\t\% Het SNPs (Contig Covered)\t\% Het SNPs (Total No. of SNPs)\n";

    for my $contig (sort keys %{$ctg_lengths}) {
        $totals{length} += $ctg_lengths->{$contig};

        if (exists $stats->{$contig}) {
            my $ctg_stats = $stats->{$contig};

            my @keys = qw/positions snps hets/;
            for (@keys) {
                $totals{$_} += $ctg_stats->{$_};
            }

            my $het_snp_total_ctg = sprintf("%.2f", 100 * $ctg_stats->{hets} / $ctg_lengths->{$contig});
            my $het_snp_ctg_cov = sprintf("%.2f", $ctg_stats->{positions} == 0 ? 0 : 100 * $ctg_stats->{hets} / $ctg_stats->{positions});
            my $het_snp_total_snp = sprintf("%.2f", $ctg_stats->{snps} == 0 ? 0 : 100 * $ctg_stats->{hets} / $ctg_stats->{snps});
            print F "$contig\t$ctg_lengths->{$contig}\t$ctg_stats->{positions}\t$ctg_stats->{snps}\t$ctg_stats->{hets}\t$het_snp_total_ctg\t$het_snp_ctg_cov\t$het_snp_total_snp\n";
        }
        else {
            print F "$contig\t$ctg_lengths->{$contig}\t0\t0\t0\t0.00\t0.00\t0.00\n";
        }
    }

    close F or Pathogens::Exception::FileIO->throw(error => "Error closing file $outfile_per_contig\n");
    $self->_write_summary_report($outfile_summary, \%totals);
}


sub run {
    my $self = shift;
    my $tmp_vcf = 'fixme';
    $self->_run_mpileup($tmp_vcf);
    my $snp_counts = $self->_filter_vcf_and_count_snps(infile, outfile);
    unlink $tmp_vcf;
    my $contig_lengths = $self->_lengths_from_fai($self->{fa_ref} . 'fai');
    $self->_write_reports(out_per_contig, out_summary, $snp_counts, $contig_lengths);
}


no Moose;
__PACKAGE__->meta->make_immutable;
1;
