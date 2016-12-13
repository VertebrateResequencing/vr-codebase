package Pathogens::QC::HetSNPCalculator;

use Moose;
use Utils;
use File::Spec;
use Pathogens::Exception;

has 'samtools'         => ( is => 'ro', isa => 'Str', required => 1 );
has 'bcftools'         => ( is => 'ro', isa => 'Str', required => 1 );
has 'min_total_depth'  => ( is => 'ro', isa => 'Int', required => 1 );
has 'min_second_depth' => ( is => 'ro', isa => 'Int', required => 1 );
has 'max_allele_freq'  => ( is => 'ro', isa => 'Num', required => 1 );
has 'fa_ref'           => ( is => 'ro', isa => 'Str', required => 1 );
has 'bam'              => ( is => 'ro', isa => 'Str', required => 1 );
has 'outprefix'        => ( is => 'ro', isa => 'Str', required => 1 );


sub _run_mpileup {
    my $self = shift;
    my $vcf_out = shift;

    my $mpileup_cmd = $self->samtools
       . ' mpileup --skip-indels -d 500 -t INFO/AD,INFO/ADF,INFO/ADR -C50 -uv'
       . " -f " . $self->fa_ref
       . ' ' . $self->bam
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

    my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format) = split(/\t/, $$line);
    return ($chrom, 0, 0) unless $alt =~/,/;

    if ($info =~ /ADF=([0-9,]+);.*ADR=([0-9,]+);/) {
        my @adf = split(/,/, $1);
        my @adr = split(/,/, $2);
        ($#adf == $#adr) or Pathogens::Exception::VcfParse->throw(error => "Mismatch in ADF and ADR lengths for this VCF line: $$line\n");

        my ($adf_sum, $adf_max) = $self->_list_sum_and_max(\@adf);
        return ($chrom, 0, 0) if $adf_sum < $self->min_total_depth;

        my ($adr_sum, $adr_max) = $self->_list_sum_and_max(\@adr);
        return ($chrom, 0, 0) if $adr_sum < $self->min_total_depth;

        my $good_indexes = 0;
        for my $i (0 .. $#adf) {
            $good_indexes++ if ($adf[$i] >= $self->min_second_depth
                                && $adf[$i] / $adf_sum <= $self->max_allele_freq
                                && $adr[$i] >= $self->min_second_depth
                                && $adr[$i] / $adr_sum <= $self->max_allele_freq);

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
                $results{$chrom} = {positions => 0, snps => 0, hets => 0};
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
    my $tmp_vcf = $self->outprefix . '.tmp.vcf';
    $self->_run_mpileup($tmp_vcf);
    my $filtered_vcf = $self->outprefix . '.vcf';
    my $snp_counts = $self->_filter_vcf_and_count_snps($tmp_vcf, $filtered_vcf);
    unlink $tmp_vcf;
    my $contig_lengths = $self->_lengths_from_fai($self->{fa_ref} . '.fai');
    my $out_summary = $self->outprefix . '_report.txt';
    my $out_per_ref_seq = $self->outprefix . '_ref_seq_breakdown.tsv';
    $self->_write_reports($out_per_ref_seq, $out_summary, $snp_counts, $contig_lengths);
    my $vcf =  $self->outprefix . '.vcf';
    my $bcftools_command = $self->bcftools . ' convert -o ' . $self->outprefix . ".bcf -O b $vcf";
    if (system("$bcftools_command")) {
        Pathogens::Exception::SystemCall->throw(error => "Error running bcftools convert: $bcftools_command\n");
    }
    unlink $vcf;
}


no Moose;
__PACKAGE__->meta->make_immutable;
1;
