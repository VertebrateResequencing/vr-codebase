package Pathogens::QC::HetSNPCalculator;

use Moose;
use File::Spec;

#executables
has 'samtools' => ( is => 'ro', isa => 'Str', required => 1 );
has 'bcftools' => ( is => 'ro', isa => 'Str', required => 1 );

#bcftools filters
has 'bcft_min_dp' => ( is => 'ro', isa => 'Str', required => 1 );
has 'bcft_min_dv' => ( is => 'ro', isa => 'Str', required => 1 );
has 'bcft_dp_dv_ratio' => ( is => 'ro', isa => 'Str', required => 1 );
has 'bcft_min_qual' => ( is => 'ro', isa => 'Str', required => 1 );
has 'bcft_dp4_ref_allele_ratio' => ( is => 'ro', isa => 'Str', required => 1 );

#Misc info
has 'fa_ref' => ( is => 'ro', isa => 'Str', required => 1 );
has 'reference_size' => ( is => 'ro', isa => 'Str', required => 1 );
has 'lane_path' => ( is => 'ro', isa => 'Str', required => 1 );
has 'lane' => ( is => 'ro', isa => 'Str', required => 1 );
has 'sample_dir' => ( is => 'ro', isa => 'Str', required => 1 );
has 'heterozygosity_report_file_name' => ( is => 'ro', isa => 'Str', required => 1 );

has 'full_path' => ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_full_path' );
has 'mpileup_command' => ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_mpileup_command' );
has 'total_number_of_snps_command' => ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_total_number_of_snps_command' );
has 'snp_call_command' => ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_snp_call_command' );
has 'bcf_query_filter_command' => ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_bcf_query_filter_command' );

sub build_full_path {

  my ($self) = @_;
  my $full_path = File::Spec->catfile($self->lane_path, $self->sample_dir);
  return($full_path);
}

sub _file_path {

  my ($self,$suffix) = @_;
  my $path = File::Spec->catfile($self->full_path, $self->{lane} . $suffix);
  return($path);
}

sub build_mpileup_command {

  my ($self) = @_;

  my $bam_file = _file_path( $self, q(.bam) );
  my $temp_vcf = _file_path( $self, q(_temp_vcf.vcf.gz) );
  my $cmd = $self->samtools . q( mpileup -d 500 -t INFO/DPR,DV -C50 -ugf );
  $cmd .= $self->fa_ref . q( ) . $bam_file . q( | bgzip > ) . $temp_vcf;

  return($cmd);
}


sub build_total_number_of_snps_command {

  my ($self) = @_;

  my $temp_vcf = _file_path( $self, q(_temp_vcf.vcf.gz) );
  my $cmd = $self->{bcftools} . q( call -m -f GQ,GP ) . $temp_vcf . q( | egrep -v "^#|DP=0" | wc -l);

  return($cmd);
}


sub build_snp_call_command {

  my ($self) = @_;

  my $temp_vcf = _file_path( $self, q(_temp_vcf.vcf.gz) );
  my $snp_called_vcf = _file_path( $self, q(_snp_called.vcf.gz) );
  my $cmd = $self->bcftools . q( call -vm -O z ) . $temp_vcf . q( > ) . $snp_called_vcf;

  return($cmd);
}


sub build_bcf_query_filter_command {

  my ($self) = @_;

  my $snp_called_vcf = _file_path( $self, q(_snp_called.vcf.gz) );
  my $filtered_snp_called_vcf = _file_path( $self, q(_filtered_snp_called.vcf) );
  my $cmd = $self->{bcftools} . q( filter -i );
  $cmd .= q{"(DP4[0]+DP4[1])/(DP4[2]+DP4[3]) > 0.3" } . $snp_called_vcf . q{ | };
  $cmd .= $self->{bcftools} . q( filter -i);
  $cmd .= q{ "MIN(DP) >= } . $self->{bcft_min_dp};
  $cmd .= q{ & MIN(DV) >= } . $self->{bcft_min_dv};
  $cmd .= q{ & MIN(DV/DP)>= } . $self->{bcft_dp_dv_ratio};
  $cmd .= q{ & QUAL >= } . $self->{bcft_min_qual};
  $cmd .= q{ & (GT='1/0' | GT='0/1' | GT='1/2')" -};
  $cmd .= q{ > } . $filtered_snp_called_vcf;

  return($cmd);
}

sub get_total_number_of_snps {

  my ($self) = @_;

}

sub get_number_of_het_snps {

  my ($self) = @_;


}

sub calculate_percentage {

  my ($self) = @_;

}

sub write_het_report {

  my ($self) = @_;

}

no Moose;
__PACKAGE__->meta->make_immutable;
1;
