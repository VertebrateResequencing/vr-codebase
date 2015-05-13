package Pathogens::QC::HetSNPCalculator;

use Moose;

has 'reference_size' => ( is => 'ro', isa => 'Str', required => 1 );
has 'lane_path' => ( is => 'ro', isa => 'Str', required => 1 );
has 'lane' => ( is => 'ro', isa => 'Str', required => 1 );
has 'sample_dir' => ( is => 'ro', isa => 'Str', required => 1 );
has 'full_path' => ( is => 'ro', isa => 'Str', lazy => 1, builder => 'build_full_path' );
has 'temp_vcf_file_path' => ( is => 'ro', isa => 'Str', lazy => 1, builder => 'build_temp_vcf_file_path' );
has 'snp_called_vcf_file_path' => ( is => 'ro', isa => 'Str', lazy => 1, builder => 'build_snp_called_vcf_file_path' );

#my $lane_path = $self->{lane_path};
#my $full_path = $self->{lane_path} . q(/) . $self->{sample_dir} . q(/);
#my $temp_vcf_file = $full_path . $self->{lane} . q(_temp_vcf.vcf.gz);
#my $snp_called_vcf_file = $full_path . $self->{lane} . q(_snp_called.vcf.gz);
#my $filtered_snp_called_vcf_file = $full_path . $self->{lane} . q(_filtered_snp_called.vcf);

sub build_full_path {

  my ($self) = @_;
  my $full_path = $self->lane_path . q(/) . $self->sample_dir . q(/);
  return($full_path);
}

sub build_temp_vcf_file_path {

  my ($self) = @_;
  my $temp_vcf_file = $self->full_path . $self->{lane} . q(_temp_vcf.vcf.gz);
  return($temp_vcf_file);

}

sub build_snp_called_vcf_file_path {

  my ($self) = @_;
  my $snp_called_vcf_file = $self->full_path . $self->lane . q(_snp_called.vcf.gz);
  return($snp_called_vcf_file);
}

sub build_filtered_snp_called_file_path {

  my ($self) = @_;
  my $filtered_snp_called_vcf_file_path = $self->full_path . $self->lane . q(_filtered_snp_called.vcf);
  return($filtered_snp_called_vcf_file_path)

}

sub build_mpileup_command {


}


sub build_total_number_of_snps_command {


}


sub build_snp_call_command {



}


sub build_bcf_query_filter_command {


}


sub build_vcf_header_command {


}


sub get_total_number_of_snps {

}


sub get_number_of_het_snps {


}

sub calculate_percentage {


}

sub write_het_report {


}

no Moose;
__PACKAGE__->meta->make_immutable;
1;
