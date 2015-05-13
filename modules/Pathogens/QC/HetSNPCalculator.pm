package Pathogens::QC::HetSNPCalculator;

use Moose;

has 'reference_size' => ( is => 'ro', isa => 'Str', required => 1 );

#my $lane_path = $self->{lane_path};
#my $full_path = $self->{lane_path} . q(/) . $self->{sample_dir} . q(/);
#my $temp_vcf_file = $full_path . $self->{lane} . q(_temp_vcf.vcf.gz);
#my $snp_called_vcf_file = $full_path . $self->{lane} . q(_snp_called.vcf.gz);
#my $filtered_snp_called_vcf_file = $full_path . $self->{lane} . q(_filtered_snp_called.vcf);



sub build_filenames {


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
