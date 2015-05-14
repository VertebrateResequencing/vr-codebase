package Pathogens::QC::HetSNPCalculator;

use Moose;
use Utils;
use File::Spec;
use Pathogens::Exception;

#executables
has 'samtools' => ( is => 'ro', isa => 'Str', required => 1 );
has 'bcftools' => ( is => 'ro', isa => 'Str', required => 1 );

#bcftools filters
has 'min_rawReadDepth' => ( is => 'ro', isa => 'Str', required => 1 );
has 'min_hqNonRefBases' => ( is => 'ro', isa => 'Str', required => 1 );
has 'rawReadDepth_hqNonRefBases_ratio' => ( is => 'ro', isa => 'Str', required => 1 );
has 'min_qual' => ( is => 'ro', isa => 'Str', required => 1 );
has 'hqRefReads_hqAltReads_ratio' => ( is => 'ro', isa => 'Str', required => 1 );

#Misc info
has 'fa_ref' => ( is => 'ro', isa => 'Str', required => 1 );
has 'reference_size' => ( is => 'ro', isa => 'Str', required => 1 );
has 'lane_path' => ( is => 'ro', isa => 'Str', required => 1 );
has 'lane' => ( is => 'ro', isa => 'Str', required => 1 );
has 'sample_dir' => ( is => 'ro', isa => 'Str', required => 1 );
has 'het_report' => ( is => 'ro', isa => 'Str', required => 1 );

has 'full_path' => ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_full_path' );
has 'het_report_path' => ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_het_report_path' );
has 'mpileup_command' => ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_mpileup_command' );
has 'total_number_of_snps_command' => ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_total_number_of_snps_command' );
has 'snp_call_command' => ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_snp_call_command' );
has 'bcf_query_command' => ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_bcf_query_command' );

has 'number_of_het_snps' => ( is => 'rw', isa => 'Str', lazy => 1, default => 0 );
has 'total_number_of_snps' => ( is => 'rw', isa => 'Str', lazy => 1, default => 0 );
has 'het_snps_genome_percentage' => ( is => 'rw', isa => 'Str', lazy => 1, default => 0 );
has 'het_snps_total_snps_percentage' => ( is => 'rw', isa => 'Str', lazy => 1, default => 0 );

sub build_full_path {

  my ($self) = @_;
  my $full_path = File::Spec->catfile( $self->lane_path, $self->sample_dir );
  return($full_path);
}

sub build_het_report_path {

  my ($self) = @_;
  my $het_report_path = File::Spec->catfile( $self->lane_path, $self->het_report );
  return($het_report_path);
}

sub _file_path {

  my ($self,$suffix) = @_;
  my $path = File::Spec->catfile( $self->full_path, $self->{lane} . $suffix );
  return($path);
}

sub build_mpileup_command {

  my ($self) = @_;

  my $bam_file = _file_path( $self, q(.bam) );
  my $temp_vcf = _file_path( $self, q(_temp_vcf.vcf.gz) );
  my $cmd = $self->samtools;
  $cmd .= q( mpileup -d 500 -t INFO/DPR,DV -C50 -ugf );
  $cmd .= $self->fa_ref;
  $cmd .= q( );
  $cmd .= $bam_file;
  $cmd .= q( | bgzip > );
  $cmd .= $temp_vcf;

  return($cmd);
}


sub build_total_number_of_snps_command {

  my ($self) = @_;

  my $temp_vcf = _file_path( $self, q(_temp_vcf.vcf.gz) );
  my $total_number_of_snps = _file_path( $self, q(_total_number_of_snps.csv) );

  my $cmd = $self->{bcftools};
  $cmd .= q( query -f "%CHROM\n");
  $cmd .= q( -i "DP > 0" );
  $cmd .= $temp_vcf;
  $cmd .= q( > );
  $cmd .= $total_number_of_snps;

  return($cmd);
}


sub build_snp_call_command {

  my ($self) = @_;

  my $temp_vcf = _file_path( $self, q(_temp_vcf.vcf.gz) );
  my $snp_called_vcf = _file_path( $self, q(_snp_called.vcf.gz) );

  my $cmd = $self->bcftools;
  $cmd .= q( call -vm -O z );
  $cmd .= $temp_vcf;
  $cmd .= q( > );
  $cmd .= $snp_called_vcf;

  return($cmd);
}


sub build_bcf_query_command {

  my ($self) = @_;

  my $snp_called_vcf = _file_path( $self, q(_snp_called.vcf.gz) );
  my $filtered_snp_called_csv = _file_path( $self, q(_filtered_snp_called_list.csv) );

  my $cmd = $self->{bcftools} . q( query -f);
  $cmd .= q{ "%CHROM %POS\n" -i};
  $cmd .= q{ "MIN(DP) >= } . $self->{min_rawReadDepth};
  $cmd .= q{ & MIN(DV) >= } . $self->{min_hqNonRefBases};
  $cmd .= q{ & MIN(DV/DP)>= } . $self->{rawReadDepth_hqNonRefBases_ratio};
  $cmd .= q{ & QUAL >= } . $self->{min_qual};
  $cmd .= q{ & (GT='1/0' | GT='0/1' | GT='1/2')};
  $cmd .= q{ & ((DP4[0]+DP4[1])/(DP4[2]+DP4[3]) > } . $self->{hqRefReads_hqAltReads_ratio} . q{)" };
  $cmd .= $snp_called_vcf;
  $cmd .= q{ > } . $filtered_snp_called_csv;

  return($cmd);
}

sub get_total_number_of_snps {

  my ($self) = @_;

  my $total_number_of_snps = _file_path( $self, q(_total_number_of_snps.csv) );

  Utils::CMD($self->total_number_of_snps_command);
  #my $exit_code = system($self->total_number_of_snps_command);

  if ( -e $total_number_of_snps) {
    open (my $fh, '<', $total_number_of_snps) or Utils::error("$total_number_of_snps: $!");
    $self->total_number_of_snps( _count_file_rows($self,$fh) );
  }
  else {
    Pathogens::Exception::HetSNPStepCommand->throw(
						   error => "A problem occured running the total number of snps command '" .
						   $self->total_number_of_snps_command .
						   "'\nThe file '$total_number_of_snps' was not created"
						  );
  }
}

sub get_number_of_het_snps {

  my ($self) = @_;

  my $temp_vcf = _file_path( $self, q(_temp_vcf.vcf.gz) );
  Utils::CMD($self->mpileup_command);

  if ( -e $temp_vcf) {
    my $snp_called_vcf = _file_path( $self, q(_snp_called.vcf.gz) );
    Utils::CMD($self->snp_call_command);

    if ( -e $snp_called_vcf) {
      my $filtered_snp_called_csv = _file_path( $self, q(_filtered_snp_called_list.csv) );
      Utils::CMD($self->bcf_query_command);

      if ( -e $filtered_snp_called_csv) {
	open (my $fh, '<', $filtered_snp_called_csv) or Utils::error("$filtered_snp_called_csv: $!");
	$self->number_of_het_snps( _count_file_rows($self,$fh) );
      }
      else {
	Pathogens::Exception::HetSNPStepCommand->throw(
						     error => "A problem occured running the bcf query command '" .
						     $self->bcf_query_command .
						     "'\nThe file '$filtered_snp_called_csv' was not created"
						    );
      }
    }
    else {
      Pathogens::Exception::HetSNPStepCommand->throw(
						   error => "A problem occured running the SNP calling command '" .
						   $self->snp_call_command .
						   "'\nThe file '$snp_called_vcf' was not created"
						  );
    }
  }
  else {
    Pathogens::Exception::HetSNPStepCommand->throw(
						   error => "A problem occured running the Mpileup command '" .
						   $self->mpileup_command .
						   "'\nThe file '$temp_vcf' was not created"
						  );
  }
}

sub get_percentages_of_het_snps {

  my ($self) = @_;
  $self->het_snps_genome_percentage( _calculate_percentage( $self,$self->number_of_het_snps,$self->reference_size ) );
  $self->het_snps_total_snps_percentage( _calculate_percentage( $self,$self->number_of_het_snps,$self->total_number_of_snps ) );

}

sub _calculate_percentage {

  my ($self,$count,$total) = @_;
  return( ($count*100)/$total ) if ($total > 0);
  Pathogens::Exception::NullDenominator->throw( error => "The value used as total is null (0). No percentage can be calculated");

}


sub write_het_report {

  my ($self) = @_;
  my $lane_path = $self->lane_path;
  my $number_of_het_snps = $self->number_of_het_snps;
  my $het_snps_genome_percentage = $self->het_snps_genome_percentage;
  my $het_snps_total_snps_percentage = $self->het_snps_total_snps_percentage;

  open(my $fh, '>', $self->het_report_path) or Utils::error($self->het_report_path . ": $!");
  print $fh "$number_of_het_snps\t$het_snps_genome_percentage\t$het_snps_total_snps_percentage\n";
  close($fh);
}

sub remove_temp_vcfs_and_csvs {

  my ($self) = @_;

  my $temp_vcf = _file_path( $self, q(_temp_vcf.vcf.gz) );
  my $snp_called_vcf = _file_path( $self, q(_snp_called.vcf.gz) );
  my $filtered_snp_called_csv = _file_path( $self, q(_filtered_snp_called_list.csv) );
  my $total_number_of_snps = _file_path( $self, q(_total_number_of_snps.csv) );
  unlink($temp_vcf);
  unlink($snp_called_vcf);
  unlink($filtered_snp_called_csv);
  unlink($total_number_of_snps);

}

sub _count_file_rows {

  my ($self,$fh) = @_;
  my $number_of_rows = 0;
  while( my $row = <$fh> ) {
    chomp($row);
    if($row && $row ne q()) {
      $number_of_rows++;
    }
  }
  close($fh);
  return $number_of_rows;
}


no Moose;
__PACKAGE__->meta->make_immutable;
1;
