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

#Misc required info
has 'fa_ref' => ( is => 'ro', isa => 'Str', required => 1 );
has 'reference_size' => ( is => 'ro', isa => 'Str', required => 1 );
has 'lane_path' => ( is => 'ro', isa => 'Str', required => 1 );
has 'lane' => ( is => 'ro', isa => 'Str', required => 1 );
has 'sample_dir' => ( is => 'ro', isa => 'Str', required => 1 );
has 'het_report' => ( is => 'ro', isa => 'Str', required => 1 );

#File path builders
has 'full_path' => ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_full_path' );
has 'bam_file' => ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_bam_file' );
has 'temp_vcf' => ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_temp_vcf' );
has 'snp_called_vcf' => ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_snp_called_vcf' );
has 'filtered_snp_called_csv' => ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_filtered_snp_called_csv' );
has 'total_number_of_snps_csv' => ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_total_number_of_snps_csv' );
has 'het_report_path' => ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_het_report_path' );

#Command string builders
has 'mpileup_command' => ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_mpileup_command' );
has 'total_number_of_snps_command' => ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_total_number_of_snps_command' );
has 'snp_call_command' => ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_snp_call_command' );
has 'bcf_query_command' => ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_bcf_query_command' );

#Calculations builders
has 'number_of_het_snps' => ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_number_of_het_snps' );
has 'total_number_of_snps' => ( is => 'rw', isa => 'Str', lazy => 1, builder => 'build_total_number_of_snps' );

sub build_full_path {

  my ($self) = @_;
  my $full_path = File::Spec->catfile( $self->lane_path, $self->sample_dir );
  return($full_path);
}

sub build_bam_file {

  my ($self) = @_;
  my $path = File::Spec->catfile( $self->full_path, $self->{lane} . q(.bam) );
  return($path);
}

sub build_temp_vcf {

  my ($self) = @_;
  my $path = File::Spec->catfile( $self->full_path, $self->{lane} . q(_temp_vcf.vcf.gz) );
  return($path);
}

sub build_snp_called_vcf {

  my ($self) = @_;
  my $path = File::Spec->catfile( $self->full_path, $self->{lane} . q(_snp_called.vcf.gz) );
  return($path);
}

sub build_filtered_snp_called_csv {

  my ($self) = @_;
  my $path = File::Spec->catfile( $self->full_path, $self->{lane} . q(_filtered_snp_called_list.csv) );
  return($path);
}

sub build_total_number_of_snps_csv {

  my ($self) = @_;
  my $path = File::Spec->catfile( $self->full_path, $self->{lane} . q(_total_number_of_snps.csv) );
  return($path);
}

sub build_het_report_path {

  my ($self) = @_;
  my $het_report_path = File::Spec->catfile( $self->lane_path, $self->lane . q(_) . $self->het_report );
  return($het_report_path);
}

sub build_mpileup_command {

  my ($self) = @_;

  my $cmd = $self->samtools;
  $cmd .= q( mpileup -d 500 -t INFO/DPR,DV -C50 -ugf );
  $cmd .= $self->fa_ref;
  $cmd .= q( );
  $cmd .= $self->bam_file;
  $cmd .= q( | bgzip > );
  $cmd .= $self->temp_vcf;

  return($cmd);
}

sub build_total_number_of_snps_command {

  my ($self) = @_;

  my $cmd = $self->{bcftools};
  $cmd .= q( query -f "%CHROM\n");
  $cmd .= q( -i "DP > 0" );
  $cmd .= $self->temp_vcf;
  $cmd .= q( > );
  $cmd .= $self->total_number_of_snps_csv;

  return($cmd);
}

sub build_snp_call_command {

  my ($self) = @_;

  my $cmd = $self->bcftools;
  $cmd .= q( call -vm -O z );
  $cmd .= $self->temp_vcf;
  $cmd .= q( > );
  $cmd .= $self->snp_called_vcf;

  return($cmd);
}

sub build_bcf_query_command {

  my ($self) = @_;

  my $cmd = $self->{bcftools} . q( query -f);
  $cmd .= q{ "%CHROM %POS\n" -i};
  $cmd .= q{ "MIN(DP) >= } . $self->{min_rawReadDepth};
  $cmd .= q{ & MIN(DV) >= } . $self->{min_hqNonRefBases};
  $cmd .= q{ & MIN(DV/DP)>= } . $self->{rawReadDepth_hqNonRefBases_ratio};
  $cmd .= q{ & QUAL >= } . $self->{min_qual};
  $cmd .= q{ & (GT='1/0' | GT='0/1' | GT='1/2')};
  $cmd .= q{ & ((DP4[0]+DP4[1])/(DP4[2]+DP4[3]) > } . $self->{hqRefReads_hqAltReads_ratio} . q{)" };
  $cmd .= $self->snp_called_vcf;
  $cmd .= q{ > } . $self->filtered_snp_called_csv;

  return($cmd);
}

sub build_total_number_of_snps {

  my ($self) = @_;

  Utils::CMD($self->total_number_of_snps_command);

  if ( -e $self->total_number_of_snps_csv) {
    open (my $fh, '<', $self->total_number_of_snps_csv) or Utils::error($self->total_number_of_snps_csv . ": $!");
    return( _count_file_rows($self,$fh) );
  }
  else {
    Pathogens::Exception::HetSNPStepCommand->throw(
						   error => "A problem occured running the total number of snps command '" .
						   $self->total_number_of_snps_command .
						   "'\nThe file '" .
						   $self->total_number_of_snps_csv .
						   "' was not created"
						  );
  }
}

sub build_number_of_het_snps {

  my ($self) = @_;

  Utils::CMD($self->mpileup_command);

  if ( -e $self->temp_vcf) {

    Utils::CMD($self->snp_call_command);

    if ( -e $self->snp_called_vcf) {
      Utils::CMD($self->bcf_query_command);

      if ( -e $self->filtered_snp_called_csv) {
	open (my $fh, '<', $self->filtered_snp_called_csv) or Utils::error($self->filtered_snp_called_csv . ": $!");
	return( _count_file_rows($self,$fh) );
      }
      else {
	Pathogens::Exception::HetSNPStepCommand->throw(
						       error => "A problem occured running the bcf query command '" .
						       $self->bcf_query_command .
						       "'\nThe file '" .
						       $self->filtered_snp_called_csv . 
						       "' was not created"
						    );
      }
    }
    else {
      Pathogens::Exception::HetSNPStepCommand->throw(
						     error => "A problem occured running the SNP calling command '" .
						     $self->snp_call_command .
						     "'\nThe file '" .
						     $self->snp_called_vcf .
						     "' was not created"
						  );
    }
  }
  else {
    Pathogens::Exception::HetSNPStepCommand->throw(
						   error => "A problem occured running the Mpileup command '" .
						   $self->mpileup_command .
						   "'\nThe file '" .
						   $self->temp_vcf .
						   "' was not created"
						  );
  }
}

sub write_het_report {

  my ($self) = @_;

  my $het_snps_genome_percentage = _calculate_percentage( $self,$self->number_of_het_snps,$self->reference_size );
  my $het_snps_total_snps_percentage = _calculate_percentage( $self,$self->number_of_het_snps,$self->total_number_of_snps );

  open(my $fh, '>', $self->het_report_path) or Utils::error($self->het_report_path . ": $!");
  print $fh "total_number_of_het_snps\t\%_of_het_snps_for_total_genome_length\t\%_of_het_snps_for_total_snps_found\n";
  print $fh ($self->number_of_het_snps,"\t$het_snps_genome_percentage\t$het_snps_total_snps_percentage\n");
  close($fh);
}

sub remove_temp_vcfs_and_csvs {

  my ($self) = @_;

  unlink($self->temp_vcf);
  unlink($self->snp_called_vcf);
  unlink($self->filtered_snp_called_csv);
  unlink($self->total_number_of_snps_csv);

}

sub _calculate_percentage {

  my ($self,$count,$total) = @_;
  return( ($count*100)/$total ) if ($total > 0);
  Pathogens::Exception::NullDenominator->throw( error => "The value used as total is null (0). No percentage can be calculated\n");

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
