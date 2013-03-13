=head1 NAME

ExpressionStatsSpreadsheet.pm   - Builds a spreadsheet of expression results

=head1 SYNOPSIS

use Pathogens::RNASeq::ExpressionStatsSpreadsheet;
my $expression_results = Pathogens::RNASeq::ExpressionStatsSpreadsheet->new(
  output_filename => '/abc/my_results.csv',
  );
$expression_results->add_result($my_rpkm_values);
$expression_results->add_result($more_rpkm_values);
$expression_results->build_and_close();


=cut
package Pathogens::RNASeq::ExpressionStatsSpreadsheet;
use Moose;
extends 'Pathogens::RNASeq::CommonSpreadsheet';

sub _result_rows
{
  my ($self) = @_;
  my @denormalised_results;
  for my $result_set (@{$self->_results})
  {
    push(@denormalised_results, 
      [
      $result_set->{seq_id},
      $result_set->{gene_id},
      $result_set->{locus_tag},
      $result_set->{mapped_reads_sense},
      $result_set->{rpkm_sense},
      $result_set->{mapped_reads_antisense},
      $result_set->{rpkm_antisense}
      ]);
  }
  return \@denormalised_results;
}

sub _header
{
  my ($self) = @_;
  my @header;
  if($self->protocol ne "StrandSpecificProtocol")
  {
    @header = ["Seq ID","GeneID","Locus Tag","Reads Mapping", "RPKM", "Antisense Reads Mapping", "Antisense RPKM"];
  }
  else
  {
    @header = ["Seq ID","GeneID","Locus Tag","Antisense Reads Mapping", "Antisense RPKM","Reads Mapping", "RPKM"];
  }
  return \@header;
}


1;
