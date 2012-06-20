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
package Pathogens::RNASeq::InsertionStatsSpreadsheet;
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
      $result_set->{forward_insert_sites},
      $result_set->{normalised_forward_insert_sites},
      $result_set->{reverse_insert_sites},
      $result_set->{normalised_reverse_insert_sites},
      $result_set->{total_insert_sites},
      $result_set->{normalised_total_insert_sites}
      ]);

  }
  return \@denormalised_results;
}

sub _header
{
  my ($self) = @_;
  my @header;
  @header = ["Seq ID","GeneID","Forward insert sites", "Normalised forward insert sites", "Reverse insert sites", "Normalised reverse insert sites", "Total insert sites", "Normalised total insert sites"];  
  return \@header;
}


1;
