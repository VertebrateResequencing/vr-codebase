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
use Text::CSV;


has 'output_filename'     => ( is => 'rw', isa => 'Str',      required   => 1 );
has '_output_file_handle' => ( is => 'rw', lazy_build => 1 );
has '_results'            => ( is => 'rw', isa => 'ArrayRef', lazy_build => 1 );
has 'protocol'            => ( is => 'rw', isa => 'Str',      default    => 'StandardProtocol' );

sub _build__output_file_handle
{
  my ($self) = @_;
  my $output_file_handle;
  open($output_file_handle, ">:encoding(utf8)",  $self->output_filename ) || Pathogens::RNASeq::Exceptions::FailedToOpenExpressionResultsSpreadsheetForWriting->throw( error => "Cant open ".$self->output_filename." for writing");
  return $output_file_handle;
}

sub _build__results
{
  my ($self) = @_;
  my @results = () ;
  return \@results;
}

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
  if($self->protocol eq "StrandSpecificProtocol")
  {
    @header = ["Seq ID","GeneID","Reads Mapping", "RPKM", "Antisense Reads Mapping", "Antisense RPKM"];
  }
  else
  {
    @header = ["Seq ID","GeneID", "Antisense Reads Mapping", "Antisense RPKM","Reads Mapping", "RPKM"];
  }
  return \@header;
}


sub add_result
{
  my ($self,$result_hash) = @_;
  push(@{$self->_results}, $result_hash );
  return 1;
}

sub build_and_close
{
  my ($self) = @_;

  my $csv = Text::CSV->new ( { binary => 1 } ) || Pathogens::RNASeq::Exceptions::FailedToOpenExpressionResultsSpreadsheetForWriting->throw( error => "Cant open ".$self->output_filename." for writing ".Text::CSV->error_diag());
  $csv->eol("\r\n");

  # print out the header
  $csv->print($self->_output_file_handle, @{$self->_header});

  # print out all the results rows
  for my $row (@{$self->_result_rows})
  {
    $csv->print($self->_output_file_handle, $row);
  }
  
  close($self->_output_file_handle);  
  return 1;
}
1;
