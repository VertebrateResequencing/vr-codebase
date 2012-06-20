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
package Pathogens::RNASeq::CommonSpreadsheet;
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
}

sub _header
{
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
