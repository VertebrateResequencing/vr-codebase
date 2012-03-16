package Pathogens::Reports::Mapping::Spreadsheet;

use Moose;
use Pathogens::Reports::Mapping::Row;
use Text::CSV;

has 'filehandle' => ( is => 'ro', isa => 'FileHandle',                                 required => 1 ); # Output filehandle
has 'rows'       => ( is => 'ro', isa => 'ArrayRef[Pathogens::Reports::Mapping::Row]', required => 1 ); # Output rows
has '_columns'   => ( is => 'ro', isa => 'ArrayRef[Str]', lazy_build => 1 ); # Row attributes to output
has '_headers'   => ( is => 'ro', isa => 'ArrayRef[Str]', lazy_build => 1 ); # Column headers for Row attributes
has '_csv_out'   => ( is => 'ro', isa => 'Text::CSV',     lazy_build => 1 ); # Output CSV


sub _build__headers
{
    my($self) = @_;
    my @headers = ('Study ID','Sample', 'Lane','Cycles','Yield (Reads)','Yield (Bases)','Type (QC/Mapping)',
		   'Reference','Reference Size','Mapper','Adapter (%)','Transposon (%)','Mapped (%)','Paired (%)',
		   'Mean Insert Size','Genome Covered (%)','Genome Covered (% >= 1X)','Genome Covered (% >= 5X)',
		   'Genome Covered (% >= 10X)','Genome Covered (% >= 50X)','Genome Covered (% >= 100X)',
		   'Depth of Coverage (X)','Duplication Rate','Error Rate','NPG QC','Manual QC');
    return \@headers;
}

sub _build__columns
{
    my($self) = @_;
    my @columns = ('study_id','sample','lanename','cycles','reads','bases','map_type',
		   'reference','reference_size','mapper','adapter_perc','transposon_perc','mapped_perc','paired_perc',
		   'mean_insert_size','genome_covered','genome_covered_1x','genome_covered_5x',
		   'genome_covered_10x','genome_covered_50x','genome_covered_100x',
		   'depth_of_coverage','duplication_rate','error_rate','npg_qc','manual_qc');          
    return \@columns;
}

sub _build__csv_out
{
    my($self) = @_;
    my $csv = Text::CSV->new( { binary => 1 } );
    $csv->eol("\r\n");
    return $csv;
}

sub _output_headers
{
    my($self) = @_;
    return $self->_csv_out->print($self->filehandle, $self->_headers);
}

sub _output_rows
{
    my($self) = @_;
    my @csv_output = ();

    # get cell values
    for my $row (@{$self->rows})
    {
	my @row_output;
	for my $col (@{$self->_columns})
	{
	    # set any undefined values to NA and output
	    push @row_output, defined($row->$col) ? $row->$col:'NA';
	}
	$self->_csv_out->print($self->filehandle, \@row_output) || return 0;
    }
    return 1;
}

sub output_csv
{
    my($self) = @_; 

    # Handle errors here...
    # No rows ?
    # Filehandle defined ?
    # 

    $self->_output_headers();
    $self->_output_rows();
    return 1;
}

1;
