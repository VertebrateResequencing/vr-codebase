=head1 NAME

ValidateFastqConversion.pm - Checks imported reads math reads from iRODS.

=head1 SYNOPSIS

use Pathogens::Import::ValidateFastqConversion;
my $validate = Pathogens::Import::ValidateFastqConversion->new(
     fastqcheck_filenames => ['1234_5_6_1.fastq.gz.fastqcheck','1234_5_6_2.fastq.gz.fastqcheck'],
     irods_filename       => '1234_5#6.bam'
    );
$validate->is_total_reads_valid();

=cut
package Pathogens::Import::ValidateFastqConversion;

use Moose;
use VertRes::Parser::fastqcheck;
use VertRes::Wrapper::iRODS;

has 'fastqcheck_filenames' => ( is => 'ro', isa => 'ArrayRef', required => 1);
has 'irods_filename'       => ( is => 'ro', isa => 'Str',      required => 1);
has '_fastq_totalreads'    => ( is => 'ro', isa => 'ArrayRef', lazy_build => 1);

sub _build__fastq_totalreads
{
    my ($self) = @_;
    my @totalreads;
    for my $filename (@{$self->fastqcheck_filenames})
    {
	my $parser = VertRes::Parser::fastqcheck->new(file => $filename);
	my $readcount = $parser->num_sequences() || 0;
	push @totalreads, $readcount;
    }
    return \@totalreads;
}

sub _sum_fastq_reads
{
    my ($self) = @_;
    my $sum = 0;
    for my $total_reads (@{$self->_fastq_totalreads})
    {
	$sum += $total_reads;
    }
    return $sum;
}

sub _sum_irods_reads
{
    my ($self) = @_;
    my $irods = VertRes::Wrapper::iRODS->new();
    my $ifile = $irods->find_file_by_name($self->irods_filename);
    my $readcount = $irods->get_total_reads($ifile) || 0;
    return $readcount;
}

sub is_total_reads_valid
{
    my ($self) = @_;
    if($self->_sum_irods_reads == $self->_sum_fastq_reads)
    {
	return 1;
    }
    return 0;
}



1;
