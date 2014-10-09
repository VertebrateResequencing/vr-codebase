=head1 NAME

=head1 SYNOPSIS

use Pathogens::Import::CompressAndValidate;
my $validator = Pathogens::Import::CompressAndValidate->new( irods_filename => $irods, fastq_filenames => \@fastqs );
$validator->is_compressed_and_validated() || $die("Compress and validate failed\n");

=head1 DESCRIPTION

Compress fastq file, add md5 checksum files and check number of reads against iRODS.

=cut

package Pathogens::Import::CompressAndValidate;
use Moose;
use Utils;
use VertRes::Wrapper::fastqcheck;
use Pathogens::Import::ValidateFastqConversion;

has 'irods_filename'  =>       ( is => 'ro', isa => 'Str',      required => 1);
has 'fastq_filenames' =>       ( is => 'ro', isa => 'ArrayRef', required => 1);

sub _compress_and_checksum
{
    my($self) = @_;

    for my $fastq (@{$self->{fastq_filenames}})
    {
      # Checksum fastqs
      Utils::CMD(qq[md5sum $fastq    > $fastq.md5]);
      # Compress fastq
      Utils::CMD(qq[gzip -9 $fastq]);
      # Checksum fastqs
      Utils::CMD(qq[md5sum $fastq.gz > $fastq.gz.md5]);
    }
}

sub _fastqcheck
{
    my($self) = @_;

    for my $fastq (@{$self->{fastq_filenames}})
    {
	# Generate fastqcheckfile
	my $fastqcheck = VertRes::Wrapper::fastqcheck->new();
	$fastqcheck->run($fastq.'.gz', $fastq.'.gz.fastqcheck.tmp');

	$fastqcheck->run_status >= 1 || return 0; 
    }
    return 1;
}

sub _confirm_valid
{
    my($self) = @_;
    my @fastqcheck_tmp = (); # fastqcheck.tmp files.

    for my $fastq (@{$self->{fastq_filenames}})
    {
	push @fastqcheck_tmp, $fastq.'.gz.fastqcheck.tmp';
    }
    
    # Validate against iRODS
    my $validate = Pathogens::Import::ValidateFastqConversion->new(
	fastqcheck_filenames => \@fastqcheck_tmp,
	irods_filename       => $self->irods_filename
    );
    return $validate->is_total_reads_valid(); # Validation check
}

sub is_compressed_and_validated
{
    my($self) = @_;

    # Compress and checksum
    $self->_compress_and_checksum();

    # Generate temp fastqcheck files
    $self->_fastqcheck() || return 0;

    # Check reads from fastqs match reads from iRODS.
    $self->_confirm_valid() || return 0;

    # Rename temp fastqcheck files if valid.
    for my $fastq (@{$self->{fastq_filenames}})
    {
	Utils::CMD(qq[mv $fastq.gz.fastqcheck.tmp $fastq.gz.fastqcheck]);
    }

    return 1;
}

1;
