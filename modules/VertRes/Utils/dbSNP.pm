=head1 NAME

VertRes::Utils::dbSNP - dbSNP submission utility functions

=head1 SYNOPSIS

use VertRes::Utils::dbSNP;

my $dbsnp_util = VertRes::Utils::dbSNP->new();

# use any of the utility functions described here, eg.
my $are_similar = $dbsnp_util->bams_are_similar(@bam_files);

=head1 DESCRIPTION

General utility functions for making dbSNP submission files

=head1 AUTHOR

Thomas Keane: thomas.keane@sanger.ac.uk

=cut

package VertRes::Utils::dbSBP;

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Getopt::Long;
use VertRes::IO;

use base qw(VertRes::Base);

our $REFERENCE_ASSEMBLY = 'NCBIM37';
our $FIVE_PRIME_REGION = 200;
our $THREE_PRIME_REGION = 200;

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Utils::dbSNP->new();
 Function: Create a new VertRes::Utils::dbSNP object.
 Returns : VertRes::Utils::dbSNP object
 Args    : n/a

=cut

sub new 
{
    my ($class, %args) = @_;
    my $self = $class->SUPER::new(%args);
    	    
	checkFields( $self );
	
    return $$self;
}

sub _checkFields
{
	my $self = $_[ 0 ];
	
	if( ! defined( $self->handle ) )
	{
		$self->throw("You must pass in a dbSNP handle at object creation!");
	}
	
	if( ! defined( $self->strain_tag ) )
	{
		$self->throw("You must pass in a strain tag at object creation!");
	}
	
	if( ! defined( $self->pop_handle ) )
	{
		$self->throw("You must pass in a population handle at object creation!");
	}
	
	if( ! defined( $self->species ) )
	{
		$self->throw("You must pass in a species name at object creation!");
	}
	
	if( $self->species ne 'Mouse' )
	{
		$self->throw( "Sorry - only valid species at the moment is Mouse\n");
	}
}

=head2 write_header

 Title   : write_header
 Usage   : my $write = $obj->write_header($self, $output_file, $name, $fax, $fone, $email, $lab, $title, $authors, $year, $status);
 Function: write the dbSNP submission file header
 Returns : n/a
 Args    : output file and various contact info about the submitter

=cut

sub write_header 
{
	my ($self, $outputFile, $name, $fax, $fone, $email, $lab, $title, $authors, $year, $status) = @_;
	
	$self->throw( "Required fields not defined - see perldoc!\n" ) unless @_ == 11;
	checkFields( $self );
    
	open( my $ofh, ">>$outputFile" ) or $self->throw("Cannot create $outputFile: $!");
    print $ofh qq[
TYPE:   CONT
HANDLE: $self->handle
NAME:   $name
FAX:    $fone
TEL:    $fax
EMAIL:  $email
LAB:    $lab
INST:   Wellcome Trust Sanger Institute
ADDR:   Wellcome Trust Genome Campus, Hinxton, Cambridge, CB10 1SA, UK
||
TYPE:    PUB
HANDLE:  $self->handle
TITLE:   $title
AUTHORS:
$authors
YEAR:    $year
STATUS:  $status
||
];
	close( $ofh );
	$self->title = $title;
}

=head2 write_methods

 Title   : write_methods
 Usage   : my $write = $obj->write_methods($self, $outputFile );
 Function: write the dbSNP submission file methods section
 Returns : n/a
 Args    : 

=cut

sub write_methods
{
	my ($self, $outputFile, $id, $methods) = @_;
	
	$self->throw( "Required fields not defined - see perldoc!\n" ) unless @_ == 2;
	checkFields( $self );
	
	open( my $ofh, ">>$outputFile" ) or $self->throw("Cannot create $outputFile: $!");
    print $ofh qq[
TYPE:                   METHOD
HANDLE:                 $self->handle
ID:                     $id
METHOD_CLASS:           Sequence
SEQ_BOTH_STRANDS:       YES
TEMPLATE_TYPE:          DIPLOID
MULT_PCR_AMPLIFICATION: NA
MULT_CLONES_TESTED:     NA
METHOD:
$methods
];
	close( $ofh );
	$self->method_id=$id;
}

=head2 write_population

 Title   : write_population
 Usage   : my $write = $obj->write_population($self, $name, $fax, $fone, $email, $lab, $title, $authors, $year, $status);
 Function: write the population section of a dbSNP submission
 Returns : n/a
 Args    : 

=cut

sub write_population
{
    my ($self, $outputFile ) = @_;
	
	$self->throw( "Required fields not defined - see perldoc!\n" ) unless @_ == 1;
	checkFields( $self );
	
	open( my $ofh, ">>$outputFile" ) or $self->throw("Cannot create $outputFile: $!");
	
    print $ofh qq[
TYPE:       POPULATION
HANDLE:     $self->handle
ID:         $self->pop_handle
POP_CLASS:
POPULATION:
||
];
	close( $ofh );
}

=head2 write_methods

 Title   : write_methods
 Usage   : my $write = $obj->write_methods($self, $file_handle, $name, $fax, $fone, $email, $lab, $title, $authors, $year, $status);
 Function: write the methods section
 Returns : n/a
 Args    : 

=cut

sub write_individual
{
	my ($self, $outputFile ) = @_;
	
	$self->throw( "Required fields not defined - see perldoc!\n" ) unless @_ == 1;
	checkFields( $self );
	
	open( my $ofh, ">>$outputFile" ) or $self->throw("Cannot create $outputFile: $!");

	# TYPE: INDIVIDUAL
	# IND: handle | loc_pop_id | loc_ind_id | tax_id | sex | breed_structure | ind_grp
	# SOURCE: src_type | source_name | src_ind_id | loc_ind_grp
	
    print $ofh qq[
TYPE:   INDIVIDUAL
IND:    $self->handle|$self->pop_handle|$self->strain_tag|10090|Unknown|I|NA
SOURCE: repository|Wellcome Trust Sanger Institute|$self->strain_tag|NA
||
];
	close( $ofh );
}

=head2 write_snpassay_section

 Title   : write_methods
 Usage   : my $write = $obj->write_methods($self, $file_handle, $name, $fax, $fone, $email, $lab, $title, $authors, $year, $status);
 Function: Check if multiple name-sorted bam files are the same, ignoring
           mapping quality 0 reads.
 Returns : boolean (false if they differ)
 Args    : list of bam files

=cut

sub write_snpassay
{
	my ($self, $outputFile, $batch_name ) = @_;
	
	$self->throw( "Required fields not defined - see perldoc!\n" ) unless @_ == 3;
	checkFields( $self );
	
	open( my $ofh, ">>$outputFile" ) or $self->throw("Cannot create $outputFile: $!");
	
	$self->throw("You must call the write_header function first to define the title!") unless defined( $self->title );
	
    print $ofh qq[
TYPE:       SNPASSAY
HANDLE:     $self->handle
BATCH:      $batch_name
MOLTYPE:    Genomic
SAMPLESIZE: 2
METHOD:     Genomic Sequencing
ORGANISM:   Mus musculus
CITATION:   $self->title
||
];
	close( $ofh );
}

=head2 write_snp_records

 Title   : write_methods
 Usage   : my $write = $obj->write_methods($self, $file_handle, $name, $fax, $fone, $email, $lab, $title, $authors, $year, $status);
 Function: Check if multiple name-sorted bam files are the same, ignoring
           mapping quality 0 reads.
 Returns : boolean (false if they differ)
 Args    : list of bam files

=cut

sub write_snp_records
{
	my ($self, $outputFile, $snp_file ) = @_;
	
	# read in the file containing the SNP data
	open( my $sfh, $snp_file ) or $self->throw("Cannot open $snp_file: $!");
	
	open( my $ofh, ">>$outputFile" ) or $self->throw("Cannot create $outputFile: $!");
	
	my $snpcnt = 0;
	# grab a connection to an EnsEMBL mouse db
	my $registry = 'Bio::EnsEMBL::Registry';
	
	$registry->load_registry_from_db(
		-host    => 'ensdb-archive',
		-user    => 'ensro',
		-port    => '5304',
		-verbose => 0
	);
	
	my $slice_adaptor = $registry->get_adaptor( $self->species, 'Core', 'Slice' );
	
	#loop over the SNPs and write out the sequence
	while ( <sfh> )
	{
		# ignore any blank lines
		if ( /\^\n/ ) 
		{
			next;
		}
		chomp;
		
		# create a simple hash record of the current snp
		my @s = split( /\t/, $_ );
		$self->throw("Cannot create $outputFile: $!") unless @s == 6;
		
		my ($chromosome, $position, $reference_base, $consensus_base, $phred_quality, $read_depth ) = @s;
		
		$self->throw("Invalid SNP entry in file: $_") unless @s == 6; #do more sanity checks here.....
		
		# grab the entire sequence once instead of having to make
		# multiple accesses to the db i.e. for upstream and downstream say
		my $start_pos = $position - $FIVE_PRIME_REGION;
		my $end_pos = $position + $THREE_PRIME_REGION;
		
		# grab the SNP and surrounding genomic region based on its 
		# chromosomal location
		my $slice = $slice_adaptor->fetch_by_region( 'chromosome',
					      $chromosome,
					      $start_pos,
					      $end_pos,
					      1,
					      $REFERENCE_ASSEMBLY);
		my $full_seq = $slice->seq();
		
		# determine the supercontig
		my $supercontig_projection = $slice->project('supercontig', $REFERENCE_ASSEMBLY);
		
		# really there should be only 1 supercontig but check for multiples just 
		# in case
		
		my $number_of_supercontigs = scalar @$supercontig_projection;
		
		if ( $number_of_supercontigs == 0 ) 
		{
			$self->throw("No supercontigs found for $_\n");
		}
		elsif ( $number_of_supercontigs < 1 ) 
		{
			$self->throw("Multiple supercontigs found for $_\n");
		}
		
		my $contig;
		my $supercontig;
		foreach my $curr_superc ( @$supercontig_projection ) 
		{
			$contig = $curr_superc->to_Slice();
			$supercontig = $contig->seq_region_name;
		}
		
		# a quick sanity check that the base in the defined snp position
		# is the identical to one defined in the read in data file

		my $snp_seq = substr($full_seq,$FIVE_PRIME_REGION,1);

		if ( $snp_seq ne $reference_base ) 
		{
			$self->throw("Snp reference_base does not match sequence retrieved from EnsEMBL: $_\n");
		}
		
		# extract the upstream and downstream sequences
		my $five_prime_seq = substr($full_seq,0,$FIVE_PRIME_REGION);

		# indexing from zero so add a 1 to the starting position
		my $three_prime_seq = substr($full_seq,$FIVE_PRIME_REGION+1,$THREE_PRIME_REGION);
		
		# local ID for marker: concatenation of handle, strain, count
		# computed location of the SNP (ACCESSION + LOCATION)
		# flanking sequence (200 bp 5' and 200 bp 3' of the variant position)
		
		my $snp_id = $self->handle.'_'.$->strain_tag.'_'.$snpcnt;
		
		print $ofh qq[
SNP:        $snp_id
ACCESSION:  $supercontig
LOCATION:   $start_pos
SAMPLESIZE: 2
LENGTH:     1
5'_FLANK:   $five_prime_seq
OBSERVED:   $reference_base
3'_FLANK:   $three_prime_seq
COMMENT:
quality: $phred_quality
read_depth: $read_depth
||
];
	}
	close( $ofh );
	close( $sfh );
}
