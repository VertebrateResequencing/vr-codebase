package QC_454_ssaha;
use strict;
use Carp;
use File::Basename;
use File::Spec;
use Cwd;

use G1KUtilities;
use QC;

my $G1K = $ENV{ 'G1K' };
my $MOUSE_REF_FA = $G1K.'/MOUSE/ref/NCBIM37_um.fa';
my $HUMAN_FEMALE_REF_FA = $G1K.'/ref/human_b36_female.fa';
my $HUMAN_MALE_REF_FA = $G1K.'/ref/human_b36_male.fa';

my $HUMAN_FEMALE_REF_FAI = $G1K.'/ref/human_b36_female.fa.fai';
my $HUMAN_MALE_REF_FAI = $G1K.'/ref/human_b36_male.fa.fai';

my $GENDERS_FILE = '/lustre/sf4/1kgenomes/ref/genders.txt';
my $REMOTE_REF_DIR= $G1K.'/ref';

my $SSAHA2 = $G1K."/bin/ssaha2";
my $SAMTOOLS = $G1K."/bin/samtools";
=pod

=head1 NAME

Mapping_454_ssaha

=head1 SYNOPSIS

=head1 REQUIRES

Perl5.8.8

=head1 DESCRIPTION

=head1 METHODS

=head2 mappingInProgress

	Arg [1]    : lane directory
	Returntype : 0/1
=cut

