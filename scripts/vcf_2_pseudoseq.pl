#!/Users/fy2/perl5/perlbrew/perls/perl-5.8.8/bin/perl

use strict;
use warnings;
use Getopt::Declare;
use Log::Log4perl qw(get_logger);
use Pathogens::Variant::Utils::PseudoReferenceMaker;


our $VERSION = 0.01;








#log4perl's log configuration setting
#Activates the logger at level "INFO", and prints the errors to SCREEN
my $logconf = q(

    #set the log category to root (i.e. Pathogens::*) and
    #set the level to INFO and higher, choose 'Screen' as appender
    log4perl.category = INFO, Screen
    
    
    #create the appender, turn off stderr (ie choose stdout) create log pattern
    log4perl.appender.Screen = Log::Log4perl::Appender::Screen
    log4perl.appender.Screen.stderr  = 0
    log4perl.appender.Screen.layout = PatternLayout
    log4perl.appender.Screen.layout.ConversionPattern=[%p][%d] %m (%C %L)%n
);
Log::Log4perl::init( \$logconf );
my $logger = get_logger();





#set the default argument values
our %args = (
     vcf_file => ''
   , out => ''
   , bam => ''
   , reference => ''
   , lane_name => ''
   , depth        => 4
   , depth_strand => 2
   , ratio        => 0.8
   , quality      => 50
   , map_quality  => 0
   , af1            => 0.95
   , af1_complement => 0.05
   , ci95           => 0.0
   , strand_bias    => 0.001
   , base_quality_bias => 0.0
   , map_bias     => 0.001
   , tail_bias    => 0.001
);

my $specification = q(

[strict]

==========================================================================================
Required parameters:

    -v[cf] <string>                    	Variation input file in vcf format [required]
                                            { 
                                            	$main::args{vcf_file} = $string; 
                                            }
    -o[ut] <string>                    	Output file name for the pseudoreference [required]
                                            { 
                                            	$main::args{out} = $string; 
                                            }
    -b[am] <string>                    	Bam file with a valid header and read alignments [required]
                                            { 
                                            	$main::args{bam} = $string; 
                                            }
    -l[ane_name] <string>              	Lane id (in 'XXXX_X#XX' format) [required]
                                            { 
                                            	$main::args{lane_name} = $string; 
                                            }
    -R[eference] <string>              	Reference sequence in fasta format (same fasta as used in the mapping) [required]
                                            { 
                                            	$main::args{reference} = $string; 
                                            }
------------------------------------------------------------------------------------------

Optional parameters:

    -u[sage]                           	Show usage information and exit
                                            {
                                            	$self->usage(0);
                                            }

    -d[epth] <integer:+i>              	(4) Minimum number of reads matching SNP 
                                            { 
                                            	$main::args{depth} = $integer; 
                                            }

    -D[epth_strand] <integer:+i>       	(2) Minimum number of reads matching SNP per strand
                                            { 
                                            	$main::args{depth_strand} = $integer; 
                                            }
                                            
    -r[atio] <float:+n>                	(0.8) Minimum ratio of first to second base call
                                            {
                                            	reject($float < 0.5 || $float > 1 => "Ratio of first to second base (-r) must be >= 0.5 and <1");
                                            	$main::args{ratio} = $float;
                                            }
                                            
    -q[uality] <integer:0+i>           	(50) Minimum variant quality (Vcf QUAL field)
                                            {
                                            	reject($integer > 99 => "Base quality (-q) must be between 0 and 99");
                                            	$main::args{quality} = $integer;
                                            }
                                            
    -m[ap_quality] <integer:0+n>       	(0) Minimum mapping quality
                                            {
                                            	reject($integer > 99 => "Mapping quality (-m) must be between 0 and 99");
                                            	$main::args{map_quality} = $integer;
                                            }
                                            
    -A[F1] <float:0+n>                 	(0.95) Minimum allele frequency (you would expect an AF of 1 for haploid SNPs). For non-SNP bases, the program will use 1- this number
                                            {
                                            	reject($float > 1 => "Minimum allele frequency for SNPs (-a) must be between 0 and 1");
                                            	$main::args{af1} = $float;
                                            }
                                            
    -C[I95] <float:0+n>                	(0.0) Maximum 95% confidence interval variation from AF
                                            {
                                            	reject($float > 1 => "Maximum 95% confidence interval of allele frequency (-c) must be between 0 and 1");
                                            	$main::args{ci95} = $float;
                                            }
                                            
    -S[trand_bias] <float:0+n>         	(0.001) P-value cutoff for strand bias
                                            {
                                            	reject($float > 1 => "p-value cutoff for strand bias (-S) must be between 0 and 1");
                                            	$main::args{strand_bias} = $float;
                                            }
                                            
    -Q[uality_bias] <float:0+n>        	(0.0) P-value cutoff for base quality bias
                                            {
                                            	reject($float > 1 => "p-value cutoff for base quality bias (-Q) must be between 0 and 1");
                                            	$main::args{base_quality_bias} = $float;
                                            }
                                            
    -M[ap_bias] <float:0+n>            	(0.001) P-value cutoff for mapping bias
                                            {
                                            	reject($float > 1 => "p-value cutoff for mapping bias (-M) must be between 0 and 1");
                                            	$main::args{map_bias} = $float;
                                            }
                                            
    -t[ail_bias] <float:0+n>           	(0.001) P-value cutoff for tail distance bias
                                            {
                                            	reject($float > 1 => "p-value cutoff for tail distance bias (-T) must be between 0 and 1");
                                            	$main::args{tail_bias} = $float;
                                            }
------------------------------------------------------------------------------------------
);

$logger->logdie("Exiting: Argument errors") 
	unless ( Getopt::Declare->new($specification) );

my $pseudo_maker = Pathogens::Variant::Utils::PseudoReferenceMaker->new(arguments => \%args);

$pseudo_maker->create_pseudo_reference();
print $pseudo_maker->get_statistics_dump();

