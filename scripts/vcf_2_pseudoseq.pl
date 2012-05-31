#!/Users/fy2/perl5/perlbrew/perls/perl-5.8.8/bin/perl

use strict;
use warnings;
use Getopt::Declare;
use Data::Dumper;
use Pathogens::Variant::Iterator::Vcf;
use Pathogens::Variant::Evaluator::Pseudosequence;
use Log::Log4perl;


my $conf = q(
    log4perl.category.Pathogens.Variant.Evaluator.Pseudosequence = FATAL, Screen
    log4perl.appender.Screen = Log::Log4perl::Appender::Screen
    log4perl.appender.Screen.layout = Log::Log4perl::Layout::SimpleLayout
);

Log::Log4perl::init( \$conf );
 
our $VERSION = 0.01;


#set the default values for commandline arguments
our %args = (
                 vcf_file     => ''
               , out_prefix   => ''
               , depth        => 4
               , depth_strand => 2
               , ratio        => 0.8
               , quality      => 50
               , map_quality  => 0
               , af1          => 0.95
               , ci95         => 0.0
               , strand_bias  => 0.001
               , base_quality_bias => 0.0
               , map_bias     => 0.001
               , tail_bias    => 0.001
);

my $specification = q(

[strict]

==========================================================================================
Required parameters:

    -v[cf_file] <string>               	Variation input file in vcf format [required]
                                            { 
                                            	$main::args{vcf_file} = $string; 
                                            }
    -o[ut_prefix] <string>             	Prefix for output file names [required]
                                            { 
                                            	$main::args{out_prefix} = $string; 
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

die "Exiting... Argument errors" 
	unless (Getopt::Declare->new($specification));

#a little more argument sanity check
my $doubled_arg_D_depth_strand = $args{depth_strand} * 2;
if ( $args{depth} < $doubled_arg_D_depth_strand ) {
    print "'-d' (depth) must be greater than '-D' (depth_strand)! Silently increasing it to " .$doubled_arg_D_depth_strand ."\n";
    $args{depth} = $doubled_arg_D_depth_strand;
}

#af1 is used to check quality of variant sites.
#we will use 1-af1 (af1's complement) to check the quality if non-variant sites
my $af1_complement = 1 - $args{af1};



#this will traverse the VCF file and return Pathogens::Variant::Event::* objects
my $iterator = Pathogens::Variant::Iterator::Vcf->new(vcf_file => $args{vcf_file});

#create an evaluator object to filter bad calls based on user criteria
my $evaluator = Pathogens::Variant::Evaluator::Pseudosequence->new(
      minimum_depth => $args{depth}
    , minimum_depth_strand => $args{depth_strand}
    , minimum_ratio => $args{ratio}
    , minimum_quality => $args{quality}
    , minimum_map_quality => $args{map_quality}
    , minimum_af1 => $args{af1}
    , af1_complement =>  $af1_complement
    , minimum_ci95 => $args{ci95}
    , minimum_strand_bias => $args{strand_bias}
    , minimum_base_quality_bias => $args{base_quality_bias}
    , minimum_map_bias => $args{map_bias}
    , minimum_tail_bias => $args{tail_bias}
);

my $c = 0;
while( my $event = $iterator->next_event() ) {
    
    $c++; print $c, "\n" unless $c % 100000;
    
    $evaluator->evaluate($event); #note: evaluator will modify heterozygous calls

    if ($event->passed_evaluation) {
        
    }
}

$iterator->close_vcf();

#print join("\n", $iterator-> get_metalines);
print $evaluator->dump_evaluation_statistics;
