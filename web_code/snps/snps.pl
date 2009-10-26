#!/usr/local/bin/perl -Tw
#
# Author:    Petr Danecek (pd3@sanger.ac.uk)    Team 145
#
#   1) No params        .. Show the form. Additional parameters to this are:
#                               error           .. display this error message
#   2) action           .. Run the SNPs (action='snps'),  indels (action='indels') query, or
#                           download cached data (action='snps_dload' or action='indels_dload')
#                           Additional parameters for this are: 
#                               cache            .. access cached results (mysql query will not be run again)
#                               cons_qual        .. consensus quality
#                               cnsq_*           .. consequence types to display (see keys %{$$params{conseqs}})
#                               depth_max        .. filter SNPs with depth less or equal than
#                               depth_min        .. filter SNPs with depth greater or equal than
#                               location         .. filter SNPs and indels from the region or the gene
#                               map_qual         .. mapping quality
#                               page             .. display this page from cache (indexed from 0)
#                               poly             .. show only SNPs polymorphic between selected strains
#                               snp_qual         .. SNP quality
#                               str_*            .. set of strains (see keys %{$$params{mouseinfo}})
#   3) ajax             .. lookup values in the db (additional prameters for this are: gene)
#
#

use strict;
use warnings;

use SangerPaths qw(core mousegenomes);
use SangerWeb;
use DBI;
use POSIX ":sys_wait_h";

use SNPs::Writer;
use SNPs::DataWriter;
use SNPs::DebugWriter;
use SNPs::QuerySNPs;
use SNPs::QueryIndels;
use SNPs::QuerySNPsData;
use SNPs::QueryIndelsData;
use SNPs::Session;

$ENV{PATH}= '/usr/local/bin';   # solely to stop taint from barfing

my $params = 
{
    'debug'        => 0,
    'myself'       => 'snps.pl',
    'lookseq'      => '/cgi-bin/modelorgs/mousegenomes/lookseq/index.pl',
    'lseq_win'     => 100,
    'lseq_width'   => 700,
    'lseq_display' => 'display=|perfect|snps|single|inversions|pairlinks|potsnps|uniqueness|gc|coverage',
    'imgs'         => '/modelorgs/mousegenomes/snps-gfx',
    'css'          => ['/modelorgs/mousegenomes/snps.css'],
    'jsfile'       => ['http://js.sanger.ac.uk/jquery-1.3.2.min.js','/modelorgs/mousegenomes/snps.js'],

    'store_hours' => 1,     # The stored data are few bytes big only so far
    'store_key'   => 'modelorgs/mousegenomes/',

    # 'db_info' => 'DBI:mysql:host=mcs4a;port=3307;database=mouse_web_snps',    # development version
    'db_info' => 'DBI:mysql:host=webdbsrv6;port=3306;database=mouse_web_snps',
    'db_user' => 'vreseq_ro', 
    'db_pass' => '',

    'mouseinfo' => {
        '129P2'       => { caption  => 'Commonly used to make embryonic stem cell lines.', 
                            bam=>'129P2.bam', },
        '129S1/SvImJ' => { caption  => 'Commonly used to make embryonic stem cell lines. 
                                              Progenitor strain of the collaborative cross.',
                            img=>'mouse-129S1_SvImJ.png',
                            bam=>'129S1_SvImJ.bam' },
#        '129S5'       => { caption  => 'Commonly used to make embryonic stem cell lines.', 
#                            bam=>'' },
        'A/J'         => { caption  => 'An asthma model. Progenitor strain of the collaborative cross and of the heterogeneous stock cross.', 
                            img=>'mouse-A_J.png', 
                            bam=>'A_J.bam' },
        'AKR/J'       => { caption  => 'High leukemie incidence. Hyporesponsive to diets containing high levels of fat and cholesterol,
                                          and resistant to aortic lesion formation. Progenitor strain of the heterogeneous stock cross.', 
                            img=>'mouse-AKR_J.png', 
                            bam=>'AKR_J.bam' },
        'BALB/cJ'     => { caption  => 'Prone to develop mammary and kidney cancer. Progenitor strain of the 
                                          collaborative cross and of the heterogeneous stock cross.', 
                            img=>'mouse-BALB_cByJ.png', 
                            bam=>'BALBc_J.bam' },
        'C3H/HeJ'     => { caption  => 'Spontaneously develops mammary tumours. Highly susceptible to Gram-negative 
                                          bacterial infections. Progenitor strain of the heterogeneous stock cross.', 
                            img=>'mouse-C3H_HeJ.png', 
                            bam=>'C3H_HeJ.bam', },
        'C57BL/6NJ'   => { caption  => 'Used in KOMP and EUCOMM programmes to knockout every gene in the mouse genome.', 
                            img=>'mouse-C57BL_6N.png', 
                            bam=>'C57BL_6N.bam', },
        'CAST/EiJ'    => { caption  => 'Resistant to cancer and infections.', 
                            img=>'mouse-CAST_EiJ.png',
                            bam=>'CAST_Ei.bam', },
        'CBA/J'       => { caption  => 'Renal tubulointerstitial lesions observed at a high frequency. Prone to exocrine 
                                          pancreatic insufficiency syndrome. Progenitor strain of the heterogeneous stock cross.', 
                            img=>'mouse-CBA_J.png',
                            bam=>'CBA_J.bam', },
        'DBA/2J'      => { caption  => 'Develops agressive early hearing loss. Extreme intolerance to alcohol and morphine. 
                                          Aging DBA/2J mice develop progressive eye abnormalities.', 
                           img=>'mouse-DBA_2J.png', 
                           bam=>'DBA_2J.bam' },
        'LP/J'        => { caption  => 'High susceptibility to audiogenic seizures. This strain is also reported to have a 
                                          fairly high incidence of tumors that develop later in life. Progenitor strain 
                                          of the heterogeneous stock cross.', 
                            img=>'mouse-LP_J.png', 
                            bam=>'LP_J.bam', },
        'NOD/ShiLtJ'  => { caption  => 'This strain is a polygenic model for type 1 (non-obese) diabetes. Progenitor strain of the collaborative cross.', 
                            img=>'mouse-NOD_ShiLtJ.png', 
                            bam=>'NOD.bam', },
        'NZO/HiLtJ'   => { caption  => 'New Zealand Obese. Susceptible to type II diabetes. Progenitor strain of the collaborative cross.', 
                            bam=>'NZO.bam', },
        'PWK/PhJ'     => { caption  => 'Susceptibility to type I diabetes and various behavioral traits. Progenitor strain of the collaborative cross.', 
                            bam=>'PWK_Ph.bam' },
        'Spretus/EiJ' => { caption  => 'Resistant to cancer and infections.', 
                            img=>'mouse-Spretus_EiJ.png',
                            bam=>'Spretus_Ei.bam', },
        'WSB/EiJ'     => { caption  => 'Displays extremely long life-span. Progenitor strain of the collaborative cross.', 
                            img=>'mouse-WSB_EiJ.png',
                            bam=>'WSB_Ei.bam', },
    },

    'conseqs' => { 
        '3PRIME_UTR'            => { label=>"3' UTR", style=>'c1' },  
        '5PRIME_UTR'            => { label=>"5' UTR", style=>'c2' },
        'ESSENTIAL_SPLICE_SITE' => { label=>'Splice site', style=>'c3' }, 
        'INTERGENIC'            => { label=>'Intergenic', style=>'c4' }, 
        'INTRONIC'              => { label=>'Intronic', style=>'c5' },
        'NON_SYNONYMOUS_CODING' => { label=>'Non-synonymous coding', style=>'c6' }, 
        'SYNONYMOUS_CODING'     => { label=>'Synonymous coding', style=>'c7' }, 
        # This is probably a mistake - the database does not contain a single consequence of this type.
        # 'NO_CODING_TRANSCRIPTS' => { label=>'Non-coding transcripts', style=>'c8' }, 
        'STOP_GAINED'           => { label=>'Stop gained', style=>'c9' }, 
        'STOP_LOST'             => { label=>'Stop lost', style=>'c10' },
        'FRAME_SHIFT'           => { label=>'Frame shift', style=>'c11' },
    },
};

my $html   = $$params{'debug'} ? SNPs::DebugWriter->new($params) : SNPs::Writer->new($params);
my $action = $html->param('action');

if ( $action )
{
    #if ( !$html->param('debug') ) { $html->error_exit("Work in progress, please come back later."); }

    my $eval_ok=1;
    if ( $action eq 'snps' )
    {
        eval 
        { 
            my $query = SNPs::QuerySNPs->new({ writer=>$html, %$params }); 
            $query->run();
            1;
        } or $eval_ok=0;
    }
    elsif ( $action eq 'indels' )
    {
        eval 
        { 
            my $query = SNPs::QueryIndels->new({ writer=>$html, %$params }); 
            $query->run();
            1;
        } or $eval_ok=0;
    }
    elsif ( $action eq 'snps_dload' ) 
    {
        eval 
        {
            my $writer = SNPs::DataWriter->new($params);
            my $query  = SNPs::QuerySNPsData->new({ writer=>$writer, %$params }); 
            $query->run();
            1;
        } or $eval_ok=0;
    }
    elsif ( $action eq 'indels_dload' ) 
    {
        eval 
        {
            my $writer = SNPs::DataWriter->new($params);
            my $query  = SNPs::QueryIndelsData->new({ writer=>$writer, %$params }); 
            $query->run();
            1;
        } or $eval_ok=0;
    }
    if ( !$eval_ok ) { $html->error_exit($@); return; }
}
else
{
    print_form($html);
}

exit;

#--------------------------------------------------------------------------

sub print_form
{
    my ($html) = @_;
    $html->print_error();
    $html->print_form();
    $html->print_footer();
    return;
}


