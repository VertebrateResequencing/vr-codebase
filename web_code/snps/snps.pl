#!/usr/local/bin/perl -Tw
#
# Author:    Petr Danecek (pd3@sanger.ac.uk)    Team 145
# Modified by : John Maslen (jm23@sanger.ac.uk)  Team 145
#
#   1) No params        .. 	Show the form. Additional parameters to this are:
#                           error display this error message
#   2) action           .. 	Run the SNPs (action='snps'),  indels (action='indels') or SVs (action='structvar') #							query, or download cached data in tab or csv format (possible actions = 			#							'snps_dload', 'indels_dload' or 'svs_dload' to download tab files or 				#							'snps_dload_csv', 'indels_dload_csv' or 'structvar_dload_csv' for csv)
#                           
#   Additional parameters for this are: 
#                               cache            ..access cached results (mysql query will not be run again)
#                               cnsq_*           ..consequence types to display (see keys %{$$params{conseqs}})
#                               location         ..filter SNPs and indels from the region or the gene
#                               page             ..display this page from cache (indexed from 0)
#                               str_*            ..set of strains (see keys %{$$params{mouseinfo}})
#								sv_*			 ..set of structural variations (see keys %{$$params{sv_types}})

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
use SNPs::QuerySNPsData;
use SNPs::QuerySNPsDataCSV;
use SNPs::QueryIndels;
use SNPs::QueryIndelsData;
use SNPs::QueryIndelsDataCSV;
use SNPs::QuerySVs;
use SNPs::QuerySVsData;
use SNPs::QuerySVsDataCSV;
use SNPs::Session;

$ENV{PATH}= '/usr/local/bin';   # solely to stop taint from barfing

my $params = 
{
    'debug'        => 0,
    'available'    => 1,
    'myself'       => 'snps.pl',
    'lookseq'      => '/cgi-bin/modelorgs/mousegenomes/lookseq/index.pl',
    'lseq_win'     => 100,
    'lseq_width'   => 700,
    'lseq_display' => 'display=|perfect|single|inversions|pairlinks|potsnps|uniqueness|gc|coverage|',
    'lseq_display_sv' => 'display=|perfect|single|inversions|pairlinks|potsnps|uniqueness|gc|coverage|orientation|',
    'imgs'         => '/modelorgs/mousegenomes/snps-gfx',
    'css'          => ['/modelorgs/mousegenomes/snps.css'],
    'jsfile'       => ['http://js.sanger.ac.uk/jquery-1.3.2.min.js','/modelorgs/mousegenomes/snps.js'],

    'store_hours' => 1,     # The stored data are few bytes big only so far
    'store_key'   => 'modelorgs/mousegenomes/',
    
    'tabix_location' => '/WWW/SANGER_docs/htdocs/modelorgs/mousegenomes/lookseq/data/tabix',
    
    'index_directory' => '/WWW/SANGER_docs/htdocs/modelorgs/mousegenomes/lookseq/data',

    'snps_location' => 'ftp://ftp-mouse.sanger.ac.uk/web_data/web_snps.vcf.gz',
    
    'indels_location' => 'ftp://ftp-mouse.sanger.ac.uk/web_data/web_indels.vcf.gz',
	
	'svs_location' => 'ftp://ftp-mouse.sanger.ac.uk/web_data/web_svs.sdp.gz',
 	
 	'gene_db_location' => '/WWW/SANGER_docs/htdocs/modelorgs/mousegenomes/lookseq/data/snps-genes.sqlite',
	
	'vcf_strain_order' => ['129P2/OlaHsd','129S1/SvImJ','129S5SvEvBrd','A/J','AKR/J','BALB/cJ','C3H/HeJ','C57BL/6NJ','CAST/EiJ','CBA/J','DBA/2J','LP/J','NOD/ShiLtJ','NZO/HlLtJ','PWK/PhJ','SPRET/EiJ','WSB/EiJ'],
    
    'vcf_strain_map' => {
    	'129P2' => { name => '129P2/OlaHsd', },
    	'129S1' => { name => '129S1/SvImJ',},
    	'129S5' => { name => '129S5SvEvBrd',},
    	'AKR'  => { name => 'AKR/J',},
    	'A_J' => { name => 'A/J',},
    	'BALB'  => { name => 'BALB/cJ',},
    	'C3H'  => { name => 'C3H/HeJ',},
    	'C57BL'  => { name => 'C57BL/6NJ',},
    	'CAST'  => { name => 'CAST/EiJ',},
    	'CBA'  => { name => 'CBA/J',},
    	'DBA'  => { name => 'DBA/2J',},
    	'LP_J'  => { name => 'LP/J',},
    	'NOD'  => { name => 'NOD/ShiLtJ',},
    	'NZO'  => { name => 'NZO/HlLtJ',},
    	'PWK'  => { name => 'PWK/PhJ',},
    	'SPRET'  => { name => 'SPRET/EiJ',},
    	'WSB'  => { name => 'WSB/EiJ',},
    },
 
    'mouseinfo' => {
        '129P2/OlaHsd'       => { 	display_name => '129P2/OlaHsd',
        					caption  => 'Commonly used to make embryonic stem cell lines.', 
                            img=>'mouse-129p2.png',
                            bam=>'129P2_OlaHsd.bam', },
        '129S1/SvImJ' => {  display_name => '129S1/SvImJ',
        					caption  => 'Commonly used to make embryonic stem cell lines. 
                                              Progenitor strain of the collaborative cross.',
                            img=>'mouse-129S1_SvImJ.png',
                            bam=>'129S1_SvImJ.bam' },
        '129S5SvEvBrd'       => { display_name => '129S5SvEvBrd',
        					caption  => 'Commonly used to make embryonic stem cell lines.', 
                            img=>'mouse-129s5.png',
                            bam=>'129S5SvEvBrd.bam' },
        'A/J'         => { display_name => 'A/J',
        					caption  => 'An asthma model. Progenitor strain of the collaborative cross and of the heterogeneous stock cross.', 
                            img=>'mouse-A_J.png', 
                            bam=>'A_J.bam' },
        'AKR/J'       => { display_name => 'AKR/J',
        					caption  => 'High leukemie incidence. Hyporesponsive to diets containing high levels of fat and cholesterol,
                                          and resistant to aortic lesion formation. Progenitor strain of the heterogeneous stock cross.', 
                            img=>'mouse-AKR_J.png', 
                            bam=>'AKR_J.bam' },
        'BALB/cJ'     => { display_name => 'BALB/cJ',
        					caption  => 'Prone to develop mammary and kidney cancer. Progenitor strain of the 
                                          collaborative cross and of the heterogeneous stock cross.', 
                            img=>'mouse-BALB_cByJ.png',
                            bam=>'BALB_cJ.bam' },
        'C3H/HeJ'     => { display_name => 'C3H/HeJ',
        					caption  => 'Spontaneously develops mammary tumours. Highly susceptible to Gram-negative 
                                          bacterial infections. Progenitor strain of the heterogeneous stock cross.', 
                            img=>'mouse-C3H_HeJ.png', 
                            bam=>'C3H_HeJ.bam', },
        'C57BL/6NJ'   => { display_name => 'C57BL/6NJ',
        					caption  => 'Used in KOMP and EUCOMM programmes to knockout every gene in the mouse genome.', 
                            img=>'mouse-C57BL_6N.png', 
                            bam=>'C57BL_6NJ.bam', },
        'CAST/EiJ'    => { display_name => 'CAST/EiJ',
        					caption  => 'Resistant to cancer and infections.', 
                            img=>'mouse-CAST_EiJ.png',
                            bam=>'CAST_EiJ.bam', },
        'CBA/J'       => { display_name => 'CBA/J',
        					caption  => 'Renal tubulointerstitial lesions observed at a high frequency. Prone to exocrine 
                                          pancreatic insufficiency syndrome. Progenitor strain of the heterogeneous stock cross.', 
                            img=>'mouse-CBA_J.png',
                            bam=>'CBA_J.bam', },
        'DBA/2J'      => { display_name => 'DBA/2J',
        					caption  => 'Develops agressive early hearing loss. Extreme intolerance to alcohol and morphine. 
                                          Aging DBA/2J mice develop progressive eye abnormalities.', 
                           img=>'mouse-DBA_2J.png', 
                           bam=>'DBA_2J.bam' },
        'LP/J'        => { display_name => 'LP/J',
        					caption  => 'High susceptibility to audiogenic seizures. This strain is also reported to have a 
                                          fairly high incidence of tumors that develop later in life. Progenitor strain 
                                          of the heterogeneous stock cross.', 
                            img=>'mouse-LP_J.png', 
                            bam=>'LP_J.bam', },
        'NOD/ShiLtJ'  => { display_name => 'NOD/ShiLtJ',
        					caption  => 'This strain is a polygenic model for type 1 (non-obese) diabetes. Progenitor strain of the collaborative cross.', 
                            img=>'mouse-NOD_ShiLtJ.png', 
                            bam=>'NOD_ShiLtJ.bam', },
        'NZO/HlLtJ'   => { display_name => 'NZO/HlLtJ',
        					caption  => 'New Zealand Obese. Susceptible to type II diabetes. Progenitor strain of the collaborative cross.', 
                            img=>'mouse-nzo.png',
                            bam=>'NZO_HlLtJ.bam', },
        'PWK/PhJ'     => { display_name => 'PWK/PhJ',
        					caption  => 'Susceptibility to type I diabetes and various behavioral traits. Progenitor strain of the collaborative cross.', 
                            img=>'mouse-pwk.png',
                            bam=>'PWK_PhJ.bam' },
        'SPRET/EiJ' => { display_name => 'SPRET/EiJ',
        					caption  => 'Resistant to cancer and infections.', 
                            img=>'mouse-Spretus_EiJ.png',
                            bam=>'SPRET_EiJ.bam', },
        'WSB/EiJ'     => { display_name => 'WSB/EiJ',
        					caption  => 'Displays extremely long life-span. Progenitor strain of the collaborative cross.', 
                            img=>'mouse-WSB_EiJ.png',
                            bam=>'WSB_EiJ.bam', },
    },

    'conseqs' => { 
        '3PRIME_UTR'            => { label=>"3' UTR", style=>'c1' },
        '5PRIME_UTR'            => { label=>"5' UTR", style=>'c2' },
        'COMPLEX_INDEL'			=> { label=>'Complex indel', style=>'c12' },
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
        'WITHIN_NON_CODING_GENE'=> { label=>'Non-coding gene', style=>'c5' },
        'WITHIN_MATURE_miRNA'   => { label=>'Mature miRNA', style=>'c7' },
    },
    
    'sv_types' => {
    	'Deletion'		=> { front_label => 'Deletion', label => 'Deletion (D)', style=>'c5'},
    	'Insertion'		=> { front_label => 'Insertion', label => 'Insertion (I)', style=>'c6' },
    	'Copy_Number_Gain'		=> { front_label => 'Copy number gains', label => 'Copy number gains (G)', style=>'c1' },
    	'Inversion'		=> { front_label => 'Inversion', label => 'Inversion (V)', style=>'c7' },
    	'Complex_events'		=> { front_label => 'Complex events', label => 'Complex events', style=>'c11' },
    },

    'sv_type_map' => {
        '^DEL$' => {display => 'DELETION', label =>'D', type => 'Deletion'},
        '^DELLINKED$' => {display => 'DELETION', label =>'D', type => 'Deletion'},
        '^DELINS$' => {display => 'DELETION+INSERTION', label =>'DI', type => 'Complex_events'},
        '^INS$' => {display => 'INSERTION', label =>'I', type => 'Insertion'},
        '^INS\|\w+' => {display => 'TE INSERTION', label =>'I', type => 'Insertion'},
        '^INSLINKED$' => {display => 'INSERTION', label =>'I', type => 'Insertion'},
        '^INSLINKED\(INV\)$' => {display => 'INSERTION', label =>'I', type => 'Insertion'},
        '^INSLINKED\(INV\)DEL$' => {display => 'DELETION+INSERTION', label =>'DI', type => 'Complex_events'},
        '^INSLINKED\(INV\)INS\|\w+' => {display => 'TE INSERTION', label =>'I', type => 'Insertion'},
        '^GAIN$' => {display => 'COPY NUMBER GAIN', label =>'G', type => 'Copy_Number_Gain'},
        '^TANDEMDUP$' => {display => 'TANDEM DUPLICATION', label =>'G', type => 'Copy_Number_Gain'},
        '^TANDEMLOWDUP$' => {display => 'TANDEM DUPLICATION', label =>'G', type => 'Copy_Number_Gain'},
        '^TANDEMDUPINV$' => {display => 'INVERTED TANDEM DUPLICATION', label =>'G', type => 'Copy_Number_Gain'},
        '^TANDEMLOWDUPINV$' => {display => 'INVERTED TANDEM DUPLICATION', label =>'G', type => 'Copy_Number_Gain'},
        '^INV$' => {display => 'INVERSION', label =>'V', type => 'Inversion'},
        '^INVDEL$' => {display => 'INVERSION+DELETION', label =>'VD', type => 'Complex_events'},
        '^INVINS$' => {display => 'INVERSION+INSERTION', label =>'VI', type => 'Complex_events'},
        '^INVINS\|\w+' => {display => 'INVERSION+TE INSERTION', label =>'VI', type => 'Complex_events'},
        '^DELINS\|\w+' => {display => 'DELETION+TE INSERTION', label =>'DI', type => 'Complex_events'},
        '^DELLINKEDINS\|\w+' => {display => 'DELETION+TE INSERTION', label =>'DI', type => 'Complex_events'},
        '^INSLINKEDINS\|\w+' => {display => 'TE INSERTION', label =>'I', type => 'Insertion'},
        '^INVDELINS\|\w+' => {display => 'INVERSION+DELETION+INSERTION', label =>'DVI', type => 'Complex_events'},
    },  
    
    'sv_combined_events_key' => {
		'DI' => {label =>'Deletion with small insertion'},
		'VD' => {label => 'Deletion contained within an inversion'},
		'VI' => {label => 'Insertion adjacent to an inversion'},
		'DVI' => {label => 'Inversion adjacent to a deletion and insertion'},		
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
    elsif ( $action eq 'structvar' )
    {
        eval 
        { 
            my $query = SNPs::QuerySVs->new({ writer=>$html, %$params }); 
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
    elsif ( $action eq 'snps_dload_csv' ) 
    {
        eval 
        {
            my $writer = SNPs::DataWriter->new($params);
            my $query  = SNPs::QuerySNPsDataCSV->new({ writer=>$writer, %$params }); 
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
    elsif ( $action eq 'indels_dload_csv' ) 
    {
        eval 
        {
            my $writer = SNPs::DataWriter->new($params);
            my $query  = SNPs::QueryIndelsDataCSV->new({ writer=>$writer, %$params }); 
            $query->run();
            1;
        } or $eval_ok=0;
    }    
    elsif ( $action eq 'structvar_dload' ) 
    {
        eval 
        {
            my $writer = SNPs::DataWriter->new($params);
            my $query  = SNPs::QuerySVsData->new({ writer=>$writer, %$params }); 
            $query->run();
            1;
        } or $eval_ok=0;
    }  
    elsif ( $action eq 'structvar_dload_csv' ) 
    {
        eval 
        {
            my $writer = SNPs::DataWriter->new($params);
            my $query  = SNPs::QuerySVsDataCSV->new({ writer=>$writer, %$params }); 
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


