#
# Author:   Petr Danecek (pd3@sanger.ac.uk)    Team 145
# Author:	John Maslen  (jm23@sanger.ac.uk)   Team 145
#
#--------------- QuerySNPsData ---------------------------------------
#
# Takes the cached data and creates a csv-delimited output file.
#

package SNPs::QuerySNPsDataCSV;

use strict;
use warnings;
use base qw(SNPs::QuerySNPs);
use POSIX qw(strftime);
use CGI::Carp qw(fatalsToBrowser);

sub new
{
    my ($class, $args) = @_;

    my $self = $class->SUPER::new($args);

    if ( !$self->cache_exists() ) { die "Error: no data in cache??\n"; }
    my $date = strftime "%Y-%m-%d", localtime;
	my $str_count =  scalar keys %{$$self{selected_strains}};
	my $loc = '['.$$self{chrm}.':'.$$self{from}.'-'.$$self{to}.']';
	my $file = 'SNPs'.$str_count."_mouse_strains_".$loc."_".$date.".csv";
    $$self{writer}->fname($file);
    return $self;
}


sub run
{
    my ($self) = @_;
	$self->print_header();
    while (my $pos=$self->cache_get_next())
    {
        $self->print_row($pos);
    }
    $self->print_footer();
    return;
}

sub print_header
{
    my ($self) = @_;
    my $strains = $$self{selected_strains};
    my $html = $$self{writer};
    if ( $$self{display_dload_params} )
    {
        $html->out($$self{display_dload_params});
    }
    $html->out("#\tBase calling key (first column for each strain):
#\t\t'A C G T'\t= High confidence SNP
#\t\t'-' (hyphen)\t= High confidence reference
#\t\t'a c t g'\t= Low confidence SNP
#\t\t'~' (tilde)\t= Low confidence reference
#\t\t'.' (period)\t= Genotype not called\n");
    $html->out("Gene,Chromosome,Position,Reference");
    for my $str (sort {$$strains{$a}<=>$$strains{$b}} keys %$strains)
    {
        $html->out(",$str,Consequence");
    }
    $html->out("\n");

}

sub print_footer
{
    my ($self) = @_;
    return;
}

sub cache_get_next
{
    my ($self) = @_;

    if ( !$$self{cached_data} ) { return 0; }

    my $data = $$self{cached_data};
    if ( $$self{icache} >= scalar @$data ) { return 0; }

    return $$data[$$self{icache}++];
}

sub print_row
{
    my ($self,$row) = @_;

    my $html = $$self{'writer'};
    my $session = $$self{'session'};

    my ($pos,$chr,$base,$gene_name) = $self->nonzero_column_data($row);

    $html->out("$gene_name,$chr,$pos,$base");

    my $ncols = scalar keys %{$$self{'selected_strains'}};
    for (my $i=0; $i<$ncols; $i++)
    {
        my $conseqs = {};
        my $atg_qual = $$row[$i]->{'atg_qual'};
        my $strain_out = '-';
                
        if ($atg_qual == -1) {
        	$strain_out = '~';
        }
		if ($atg_qual == -8) {
        	$strain_out = '.';
        }
        for my $type (@{$$row[$i]->{'consequence'}})
        {
            if ( !$type || $type eq 'SPLICE_SITE' ) { next }  # ignore these - according to Dave these are rubbish
            $$conseqs{$type} = 1;
        }

        if ( $$row[$i]->{'ref_base'} )
        {	
        	$strain_out = $$row[$i]->{'snp_base'};
        	if ($atg_qual != 1) {
            	$strain_out = lc($strain_out);
            } 
            $html->out(",$strain_out");
            $html->out("," . join(';',sort keys %$conseqs));
        }
        else
        {
            $html->out(",$strain_out,-");
        }
    }
    $html->out("\n");
}

1;

