#
# Author:		John Maslen  (jm23@sanger.ac.uk)   Team 145
# Author:    	Petr Danecek (pd3@sanger.ac.uk)    Team 145
#
#--------------- QueryIndelsDataCSV ---------------------------------
#
# Takes the cached data and creates a cvs-delimited output file. Same as QuerySNPsDataCSV, but for Indels
#

package SNPs::QueryIndelsDataCSV;

use strict;
use warnings;

use POSIX qw(strftime);
use base qw(SNPs::QuerySNPsDataCSV);

sub new
{
    my ($class, $args) = @_;

    my $self = $class->SUPER::new($args);
    if ( !$self->cache_exists() ) { die "Error: no data in cache??\n"; }
    my $date = strftime "%Y-%m-%d", localtime;
	my $str_count =  scalar keys %{$$self{selected_strains}};
	my $loc = '['.$$self{chrm}.':'.$$self{from}.'-'.$$self{to}.']';
	my $file = 'Indels'.$str_count."_mouse_strains_".$loc."_".$date.".csv";
    $$self{writer}->fname($file);
    return $self;
}

sub nonzero_column_data
{
    my ($self,$row) = @_;

    my ($pos,$chr,$ref,$gene_name);

    my $ncols = @$row;
    for (my $i=0; $i<$ncols; $i++)
    {
        if ( !$$row[$i] || !$$row[$i]->{'pos'} ) { next; }

        $pos  = $$row[$i]->{'pos'};
        $chr  = $$row[$i]->{'chr'};
        $ref = $$row[$i]->{'ref_base'};
        if ( exists($$row[$i]->{'consequence'}) )
        {
            $gene_name = $$row[$i]->{'gene_name'};
            if ( $gene_name && $pos ) { last; }
        }
    }

    return ($pos,$chr,$ref,$gene_name);
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
    $html->out("Gene,Chromosome,Position,Reference");
    for my $str (sort {$$strains{$a}<=>$$strains{$b}} keys %$strains)
    {
        $html->out(",$str,Consequence");
    }
    $html->out("\n");
    return;
}

sub print_row
{
    my ($self,$row) = @_;

    my $html = $$self{'writer'};
    my $session = $$self{'session'};

    my ($pos,$chr,$ref,$gene_name) = $self->nonzero_column_data($row);

    $html->out("$gene_name,$chr,$pos,$ref");

    my $ncols = scalar keys %{$$self{'selected_strains'}};
    for (my $i=0; $i<$ncols; $i++)
    {
        my $conseqs = {};
        for my $type (@{$$row[$i]->{'consequence'}})
        {
            if ( !$type || $type eq 'SPLICE_SITE' ) { next }  # ignore these - according to Dave these are rubbish
            $$conseqs{$type} = 1;
        }
        
        if ($$row[$i]->{'sequence'}) {
        	my @alts = split('/', $$row[$i]->{'alt_index'});
        	my $seq_out;
        	foreach (@alts) {
        		$seq_out = $seq_out ? $seq_out.'/' : '';
        		$seq_out = $seq_out.$$row[$i]->{'alt_sequences'}{$_};
        	}
        	$html->out("," . $seq_out . "," . join(';',sort keys %$conseqs));
        } 
        else {
        	$html->out(",-,-");
        }	
    }
    $html->out("\n");
}


1;

