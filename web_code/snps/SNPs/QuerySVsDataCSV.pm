#
# Author:    John Maslen (jm23@sanger.ac.uk)    Team 145
#
#--------------- QuerySVsData ---------------------------------
#
# Takes the cached data and creates a tab-delimited output file. Same as QuerySNPsData, but for SVs
#

package SNPs::QuerySVsDataCSV;

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
	my $file = 'SVs'.$str_count."_mouse_strains_".$loc."_".$date.".csv";
    $$self{writer}->fname($file);
    return $self;
}

sub print_header
{
    my ($self) = @_;
    my $strains = $$self{selected_strains};
    my $html = $$self{writer};
    if ( $$self{display_dload_svparams} )
    {
        $html->out($$self{display_dload_svparams});
    }
    $html->out("Chromosome,Grouped position");
    for my $str (sort {$$strains{$a}<=>$$strains{$b}} keys %$strains)
    {
        $html->out(",$str");
    }
    $html->out("\n");
    return;
}

sub print_row
{
    my ($self,$row) = @_;

    my $html = $$self{'writer'};
    my $session = $$self{'session'};

    my ($pos,$endpos,$chr,$sv_type, $display_token, $display_type);
    my $ncols = @$row;
    for (my $i=0; $i<$ncols; $i++)
    {
        if ( !$$row[$i] || !$$row[$i]->{pos} ) { next; }

        $pos = $$row[$i]->{'pos'};
        $endpos = $$row[$i]->{'endpos'};
        $chr = $$row[$i]->{'chr'};
    }
    $html->out(qq[$chr,$pos-$endpos]);

	my $sv_out;
    $ncols = scalar keys %{$$self{'selected_strains'}};
    for (my $i=0; $i<$ncols; $i++)
    {
        my $sv_out = '-';
        if ( $$row[$i]->{'sv_type'} ) {
        	$sv_out = $$row[$i]->{'display'};
        	$sv_out =~ s/\W+/_/g;
        	$sv_out = $sv_out.';'.$$row[$i]->{'strain_loc'};
        }	
        $html->out(qq[,$sv_out]);
    }
    $html->out("\n");
}

1;

