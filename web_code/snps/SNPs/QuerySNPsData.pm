#
# Author:    Petr Danecek (pd3@sanger.ac.uk)    Team 145
#
#--------------- QuerySNPsData ---------------------------------------
#
# Takes the cached data and creates a tab-delimited output file.
#

package SNPs::QuerySNPsData;

use strict;
use warnings;
use base qw(SNPs::QuerySNPs);

sub new
{
    my ($class, $args) = @_;

    my $self = $class->SUPER::new($args);

    if ( !$self->cache_exists() ) { die "Error: no data in cache??\n"; }
    $$self{writer}->fname($$self{writer}{cgi}->param('cache') . '.tab');

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
    $html->out("Gene\tChromosome\tPosition\tReference");
    for my $str (sort {$$strains{$a}<=>$$strains{$b}} keys %$strains)
    {
        $html->out("\t$str\tConsequence");
    }
    $html->out("\n");
    return;
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

    my ($pos,$chr,$base,$gene_name,$gene_id) = $self->nonzero_column_data($row);

    $html->out("$gene_name\t$chr\t$pos\t$base");

    my $ncols = scalar keys %{$$self{'selected_strains'}};
    for (my $i=0; $i<$ncols; $i++)
    {
        my $conseqs = {};
        for my $cons (@{$$row[$i]->{'_conseqs'}})
        {
            my $type = $$cons{'consequence'};
            if ( !$type || $type eq 'SPLICE_SITE' ) { next }  # ignore these - according to Dave these are rubbish
            $$conseqs{$type} = 1;
        }

        if ( $$row[$i]->{'depth'} )
        {
            $html->out("\t" . $$row[$i]->{'snp_base1'});
            if ( $$row[$i]->{'snp_base2'} ) { $html->out('/'.$$row[$i]->{'snp_base2'}); }
            $html->out("\t" . join(',',sort keys %$conseqs));
        }
        else
        {
            $html->out("\t-\t");
        }
    }
    $html->out("\n");
}

1;

