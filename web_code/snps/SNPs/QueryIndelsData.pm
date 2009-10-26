#
# Author:    Petr Danecek (pd3@sanger.ac.uk)    Team 145
#
#--------------- QueryIndelsData ---------------------------------
#
# Takes the cached data and creates a tab-delimited output file. Same as QuerySNPsData, but for Indels
#

package SNPs::QueryIndelsData;

use strict;
use warnings;

use base qw(SNPs::QuerySNPsData);

sub new
{
    my ($class, $args) = @_;

    my $self = $class->SUPER::new($args);

    return $self;
}

sub nonzero_column_data
{
    my ($self,$row) = @_;

    my ($pos,$chr,$sequence,$type,$gene_name,$gene_id);

    $sequence = '*';
    my $ncols = @$row;
    for (my $i=0; $i<$ncols; $i++)
    {
        if ( !$$row[$i] || !$$row[$i]->{'pos'} ) { next; }

        $pos  = $$row[$i]->{'pos'};
        $chr  = $$row[$i]->{'chr'};
        $type = $$row[$i]->{'type'};
        if ( $type && $type eq 'D' ) 
        { 
            $sequence=$$row[$i]->{'sequence1'}; 
        }
        if ( exists($$row[$i]->{'_conseqs'}) )
        {
            my $conseqs = $$row[$i]->{'_conseqs'};
            $gene_id   = $$conseqs[0]->{'gene_ensid'} || '';
            $gene_name = $$conseqs[0]->{'gene_name'} || $gene_id;
            if ( $gene_id && $pos ) { last; }
        }
    }

    return ($pos,$chr,$sequence,$gene_name,$gene_id);
}



sub print_row
{
    my ($self,$row) = @_;

    my $html = $$self{'writer'};
    my $session = $$self{'session'};

    my ($pos,$chr,$ref,$gene_name,$gene_id) = $self->nonzero_column_data($row);

    $html->out("$gene_name\t$chr\t$pos\t$ref");

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

        if ( $$row[$i]->{'sequence1'} )
        {
            my $sequence;

            if ( $$row[$i]->{type} eq 'I' ) { $sequence=$$row[$i]->{sequence1}; }
            else { $sequence='*'; }
            if ( $$row[$i]->{sequence2} )
            {
                if ( $$row[$i]->{sequence2} eq '*' )
                {
                    $sequence .= '/' . $ref;
                }
                else
                {
                    $sequence .= '/' . $$row[$i]->{sequence2};
                }
            }

            $html->out("\t" . $sequence);
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

