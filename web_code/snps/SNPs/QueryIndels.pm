#
# Author:    Petr Danecek (pd3@sanger.ac.uk)    Team 145
#
#--------------- QueryIndels ---------------------------------------
#
# Builds and executes the SQL query for indels.
#

package SNPs::QueryIndels;

use strict;
use warnings;
use base qw(SNPs::QuerySNPs);

sub new
{
    my ($class, $args) = @_;

    $$args{title}        = 'Indels' unless exists($$args{'title'});
    $$args{print_legend} = 1 unless exists($$args{'print_legend'});

    my $self = $class->SUPER::new($args);

    $$self{'action'}     = 'indels';

    return $self;
}

sub sql_query
{
    my ($self,$filters) = @_;

    my $query = qq[
        SELECT 
            r.value AS strain_name, 
            s1.value AS sequence1, s2.value AS sequence2,
            chr.value AS chr,
            s.pos, s.type, s.length, s.depth, s.snp_qual, s.map_qual, s.cons_qual,
            g.value AS gene_name,
            c.*
        FROM 
            strains r, chroms chr, indels s
            LEFT JOIN indel_conseqs c ON s.indel_id=c.indel_id
            LEFT JOIN genes g ON c.gene=g.id
            LEFT JOIN strings s1 ON s.sequence1=s1.id
            LEFT JOIN strings s2 ON s.sequence2=s2.id
        WHERE
            ] .join("\n AND ",@$filters). q[ ORDER BY s.chr ASC, s.pos ASC];

    return $query;
}


sub print_no_match
{
    my ($self) = @_;
    die("Sorry, no matching indels.\n");
}


sub print_row
{
    my ($self,$row) = @_;

    my $html = $$self{'writer'};
    my $session = $$self{'session'};

    # Find any non-empty column to get the information common to all columns: Gene, chromosome, position,
    #   the reference base.
    #
    my ($pos,$chr,$seq,$gene_id,$insertion);
    my $gene_name = '';
    my $ncols = @$row;
    for (my $i=0; $i<$ncols; $i++)
    {
        if ( !$$row[$i] || !$$row[$i]->{pos} ) { next; }

        $pos = $$row[$i]->{'pos'};
        $chr = $$row[$i]->{'chr'};
        $seq = $$row[$i]->{'sequence1'};
        $insertion = $$row[$i]->{'type'} eq 'I' ? 1 : 0;
        if ( exists($$row[$i]->{'_conseqs'}) )
        {
            my $conseqs = $$row[$i]->{'_conseqs'};
            $gene_id   = $$conseqs[0]->{'gene_ensid'} || '';
            $gene_name = $$conseqs[0]->{'gene_name'} || $gene_id;
            if ( $gene_id ) { last; }
        }
    }

    # Alternate background for each gene and print gene name only if different from the previous row
    #
    my $print_gene = 0;
    my $style      = 0;
    if ( !$gene_id ) { $style=-1; }
    elsif ( !$$self{'last_gene_id'} ) { $style=0; $print_gene=1; }
    elsif ( $$self{'last_gene_id'} eq $gene_id ) { $style=$$self{'last_gene_style'}; }
    else { $style = $$self{'last_gene_style'} ? 0 : 1; $print_gene=1; }

    $$self{'last_gene_id'}    = $gene_id;
    $$self{'last_gene_style'} = $style;
    
    my @styles = ('class="gene1"','class="gene2"','');
    $style  = $styles[$style];
    my $ref = $insertion ? '*' : $seq;

    $html->out("<tr><td $style>");
    if ( $print_gene ) { $html->out(sprintf qq[<a href="http://www.ensembl.org/Mus_musculus/Gene/Summary?g=ENSMUSG%.11d">$gene_name</a>], $gene_id); }
    $html->out(qq[</td><td $style>$chr</td><td $style>$pos</td><td class="resultsbar">$ref</td>]);
    # $session->append("$gene_name\t$chr\t$pos\t$ref");

    $ncols = scalar keys %{$$self{'selected_strains'}};
    for (my $i=0; $i<$ncols; $i++)
    {
        # Selected from all the possible consequence types one for coloured display.
        #   and the rest will be reported in a popup div as a list.
        #
        my $reported     = {};
        my $conseq_type  = '';
        my $details      = '';
        for my $cons (@{$$row[$i]->{'_conseqs'}})
        {
            my $type = $$cons{'consequence'};

            if ( !$type ) { next }
            if ( $type eq 'SPLICE_SITE' ) { next }  # ignore these - according to Dave these are rubbish

            if ( exists($$reported{$type}) ) { next } # report each type only once
                $$reported{$type} = 1;

            # UTR consequence types have lower priority, but we still want to keep it if there is nothing else
            #
            if ( $conseq_type && $type =~ /UTR$/ ) { next }
            $conseq_type = $type;
        }

        # The SNPs will be coloured according to its consequence type.
        #
        my $style = '';
        if ( $conseq_type )
        {
            my @types = ();
            for my $key (sort keys %$reported)
            {
                if ( !exists($$self{'conseqs'}{$key}) ) { die("FIXME: The consequence type \"$key\" not covered.\n"); }
                push @types, qq[<span class="$$self{'conseqs'}{$key}{'style'}">$$self{'conseqs'}{$key}{'label'}</span>];
            }
            $details = '<tr><td>Type:</td><td>' . join('<br />',@types) . '</td></tr>';
            $style   = $$self{'conseqs'}{$conseq_type}{'style'};
        }

        my $onclick='';
        my $div='';
        if ( $$row[$i]->{'depth'} )
        {
            # If we are here, there is some info, not an empty row.
            #
            $details = 
                "<tr><td>Strain:</td><td>" .$$row[$i]->{strain_name}. "</td></tr>" .
                "<tr><td>Location:</td><td>$chr:$pos</td></tr>" . $details;
            $details .= '<tr><td>Depth:</td><td> ' . $$row[$i]->{'depth'} . '</td></tr>';
            $details .= '<tr><td>Quality:</td><td> ' . $$row[$i]->{'snp_qual'} .'/'. $$row[$i]->{'map_qual'} .'/'. $$row[$i]->{'cons_qual'} . '</td></tr>';
            $details .= qq[<tr><td colspan="2"><a href="http://www.ensembl.org/Mus_musculus/Location/View?db=core;r=$chr:$pos-$pos;">View in Ensembl</a></td></tr>];
            if ( $$row[$i]->{bam} )
            {
                my $from = $pos - int($$self{lseq_win}*0.5);
                if ( $from < 0 ) { $from = 0; }

                my $lookseq =
                    "$$self{lookseq}?show=$chr:$from-$from,paired_pileup&amp;win=$$self{lseq_win}&amp;width=$$self{lseq_width}" .
                    "&amp;$$self{lseq_display}&amp;lane=" . $$row[$i]->{bam};
                $details .= qq[<tr><td colspan="2"><a href="$lookseq">View in LookSeq</a></td></tr>];
            }


            my $id = 'ud' . (++$$html{'unique_div'});
            $onclick = qq[onclick="toggle_one_element('#$id')"];
            $div = qq[
                <div class="hidden" id="$id">
                <div style="text-align:right;font-size:smaller;margin-top:-0.5em; margin-right:-0.5em;">
                        <span class="button">[x]</span></div>
                <table class="details">
                    $details
                </table>
                </div>];

            $style = qq[class="$style button"];
        }

        my $sequence = '';
        if ( !$$row[$i]->{'sequence1'} ) { $sequence='-'; }
        else
        {
            if ( $insertion ) { $sequence=$seq; }
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
        }
        if ( keys %$reported > 1 ) { $sequence .= '<sup>+</sup>'; }
        $html->out(qq[<td $style $onclick>$sequence$div</td>]);
    }
    $html->out("</tr>\n");
}

1;

