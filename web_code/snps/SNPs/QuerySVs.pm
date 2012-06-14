#
# Author:    John Maslen (jm23@sanger.ac.uk)    Team 145
#
#--------------- QuerySVs---------------------------------------
#
# Builds and executes the TABIX command to retrieve SV data via FTP
#

package SNPs::QuerySVs;

use strict;
use warnings;
use base qw(SNPs::QuerySNPs);

sub new
{
    my ($class, $args) = @_;

    $$args{title}        = 'Structural Variations' unless exists($$args{'title'});
    $$args{print_legend} = 1 unless exists($$args{'print_legend'});

    my $self = $class->SUPER::new($args);

    $$self{'action'}     = 'structvar';

    return $self;
}


sub run_tabix_command
{
	my ($self) = @_;
	my $tabix_cmd = "cd $$self{'index_directory'} ; $$self{'tabix_location'} -h $$self{'svs_location'}  $$self{'pos_str'}";
	my $mouse_strains = $$self{'vcf_strain_order'};
	my $strains = $$self{'selected_strains'};	
	my @sort_strains = (sort {$$strains{$a}<=>$$strains{$b}} keys %$strains);	
	my $input = `$tabix_cmd`;
	my @lines = split('\n', $input); 
	my %output = ();
	my %sv_types_selected = %{$$self{'selected_svtypes'}};	
	my %strains_selected = %{$$self{'selected_strains'}};
	my $linecount = 1;	
	
	foreach ( @lines ) {
		my @line = split('\t', $_);
		my @lineoutput = ();
    	my ($display_label, $display_type, $display);
    	my $chr = $line[0];
    	my $from = $line[1];
    	my $to = $line[2];
    	my @strain_data = @line[3 .. 19];
    	die "Mouse strains in vcf file do not match expected strains" unless scalar(@strain_data) == scalar( @{$mouse_strains}); 
    	for my $i ( 0 .. $#strain_data ) {
    		my %linehash;
			my $strain = $$mouse_strains[$i];
			my @strain_sv = split(';', $strain_data[$i]);
				for my $key (keys %{$$self{'sv_type_map'}})
    			{
        			if ($strain_sv[1] =~ /$key/) {
        				my $svtype = $$self{'sv_type_map'}{$key};
        			    $display_label = $$svtype{'label'};
        				$display_type=$$svtype{'type'};
        				$display=$$svtype{'display'};
        			}
    			}
    				$linehash{'display_label'} = $display_label;
        			if($sv_types_selected{$display_type}) {
        				$linehash{'display_type'} = $display_type;
        			}	
        			$linehash{'display'} = $display;    
					$linehash{'strain_name'} = $strain;
	    			$linehash{'chr'} = $chr;
        			$linehash{'pos'} = $from+1;
        			$linehash{'endpos'} = $to;
        			if (scalar @strain_sv == 4) {
        				my @chr_pos = split(':', $strain_sv[0]);
        				my @pos = split('-', $chr_pos[1]);
        				$linehash{'strain_loc'} = $chr_pos[0].':'.($pos[0]+1).'-'.$pos[1];
        				$linehash{'sv_type'} = $strain_sv[1];
        				$linehash{'brkpt'} = $strain_sv[2];
        				$linehash{'pem_pattern'} = $strain_sv[3];
        			}
        			if ( $$self{mouseinfo}{$strain}{bam} ) {
    					$linehash{'bam'} = $$self{mouseinfo}{$strain}{bam};
        			}
        			push @lineoutput, {%linehash};
  		}
    	if (@lineoutput) {
    		my @strain_output = ();
    		my $sv_found = 0;
    	    for my $str ( @sort_strains ) {
    	    	for my $str_res ( @lineoutput ) {
    	    		if ( $$str_res{'strain_name'} eq $str ) {
    	    			push @strain_output, $str_res;
    	    			if ( !$sv_found && exists $$str_res{'sv_type'} && exists $$str_res{'display_type'} ) {
    	    				$sv_found = 1;
    	    			}
    	    		}
    	    	}
    	    }
    	    if ( $sv_found ) {
    			$output{$linecount} = [@strain_output];
    			$linecount++;
    		}	
    	}  			
	}
	
	return \%output;
}	


sub print_no_match
{
    my ($self) = @_;
    die("SVs not found ;;; Sorry, unable to find SVs for $$self{loc}.\n");
}

sub print_legend
{
    my ($self) = @_;

    if ( !$$self{print_legend} ) { return; }

    # TODO: position taken from the first result on the page??
    my $html = $$self{'writer'};
    $html->out(qq[<div id="legend">
            <b>Try:</b>
                <div style="margin-left:1em;"> <a href="$$self{'myself'}">New search</a> </div>
                <div style="margin-left:1em;"> <a href="http://www.ensembl.org/Mus_musculus/Location/View?db=core;r=$$self{chrm}:$$self{from}-$$self{to};" target="_ensembl_snp">View in Ensembl</a> </div>
                <div style="margin-left:1em;"> Click on SVs for details</div>
            <div style="padding-top:1em;"><b>Download:</b></div>
        ]);
    my $session = $$self{session}->id();
    $html->out(qq[<div style="margin-left:1em;">Formats: <a href="$$self{'myself'}?cache=$session&action=$$self{action}_dload">tab</a> or <a href="$$self{'myself'}?cache=$session&action=$$self{action}_dload_csv">csv</a>.\n </div>]);


	$html->out(qq[<div style="padding-top:1em;"><b>Reference:</b></div>
				  <div style="margin-left:.2em;">Yalcin, B. <i>et al.</i> Sequence-based characterization of structural variation in the mouse genome.  <a href="http://dx.doi.org/10.1038/nature10432" target="_blank"><i>Nature</i></a> 477, 326-329 (2011).</div>]);
    
    $html->out(qq[
        <div style="padding-top:1em;">
        <b>Legend:</b>
        <table>
        ]);
    for my $key (sort keys %{$$self{'selected_svtypes'}})
    {
        my $svtype = $$self{'sv_types'}{$key};
        $html->out(qq[
            <tr><td><img class="$$svtype{'style'}" width="10" height="10" alt="" src="/icons/blank.gif"></td>
                <td>$$svtype{'label'}</td>\n]);
    }
    $html->out(qq[</table></div>]);
    
	$html->out(qq[<div style="padding-top:1em;"><b>Complex events key:</b></div><table rules="rows">]);	
    for my $key (sort keys %{$$self{'sv_combined_events_key'}})
    {
        my $svkey = $$self{'sv_combined_events_key'}{$key};
        $html->out(qq[<tr><td><b><i>$key<ii></b></td><td>$$svkey{'label'}</td></tr>]);
    }					
    $html->out(qq[</table>]);    
    
    $html->out(qq[<div style="padding-top:1em;"><b>Breakpoint definitions:</b></div>
				  <div style="margin-left:.4em;"><b><i>REF</i></b> = 'refined by local assembly analysis'<br />
					<b><i>RAW</i></b> = 'taken from the original SV calls'</div>]);				  
}

sub print_header
{
    my ($self) = @_;

    # Run do-labels.pl to generate all column labels using the cimg.pl script.
    # No reason to run this on fly at each request: The script do-labels will
    # create them by calling:
    #       http://wwwdev.sanger.ac.uk/cgi-bin/modelorgs/mousegeneomes/cimg.pl?t=Gene

    my $url_imgs   = $$self{'imgs'};
    my $strains    = $$self{'selected_strains'};
    my $html       = $$self{'writer'};
    my $session    = $$self{'session'};

    # $session->append("Gene\tChromosome\tPosition\tReference");
    $html->out(qq[
        <table class="results"><thead><tr>
        <th class="results"><img src="$url_imgs/chromosome.png" alt="Chromosome" /></th>
        <th class="results"><img src="$url_imgs/position.png" alt="Position" /></th>
        ]);
    for my $str (sort {$$strains{$a}<=>$$strains{$b}} keys %$strains)
    {
        my $img_name = lc($str);
        $img_name =~ s{/}{_}g;
        $html->out(qq[<th class="results"><img src="$url_imgs/$img_name.png" alt="$str" /></th>]);
        # $session->append("\t$str");
    }
    $html->out("</tr></thead><tbody>\n");
    # $session->append("\n");
}

sub print_row
{
    my ($self,$row) = @_;

    my $html = $$self{'writer'};
    my $session = $$self{'session'};
    my ($pos,$endpos,$chr);
    my $ncols = @$row;
    for (my $i=0; $i<$ncols; $i++)
    {
        if ( !$$row[$i] || !$$row[$i]->{pos} ) { next; }

        $pos = $$row[$i]->{'pos'};
        $endpos = $$row[$i]->{'endpos'};
        $chr = $$row[$i]->{'chr'};
    }

    my $print_gene = 0;
    my $style      = 0;
    my @styles = ('class="gene1"','class="gene2"','');
    $style  = $styles[$style];
    my $ref = '*';

    
    $html->out(qq[<tr><td $style>$chr</td><td $style>$pos-$endpos</td>]);
    $ncols = scalar keys %{$$self{'selected_strains'}};
    for (my $i=0; $i<$ncols; $i++)
    {
        my $reported     = {};
        my $details      = '';
        my $style = '';
        my $onclick='';
        my $div='';
        if ( $$row[$i]->{'sv_type'} )
        {
            $style   = $$self{'sv_types'}{$$row[$i]->{'display_type'}}{'style'};
            $details = 
                "<tr><td>Strain:</td><td>" .$$row[$i]->{strain_name}. "</td></tr>" .
                "<tr><td>Location:</td><td>".$$row[$i]->{strain_loc}."</td></tr>" . $details;
            $details .= '<tr><td>SV type:</td><td> ' . $$row[$i]->{'display'}. '</td></tr>';
            $details .= '<tr><td>Breakpoint:</td><td> ' . $$row[$i]->{'brkpt'} . '</td></tr>';
            $details .= qq[<tr><td colspan="2"><a href="http://www.ensembl.org/Mus_musculus/Location/View?db=core;r=$chr:$pos-$endpos;" target="_ensembl">View in Ensembl</a></td></tr>];
            if ( $$row[$i]->{bam} )
            {
                my $from = $pos - 1000;
                my $to = $endpos + 1000;
                my $winsize = $to - $from;
                my $insert_size = ($endpos - $pos) + 500;
                if ( $from < 0 ) { $from = 0; }
                my $lookseq =
                    "$$self{lookseq}?show=$chr:$from-$to,indel&amp;lane=" . $$row[$i]->{bam} . "&amp;win=$winsize&amp;width=$$self{lseq_width}&amp;$$self{lseq_display_sv}&amp;maxdist=$insert_size";
                $details .= qq[<tr><td colspan="2"><a href="$lookseq" target="_lookseq">View in LookSeq</a></td></tr>];
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

        my $var_type = '';
        if ( !$$row[$i]->{'sv_type'} ) { $var_type='-'; }
        else
        {
            $var_type = $$row[$i]->{'display_label'};
        }
        $html->out(qq[<td $style $onclick>$var_type$div</td>]);
    }
    $html->out("</tr>\n");
}

1;

