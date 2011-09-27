#
# Author:    Petr Danecek (pd3@sanger.ac.uk)    Team 145
# Modified:		John Maslen  (jm23@sanger.ac.uk)   Team 145
#
#--------------- QueryIndels ---------------------------------------
#
# Builds and executes the TABIX command to retrieve indel data from a vcf file via FTP
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

sub run_tabix_command
{
	my ($self) = @_;
	my $tabix_cmd = 'cd '.$$self{'index_directory'}.' ; '.$$self{'tabix_location'}.' -h '.$$self{'indels_location'}.' '.$$self{'pos_str'};
	my @mouse_strains2 = ();
	my $strains = $$self{'selected_strains'};	
	my @sort_strains = (sort {$$strains{$a}<=>$$strains{$b}} keys %$strains);	
	my $input = `$tabix_cmd`;
	my @lines = split('\n', $input); 
	my %output = ();
	my %conseq_selected = %{$$self{'selected_conseq'}};
	my $linecount = 1;	

	for my $line_in ( @lines ) {
		if ( $line_in =~ /^##\w+/ ) {
			next;
		}
		elsif ( $line_in =~ /^#CHROM/ ) {
			my @strain_line = split('\t', $line_in); 
			for my $vcf_str (@strain_line) {
				if ($$self{vcf_strain_map}{$vcf_str}{name}) {
					push @mouse_strains2, $$self{vcf_strain_map}{$vcf_str}{name};
				}
			}
		}
		else{
		my @line = split('\t', $line_in);
		my @lineoutput = ();	
		my %alt_sequences = ();	
		my $info_csq = return_info_field($line[7], 'CSQ');
		my (%types, %gene_conseq);
    	if ( $info_csq ) {
        	for my $cons (split('\+', $info_csq)) {
        		my @gene_type = split(/[:,]/, $cons);
        		if (!$types{$gene_type[2]} && $conseq_selected{$gene_type[2]}) {
        			$types{$gene_type[2]} = 1;
        			push @{$gene_conseq{$gene_type[1]}}, $gene_type[2];
        		}
        	}
    	}
    	elsif ($conseq_selected{'INTERGENIC'}) {
    		push @{$gene_conseq{'-'}}, 'INTERGENIC';
    	} 
    		
	    my @strain_data = @line[9 .. 25];
    	my %alleles = ();
    	my $allele_count = 1;
    	my $ref_seq = $line[3];
    	for my $alt_seq ( split(',', $line[4] )) {
    		$alt_sequences{$allele_count} = $alt_seq;	
            my $seq_out = '';
    		if (length($alt_seq) < length($ref_seq)) {
    			$seq_out = '-'.substr($ref_seq, length($alt_seq));
    		}
    		elsif (length($alt_seq) > length($ref_seq)) {
    			$seq_out = '+'.substr($alt_seq, length($ref_seq));
    		}
    		$alleles{$allele_count} = $seq_out;
    		$allele_count++;
    	}
    	die "Mouse strains in vcf file do not match expected strains" unless scalar(@strain_data) == scalar( @mouse_strains2); 
    	for my $i ( 0 .. $#strain_data ) {
    		my %linehash;
			my @data = split (':', $strain_data[$i]);
			my $strain = $mouse_strains2[$i];
    		for my $gene (keys %gene_conseq) {
    			my $cons_qual = $data[1];
    			my %allele_hash;
    			my @alt_indices;
    			my @alt_seq;
    			foreach (split('/', $data[0])) {
    				if (!$allele_hash{$alt_sequences{$_}}) {
    					$allele_hash{$alt_sequences{$_}} = 1;
    					push @alt_indices,$_;
    					push @alt_seq, $alleles{$_};
    				}
    			}
    			$linehash{'strain_name'} = $strain;
    			$linehash{'consequence'} = $gene_conseq{$gene};
    			$linehash{'gene_name'} = $gene unless $gene eq '-';
    			$linehash{'chr'} = $line[0];
    			$linehash{'pos'} = $line[1];
    			$linehash{'ref_base'} = $line[3];
    			$linehash{'alt_index'} = join('/',@alt_indices);
    			$linehash{'alt_sequences'} = \%alt_sequences; 
    			$linehash{'sequence'} = join('/',@alt_seq);
    			$linehash{'type'} = (length($linehash{'sequence'})>length($linehash{'ref_base'}) ? 'I' : 'D');
    			$linehash{'seq_size'} = $linehash{'sequence'} ? length($linehash{'sequence'}) : -1;
    			if ($data[1]) {
    				$linehash{'cons_qual'} = $cons_qual;
    			}	
    			if ( $$self{mouseinfo}{$strain}{bam} ) {
    				$linehash{'bam'} = $$self{mouseinfo}{$strain}{bam};
        		}    			
    		}
    		if (%linehash) {
    			push @lineoutput, {%linehash};
    		}	
    	}
    	if (@lineoutput) {
    		my @strain_output = ();
    		my $alt_idx_found = 0;
    	    for my $str ( @sort_strains ) {
    	    	for my $str_res ( @lineoutput ) {
    	    		if ( $$str_res{'strain_name'} eq $str ) {
    	    			push @strain_output, $str_res;
    	    			if ( !$alt_idx_found && $$str_res{'seq_size'} > 0) {
    	    				$alt_idx_found = 1;
    	    			}    	    		
    	    		}
    	    	}
    	    }
    	    if ( $alt_idx_found ) {
    			$output{$linecount} = [@strain_output];
    			$linecount++;
    		}
    	}    	
    }
    }
	return \%output;
}	

sub print_no_match
{
    my ($self) = @_;
    die("Sorry, no matching indels.\n");
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
        <th class="results"><img src="$url_imgs/gene.png" alt="Gene" /></th>
        <th class="results"><img src="$url_imgs/chromosome.png" alt="Chromosome" /></th>
        <th class="results"><img src="$url_imgs/position.png" alt="Position" /></th>
        <th class="results resultsbar"><img src="$url_imgs/reference.png" alt="Reference" /></th>
        <th class="results resultsbar"><img src="$url_imgs/alternate.png" alt="Alternate" /></th>
        ]);
    for my $str (sort {$$strains{$a}<=>$$strains{$b}} keys %$strains)
    {
        my $img_name = lc($str);
        $img_name =~ s{/}{_}g;
        $html->out(qq[<th class="results"><img src="$url_imgs/$img_name.png" alt="$str" /></th>]);
    }
    $html->out("</tr></thead><tbody>\n");
}

sub print_legend
{
    my ($self) = @_;

    if ( !$$self{print_legend} ) { return; }

    # TODO: position taken from the first result on the page??
    my $html = $$self{'writer'};
    my $lookseq_win = $$self{to}-$$self{from};
    if ($lookseq_win > 100000) {
    	$lookseq_win = 100000;
    }
    my $lookseq_params = 
        "show=$$self{chrm}:$$self{from}-$$self{to},paired_pileup" .
        "&amp;win=$lookseq_win" .
        "&amp;width=$$self{lseq_width}" .
        "&amp;$$self{lseq_display}" .
        "&amp;lane=DBA.bam";
    $html->out(qq[<div id="legend">
            <b>Try:</b>
                <div style="margin-left:1em;"> <a href="$$self{'myself'}">New search</a> </div>
                <div style="margin-left:1em;"> <a href="http://www.ensembl.org/Mus_musculus/Location/View?db=core;r=$$self{chrm}:$$self{from}-$$self{to};" target="_ensembl_snp">View in Ensembl</a> </div>
                <div style="margin-left:1em;"> <a href="$$self{'lookseq'}?$lookseq_params" target="_lookseq_snp">View in LookSeq</a> </div>
                <div style="margin-left:1em;"> Click on Indels for details</div>
            <div style="padding-top:1em;"><b>Download:</b></div>
        ]);

    my $session = $$self{session}->id();
	$html->out(qq[<div style="margin-left:1em;">Formats: <a href="$$self{'myself'}?cache=$session&action=$$self{action}_dload">tab</a> or <a href="$$self{'myself'}?cache=$session&action=$$self{action}_dload_csv">csv</a>.\n </div>]);

    $html->out(qq[
        <div style="padding-top:1em;">
        <b>Legend:</b><br/>
         &nbsp;&nbsp;<b>Indel Consequences:</b>
        <table>]);
    for my $key (sort keys %{$$self{'selected_conseq'}})
    {
        my $conseq = $$self{'conseqs'}{$key};
        $html->out(qq[
            <tr><td><img class="$$conseq{'style'}" width="10" height="10" alt="" src="/icons/blank.gif"></td>
                <td>$$conseq{'label'}</td>\n]);
    }
    $html->out(qq[<tr><td><sup>+</sup></td><td>Multiple consequences</td></tr>]);
    $html->out(qq[</table></div></div>]);
}

sub print_row
{
    my ($self,$row) = @_;

    my $html = $$self{'writer'};
    my $session = $$self{'session'};

    # Find any non-empty column to get the information common to all columns: Gene, chromosome, position,
    #   the reference base.
    #
    my ($pos,$chr,$seq,$ref,$alt_seq,$insertion);
    my $gene_name = '';
    my $ncols = @$row;
    for (my $i=0; $i<$ncols; $i++)
    {
        if ( !$$row[$i] || !$$row[$i]->{pos} || !$$row[$i]{'alt_index'}) { next; }

        $pos = $$row[$i]->{'pos'};
        $chr = $$row[$i]->{'chr'};
        $seq = $$row[$i]->{'sequence'};
        $ref = $$row[$i]->{'ref_base'};
        $alt_seq = $$row[$i]->{'alt_sequences'};
        $insertion = $$row[$i]->{'type'} eq 'I' ? 1 : 0;
        if ( exists($$row[$i]->{'consequence'}) )
        {
            $gene_name = $$row[$i]->{'gene_name'};
            if ( $gene_name && $pos && $alt_seq ) { last; }
        }
    }

    # Alternate background for each gene and print gene name only if different from the previous row
    #
    my $print_gene = 0;
    my $style      = 0;
    if ( !$gene_name ) { $style=-1; }
    elsif ( !$$self{'last_gene_name'} ) { $style=0; $print_gene=1; }
    elsif ( $$self{'last_gene_name'} eq $gene_name ) { $style=$$self{'last_gene_style'}; }
    else { $style = $$self{'last_gene_style'} ? 0 : 1; $print_gene=1; }

    $$self{'last_gene_name'}    = $gene_name;
    $$self{'last_gene_style'} = $style;
    
    my @styles = ('class="gene1"','class="gene2"','');
    $style  = $styles[$style];
    $html->out("<tr><td $style>");
    
    my $alt_seq_size = scalar keys %{$alt_seq};
    my $current_alt = 1;
    if ( $print_gene ) { $html->out(sprintf qq[<a href="http://www.ensembl.org/Mus_musculus/Gene/Summary?g=$gene_name" target="_ensembl_snp">$gene_name</a>]); }    
    
   	$html->out(qq[</td><td $style>$chr</td><td $style>$pos</td><td style="text-align: left;" ]);
   	
   	if ($alt_seq_size == $current_alt) {
   		$html->out(qq[class="resultsbar resultsdemarc">$ref</td><td style="white-space:nowrap; text-align: left;" class="resultsbar resultsdemarc alt_alleles">]);
   	}
   	else {
   		$html->out(qq[class="resultsbar alt_alleles">$ref</td><td style="white-space:nowrap; text-align: left;" class="resultsbar">]);
   	}
    if ($$alt_seq{$current_alt}) {
    	$html->out(qq[$current_alt = $$alt_seq{$current_alt}</td>]);
    }

    $ncols = scalar keys %{$$self{'selected_strains'}};
    for (my $i=0; $i<$ncols; $i++)
    {
        # Selected from all the possible consequence types one for coloured display.
        #   and the rest will be reported in a popup div as a list.
        #
        my $reported     = {};
        my $conseq_type  = '';
        my $details      = '';
        for my $type (@{$$row[$i]->{'consequence'}})
        {
            #my $type = $$cons{'consequence'};

            if ( !$type || $type eq 'SPLICE_SITE' ) { next }  # ignore these - according to Dave these are rubbish

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
        if ( $$row[$i]->{'sequence'} )
        {
            # If we are here, there is some info, not an empty row.
            #
            $details = 
                "<tr><td>Strain:</td><td>" .$$row[$i]->{strain_name}. "</td></tr>" .
                "<tr><td>Location:</td><td>$chr:$pos</td></tr>" . $details;

            $details .= '<tr><td>Alleles:</td><td> ' . $$row[$i]->{'sequence'}. '</td></tr>';
			$details .= '<tr><td>GT quality:</td><td> ' . $$row[$i]->{'cons_qual'} . '</td></tr>';            
            $details .= qq[<tr><td colspan="2"><a href="http://www.ensembl.org/Mus_musculus/Location/View?db=core;r=$chr:$pos-$pos;" target="_ensembl_snp">View in Ensembl</a></td></tr>];
            if ( $$row[$i]->{bam} )
            {
                my $from = $pos - int($$self{lseq_win}*0.5);
                if ( $from < 0 ) { $from = 0; }

                my $lookseq =
                    "$$self{lookseq}?show=$chr:$from-$from,paired_pileup&amp;win=$$self{lseq_win}&amp;width=$$self{lseq_width}" .
                    "&amp;$$self{lseq_display}&amp;lane=" . $$row[$i]->{bam};
                $details .= qq[<tr><td colspan="2"><a href="$lookseq" target="_lookseq_snp">View in LookSeq</a></td></tr>];
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

        my $sequence = $$row[$i]->{'sequence'} ? $$row[$i]->{'alt_index'} : '-';
        if ( keys %$reported > 1 && $$row[$i]->{'sequence'}) { $sequence .= '<sup>+</sup>'; }
        if ($alt_seq_size == $current_alt) {
        	$html->out(qq[<td $style $onclick class="resultsdemarc">$sequence$div</td>]);
        }
        else {
        	$html->out(qq[<td $style $onclick>$sequence$div</td>]);
        }	
    }
    $html->out("</tr>\n");
    
    if ($alt_seq_size > 1) {
    	for my $i ( 2 .. $alt_seq_size ) {
    		$current_alt++;
        	
        	$html->out(qq[<tr><td $style></td><td $style></td><td $style></td>]);
        	if ($alt_seq_size == $current_alt) {
        		$html->out(qq[<td class="resultsbar resultsdemarc"></td><td style="text-align: left;" class="resultsbar resultsdemarc">$i = $$alt_seq{$i}</td><td colspan="$ncols"  class="resultsdemarc"></td></tr>\n]);
        	}
        	else {
   				$html->out(qq[<td class="resultsbar"></td><td style="text-align: left;" class="resultsbar">$i = $$alt_seq{$i}</td><td colspan="$ncols"></td></tr>\n]);
   			}	
        }
    }
    
}

sub return_info_field
{
    my ($info,$search) = @_;

    my $out;

    # Split the info string on the ';' first, then extract the field of interest by splitting on equals
    for my $field (split(/;/,$info))
    {
    my ($key,$value) = split(/=/,$field);
    if ( $key eq $search ) { return $value; }
    }

    # Field not found, return 0
    return 0;
}

1;

