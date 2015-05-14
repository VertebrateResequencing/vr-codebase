#
# Author:    Petr Danecek (pd3@sanger.ac.uk)    Team 145
# Modified:		John Maslen  (jm23@sanger.ac.uk)   Team 145
#
#--------------- QuerySNPs ------------------------------------------
#
# Builds and executes the TABIX command to retrieve SNP data via FTP
#

package SNPs::QuerySNPs;

use strict;
use warnings;

use SNPs::Session;
use SNPs::Writer;
use POSIX;
use DBI;
use Storable qw(freeze thaw);
use CGI::Carp qw(fatalsToBrowser);
use Data::Dumper;

sub new
{
    my ($class,$args) = @_;

    my $self = $args ? $args : {};
    if ( !exists($$self{'writer'}) ) { die "Missing the 'writer' option.\n"; }
    bless $self, ref($class) || $class;

    if ( !$$self{'session'} ) 
    { 
        my $args = { 'sw'=>$$self{'writer'}->{'sw'}, 'prefix'=>$$self{store_key}, 'hours'=>$$self{store_hours} };

        my $cache_id = $$self{writer}{cgi}->param('cache');
        if ( $cache_id ) { $$args{id}=$cache_id; }

        $$self{'session'} = SNPs::Session->new($args); 
    }
    
	#$$self{'writer'}->error_exit("Service Temporarily Unavailable ;;; The Mouse Genome Project SNP/Indel/SV viewer is currently being updated. We apologise for any inconvenience caused.\n") unless $$self{'available'};    
    
    if ( !$self->cache_exists() ) 
    { 
        $self->validate_cgi_params(); 
    }

    $$self{title}        = 'SNPs' unless exists($$self{'title'});
    $$self{print_legend} = 1 unless exists($$self{'print_legend'});
    $$self{'action'}     = 'snps';
    
    # The cache_page_size can be overriden by pg_size - user can browse his old results with different 
    #   page size than the new ones.
    if ( $$self{writer}{cgi}->param('pg_size') )
    {
        $$self{cache_page_size} = $self->validate_int($$self{writer}{cgi}->param('pg_size'));
    }
    if ( !$$self{cache_page_size} ) 
    { 
        # We should never get here, unless the user builds the URL himself.
        $$self{cache_page_size}=50; 
    }
    return $self;
}

sub validate_int
{
    my ($self,$value) = @_;
    if ( !$value ) { return ''; }
    $value =~ s/^\s*//;
    $value =~ s/\s*$//;
    if ( !($value=~/^\d+$/) ) { die ("Could not parse the number \"$value\".\n"); }
    return $value;
}


# The following locations are valid:
#   1:10000-20000
#   1 : 1,000,000 -2,000,000
#
sub validate_location
{
    my ($self,$location) = @_;

    if ( !$location ) { return ''; }

    my $loc = $location;
    $loc =~ s/,//g;
    $loc =~ s/\s+//g;

    if ( $loc=~/^([A-Za-z0-9.\-_]+):(\d+)-(\d+)$/i )
    {
        my $chrm = uc($1);
        my $from = $2;
        my $to   = $3;

        if ( $chrm=~/^\d+$/ )
        {
            if ($chrm<1 || $chrm>19) { $chrm=0; }
        }
        elsif ($chrm ne 'X') { $chrm=0; }

        if ( !$chrm ) { die("Chromosome not available ;;; Chromosomes 1-19 and X are available for searching.\n"); }

        if ( $to < $from )
        {
            die ("Error in search coordinates ;;; The end coordinate ($to) is smaller than the start coordinate ($from).\n");
        }

        if ( $to - $from > 20_000_000 )
        {
            die (qq[Search region too large ;;; The current maximum region search limit is 20Mb. For larger queries - the raw data can be obtained from our <a href=&quot;ftp://ftp-mouse.sanger.ac.uk/&quot;>ftp site</a>.\n]);
        }

        $$self{chrm}    = $chrm;
        $$self{from}    = $from;
        $$self{to}      = $to;
        $$self{loc}		= $chrm.':'.$from.'-'.$to;	
        return "$chrm:$from-$to";
    }
    elsif ( $loc=~/^[A-Za-z0-9.\-_]+$/ )
    {
        my $dbfile = $$self{'gene_db_location'};
		my $dbh = DBI->connect(
			"dbi:SQLite:dbname=$dbfile", # DSN: dbi, driver, database file
			"",                          # no user
			"",                          # no password
			{ RaiseError => 1 },         # complain if something goes wrong
		) or die "ERROR\n".$DBI::errstr;
		my $select_sql = "select chromosome, from_pos, to_pos, ensembl_id from gene_position where lower(gene_name) = ?";
		my $sth = $dbh->prepare($select_sql);
  		my $location;
		my ($chrm, $from, $to, $ensembl);
		$sth->execute(lc($loc)) 
			or die "Couldn't execute statement: " . $sth->errstr;
		$sth->bind_columns(undef, \$chrm, \$from, \$to, \$ensembl);
  		while ($sth->fetch) {
  	    	$location = $chrm.':'.$from.'-'.$to;
  	    	$$self{chrm}    = $chrm;
        	$$self{from}    = $from;
        	$$self{to}      = $to;
        	$$self{ens_id}  = $ensembl;
        	$$self{loc}		= $location;	
		}	
		$sth->finish;
    	$dbh->disconnect;
    	$location ? return $location : $$self{'writer'}->error_exit("Gene name not found ;;; No chromosome location could be found for $loc.\n");
    }
	else
    {
        die(qq[Sorry, could not parse the region '$location'.\n]);
    }
}


# Validate strains:
#   1) if none strain selected, all will be selected.
#   2) if 'All' is selected but also some other strain is selected, all will be selected.
#
# Returned is a a hash containing the selected strains and a SQL query fragment.
#
sub validate_strains
{
    my ($self,$list) = @_;

    my @selected = ();
    for my $str (@$list)
    {
        if ( $str eq 'All' ) { @selected=(); last; }
        if ( !exists($$self{mouseinfo}{$str}) ) { next }
        push @selected, $str;
    }

    @selected = @{$self->{'strains'}} unless @selected;

    my $strains={};
    my $i=0;
    for my $key (sort @selected) { 
    	$$strains{$key} = $i;
    	$i++; }

    return ($strains);
}



# Similar to validate_strains, but the selection can be void.
#
sub validate_conseq_type
{
    my ($self,$list) = @_;

    my $sql_filter = '';
    my %selected;
    my $intergenic = 0;
    for my $type (@$list)
    {
        if ( !($type=~/^[A-Za-z0-9_-]+$/) ) { die("Could not parse the consequence type: [$type]\n"); }
        $selected{$type} = 1;
    }

    return \%selected;
}

sub validate_sv_type
{
    my ($self,$list) = @_;

    my %selected;
    for my $type (@$list)
    {
        if ( !($type=~/^[A-Za-z0-9_-]+$/) ) { die("Could not parse the sv type: [$type]\n"); }
        $selected{$type} = 1;
    }

    return \%selected;
}

sub validate_cgi_params
{
    my ($self) = @_;

    my $cgi = $$self{'writer'}->{'cgi'};
	
	$$self{'writer'}->error_exit("Service Temporarily Unavailable ;;; The Mouse Genome Project SNP/Indel/SV viewer is currently being updated. We apologise for any inconvenience caused.\n") unless $$self{'available'};
    
    my @strains = ();
    for my $strain (keys %{$$self{mouseinfo}})
    {
        my $id = $strain;
        $id =~ s{/}{_}g;
        if ( $cgi->param('str_'.$strain) ) { push @strains, $strain; }
    }
	
    my @conseqs = ();
    for my $conseq (keys %{$$self{conseqs}})
    {
        if ( $cgi->param('cnsq_'.$conseq) ) { push @conseqs, $conseq; }
    }
		
    my @svtypes = ();
    for my $svtype (keys %{$$self{sv_types}})
    {
        if ( $cgi->param('svt_'.$svtype) ) { push @svtypes, $svtype; }
    }
		
    $$self{'writer'}->set_cookie('strain',join(',',@strains));
    $$self{'writer'}->set_cookie('conseqs',join(',',@conseqs));
    $$self{'writer'}->set_cookie('svtypes',join(',',@svtypes));
    $$self{'writer'}->set_cookies('location','action','rows');

    $$self{cache_page_size} = $self->validate_int($$self{writer}{cgi}->param('rows'));
    ($$self{'selected_strains'}) = $self->validate_strains(\@strains);
    $$self{'selected_conseq'} = $self->validate_conseq_type(\@conseqs);
    $$self{'selected_svtypes'} = $self->validate_sv_type(\@svtypes);
    $$self{'pos_str'}    = $self->validate_location($cgi->param('location'));

    if ( !$$self{'pos_str'} ) { die "Chromosome location must be a region (format = 1:10000000-10040000).\n"; }


    # Save the parameters for display, both on the web page and the downloadable file.
    my @display = ();
    my @svdisplay = ();
    if ( $$self{'pos_str'} ) { 
    	push @display, ['Location', $cgi->param('location')]; 
    	push @svdisplay, ['Location', $cgi->param('location')]; 
    }
    if ( $$self{'selected_conseq'} ) 
    { 
        push @display, ['Consequence', \@conseqs]; 
    }
    if ( $$self{'selected_svtypes'} ) 
    { 
        push @svdisplay, ['SV types', \@svtypes]; 
    }
    
    my ($dload_params,$html_params,$svdload_params,$svhtml_params);
    for my $param (@display)
    {
        my $key   = $$param[0];
        my $value = $$param[1];
        $html_params .= qq[<tr><td style="text-align: left;"><b><u>$key:</u></b></td></tr><tr><td style="padding-left: 2em;">];
        $html_params .= (ref($value) eq 'ARRAY' ) ? join('<br>', @$value) : $value;
        $html_params .= qq[</td></tr>];

        $dload_params .= qq[#\t$key .. ];
        $dload_params .= (ref($value) eq 'ARRAY' ) ? join(',', @$value) : $value;
        $dload_params .= "\n";
    }
    for my $param (@svdisplay)
    {
        my $key   = $$param[0];
        my $value = $$param[1];
        $svhtml_params .= qq[<tr><td style="text-align: left;"><b><u>$key:</u></b></td></tr><tr><td style="padding-left: 1em;">];
        $svhtml_params .= (ref($value) eq 'ARRAY' ) ? join('<br>', @$value) : $value;
        $svhtml_params .= qq[</td></tr>];

        $svdload_params .= qq[#\t$key .. ];
        $svdload_params .= (ref($value) eq 'ARRAY' ) ? join(',', @$value) : $value;
        $svdload_params .= "\n";
    }    
    if ( $html_params )
    {
        $html_params  = q[<table id="params">] . $html_params . q[</table>];
        $dload_params = qq[# The filters used:\n] . $dload_params . "#\n";
        $$self{display_html_params}  = $html_params;
        $$self{display_dload_params} = $dload_params;
    }
    if ( $svhtml_params )
    {
        $svhtml_params  = q[<table id="params">] . $svhtml_params . q[</table>];
        $svdload_params = qq[# The filters used:\n] . $svdload_params . "#\n";
        $$self{display_html_svparams}  = $svhtml_params;
        $$self{display_dload_svparams} = $svdload_params;
    }    
    return;
}

sub run_tabix_command
{
	my ($self) = @_;
	my $tabix_cmd = "cd $$self{'index_directory'} ; $$self{'tabix_location'} -h $$self{'snps_location'}  $$self{'pos_str'}";
	my @mouse_strains2 = ();
	my $strains = $$self{'selected_strains'};
	my @sort_strains = (sort {$$strains{$a}<=>$$strains{$b}} keys %$strains);
	my $input = `$tabix_cmd`;
	my @lines = split('\n', $input); 
	my %output = ();
	my %conseq_selected = %{$$self{'selected_conseq'}};
	
	my ($gt, $atg, $mq, $gq, $dp);
	my $linecount = 0;

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
		else {
		my @line = split('\t', $line_in);
		my @lineoutput = ();		
		if ($linecount==0) {
			($gt, $atg, $mq, $gq, $dp) = get_filter_order($line[8]);
			$linecount++;
		}		
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
    	my $allele_count = 0;
    	$alleles{$allele_count++} = $line[3];
    	foreach ( split(',', $line[4] )) {
    		$alleles{$allele_count++} = $_;	
    	}
    	die "Mouse strains in vcf file do not match expected strains" unless scalar(@strain_data) == scalar( @mouse_strains2);
    	 
    	for my $i ( 0 .. $#strain_data ) {
    		my %linehash;
    		my $high_qual_filter;
    		my $strain = $mouse_strains2[$i];   	
			my @data = split (':', $strain_data[$i]);
			for my $gene (keys %gene_conseq) {
    			my $depth = $data[$dp];
    			my $cons_qual = $data[$gq]-10;
    			my $map_qual = $data[$mq];
    			my $gt_qual = $data[$gq];
				$high_qual_filter = $data[$atg];
    			if ($high_qual_filter == 0 || $high_qual_filter == -1 || $high_qual_filter == -8) {
    				$linehash{'strain_name'} = $strain;
    				$linehash{'chr'} = $line[0];
    				$linehash{'pos'} = $line[1];
    				$linehash{'atg_qual'} = $high_qual_filter;
    			}
    			else {					
    				my %allele_hash;
    				my @allele_out;
    				foreach (split('/', $data[$gt])) {
    					if (!$allele_hash{$alleles{$_}}) {
    						$allele_hash{$alleles{$_}} = 1;
    						push @allele_out, $alleles{$_};
    					} 
    				}
    				$linehash{'strain_name'} = $strain;
    				$linehash{'consequence'} = $gene_conseq{$gene};
    				$linehash{'gene_name'} = $gene unless $gene eq '-';
    				$linehash{'chr'} = $line[0];
    				$linehash{'pos'} = $line[1];
    				$linehash{'ref_base'} = $line[3];
    				$linehash{'snp_base'} = join('/',@allele_out);
    				$linehash{'depth'} = $depth;
    				$linehash{'map_qual'} = $map_qual;
    				$linehash{'gt_qual'} = $gt_qual;
    				$linehash{'cons_qual'} = $cons_qual;
    				$linehash{'atg_qual'} = $high_qual_filter; 
    				if ( $$self{mouseinfo}{$strain}{bam} ) {
    					$linehash{'bam'} = $$self{mouseinfo}{$strain}{bam};
        			}
    			}
    		}
    		if (%linehash) {
    			push @lineoutput, {%linehash};
    		}	
    	}
    	if (@lineoutput) {
    		my @strain_output = ();
    		my $conseq_found = 0;
    	    for my $str ( @sort_strains ) {
    	    	for my $str_res ( @lineoutput ) {
    	    		if ( $$str_res{'strain_name'} eq $str ) {
    	    			push @strain_output, $str_res;
    	    			if ( !$conseq_found && exists $$str_res{'consequence'} ) {
    	    				$conseq_found = 1;
    	    			}
    	    		}
    	    	}
    	    }
    	    if ( $conseq_found ) {
    			$output{$linecount} = [@strain_output];
    			$linecount++;
    		}	
    	}
    	}
    }
	return \%output;
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
                <div style="margin-left:1em;"> <a href="http://may2012.archive.ensembl.org/Mus_musculus/Location/View?db=core;r=$$self{chrm}:$$self{from}-$$self{to};" target="_ensembl_snp">View in Ensembl</a> </div>
                <div style="margin-left:1em;"> Click on SNPs for details</div>
            <div style="padding-top:1em;"><b>Download:</b></div>
        ]);
        
#lookseq details removed from above:
#        
#     my $lookseq_win = $$self{to}-$$self{from};
#     if ($lookseq_win > 100000) {
#     	$lookseq_win = 100000;
#     }
#     my $lookseq_params = 
#         "show=$$self{chrm}:$$self{from}-$$self{to},paired_pileup" .
#         "&amp;win=$lookseq_win" .
#         "&amp;width=$$self{lseq_width}" .
#         "&amp;$$self{lseq_display}" .
#         "&amp;lane=DBA.bam";
#link to lookseq:
#   <div style="margin-left:1em;"> <a href="$$self{'lookseq'}?$lookseq_params" target="_lookseq_snp">View in LookSeq</a> </div>
                

    my $session = $$self{session}->id();
    $html->out(qq[<div style="margin-left:1em;">Formats: <a href="$$self{'myself'}?cache=$session&action=$$self{action}_dload">tab</a> or <a href="$$self{'myself'}?cache=$session&action=$$self{action}_dload_csv">csv</a>.\n </div>]);

    $html->out(qq[
        <div style="padding-top:1em;">
        <b>Legend:</b><br/>
         &nbsp;&nbsp;<b>SNP Consequences:</b>
        <table>]);
    #for my $key (sort keys %{$$self{'conseqs'}})
    for my $key (sort keys %{$$self{'selected_conseq'}})
    {
        my $conseq = $$self{'conseqs'}{$key};
        $html->out(qq[
            <tr><td><img class="$$conseq{'style'}" width="10" height="10" alt="" src="/icons/blank.gif"></td>
                <td>$$conseq{'label'}</td>\n]);
    }
    $html->out(qq[<tr><td><sup>+</sup></td><td>Multiple consequences</td></tr>]);
    $html->out(qq[</table><br/>
            &nbsp;&nbsp;<b>Base Calling examples:</b>
        <table>
        <tr><td>T</td><td>High confidence SNP</td></tr>\n
        <tr><td>t</td><td>Low confidence SNP</td></tr>\n        
        <tr><td>-</td><td>High confidence reference</td></tr>\n
        <tr><td class="c14">-</td><td>Low confidence reference</td></tr>\n
        <tr><td class="c14"> </td><td>Genotype not called</td></tr> 
        </table></div>]);
    
    $html->out(qq[</div>]);
}


sub print_no_match
{
    my ($self) = @_;
    die("SNPs not found ;;; Sorry, unable to find SNPs for $$self{loc}.\n");
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
        ]);
    for my $str (sort {$$strains{$a}<=>$$strains{$b}} keys %$strains)
    {
        my $img_name = lc($str);
        $img_name =~ s{/}{_}g;
        $html->out(qq[<th class="results"><img src="$url_imgs/$img_name.png" alt="$str" /></th>]);
    }
    $html->out("</tr></thead><tbody>\n");
}

sub print_footer
{
    my ($self) = @_;
    $$self{'writer'}->out("</tbody></table>");
}


# Find any non-empty column to get the information common to all columns: Gene, chromosome, position,
#   the reference base.
#
sub nonzero_column_data
{
    my ($self,$row) = @_;

    my ($pos,$chr,$base,$gene_name);

    my $ncols = @$row;
    for (my $i=0; $i<$ncols; $i++)
    {
        if ( !$$row[$i] || !$$row[$i]->{'pos'} ) { next; }

        $pos  = $$row[$i]->{'pos'};
        $chr  = $$row[$i]->{'chr'};
        if ( exists($$row[$i]->{'ref_base'}) ) {
        	$base = $$row[$i]->{'ref_base'};
        }
        if ( exists($$row[$i]->{'consequence'}) )
        {
            $gene_name = $$row[$i]->{'gene_name'};
            if ( $gene_name && $pos ) { last; }
        }
    }
	
    return ($pos,$chr,$base,$gene_name);
}


sub print_row
{
    my ($self,$row) = @_;

    my ($pos,$chr,$base,$gene_name) = $self->nonzero_column_data($row);
    my $strains    = $$self{'selected_strains'};
    # Alternate background for each gene and print gene name only if different from the previous row
    #
    
	# 1. ATG = 1 (high qual ref) : can be left as '-'
	# 2. ATG = -8 (missing or het) : no dash '-'
	# 3. ATG = -1 (low qual ref) : '-' with light grey box
	# 4. all other ATGs (low qual SNPs): The SNP base in italics (and the background colored according to consequence, if any) 
    
    my $print_gene = 0;
    my $style      = 0;
    if ( !$gene_name ) { $style=-1; }
    elsif ( !$$self{'last_gene_name'} ) { $style=0; $print_gene=1; }
    elsif ( $$self{'last_gene_name'} eq $gene_name ) { $style=$$self{'last_gene_style'}; }
    else { $style = $$self{'last_gene_style'} ? 0 : 1; $print_gene=1; }

    $$self{'last_gene_name'} = $gene_name;
    $$self{'last_gene_style'} = $style;
    
    my @styles = ('class="gene1"','class="gene2"','');
    $style = $styles[$style];

    my $html = $$self{writer};
    $html->out("<tr><td $style>");
    if ( $print_gene ) { 
        if ($$self{ens_id}) {
        	my $ens_id = $$self{ens_id};
        	$html->out(sprintf qq[<a href="http://may2012.archive.ensembl.org/Mus_musculus/Gene/Summary?g=$ens_id" target="_ensembl_snp">$gene_name</a>]); 
        }
        else {
        	$html->out(sprintf qq[<a href="http://may2012.archive.ensembl.org/Mus_musculus/Gene/Summary?g=$gene_name" target="_ensembl_snp">$gene_name</a>]); 
        }	
    }
    $html->out(qq[</td><td $style>$chr</td><td $style>$pos</td><td class="resultsbar">$base</td>]);

    my $ncols = scalar keys %{$$self{'selected_strains'}};
    for (my $i=0; $i<$ncols; $i++)
    {
        # Selected from all the possible consequence types one for coloured display.
        #   and the rest will be reported in a popup div as a list.
        #
        my $reported     = {};
        my $conseq_type  = '';
        my $details      = '';
        my $codons       = '';
        my $atg_qual = $$row[$i]->{'atg_qual'};
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
                if ( !exists($$self{'conseqs'}{$key}) ) { die("FIX ME: The consequence type \"$key\" not covered.\n"); }
                push @types, qq[<span class="$$self{'conseqs'}{$key}{'style'}">$$self{'conseqs'}{$key}{'label'}</span>];
            }
            $style   = $$self{'conseqs'}{$conseq_type}{'style'};
            $details = '<tr><td>Type:</td><td>' . join('<br />',@types) . '</td></tr>';
        }
        if ($atg_qual == -1) {
        	$style = 'class="c14"';
        }

        if ( $codons ) { $details .= $codons; }

        my $onclick='';
        my $div='';
        my $snp1 = '-';
        if ($atg_qual == -8) {
        	$style = 'class="c14"';
        	$snp1 = ' ';
        }
        if ( $$row[$i]->{'ref_base'} )
        {
            # If we are here, there is a SNP, not an empty row.
            #
            $details = "<tr><td>Location:</td><td>$chr:$pos</td></tr>" . $details;
            $details = "<tr><td>Strain:</td><td>" .$$row[$i]->{strain_name}. "</td></tr>" . $details;
            $details .= '<tr><td>Depth:</td><td> ' . $$row[$i]->{'depth'} . '</td></tr>';
            $details .= '<tr><td>Quality:</td><td> ' . $$row[$i]->{'gt_qual'} .'/'. $$row[$i]->{'map_qual'} .'/'. $$row[$i]->{'cons_qual'} 
                            . '<br><span style="font-size:xx-small">(GT/map/cons)</span></td></tr>';
            $details .= qq[<tr><td colspan="2"><a href="http://may2012.archive.ensembl.org/Mus_musculus/Location/View?db=core;r=$chr:$pos-$pos;" target="_ensembl_snp">View in Ensembl</a></td></tr>];
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
            $onclick = qq[onclick="toggle_one_element('#$id')" ];
            $div = qq[
                <div class="hidden" id="$id">
                <div style="text-align:right;font-size:smaller;margin-top:-0.5em; margin-right:-0.5em;">
					<span class="button">[x]</span></div>
                <table class="details">
                    $details
                </table>
                </div>];
            if ( !$$row[$i]->{'snp_base'} ) { die "FIXME: yes, it can happen\n"; }
            $snp1 = $$row[$i]->{'snp_base'};
            if ( keys %$reported > 1 ) { $snp1 .= '<sup>+</sup>'; }
            if ($atg_qual != 1) {
            	$snp1 = lc($snp1);
            }  

            $style = qq[class="$style button"];
        }

        $html->out(qq[<td $style $onclick>$snp1$div</td>]);
    }
    $html->out("</tr>\n");
}


sub cache_exists
{
    my ($self) = @_;

    my $cache_id = $$self{writer}{cgi}->param('cache');
    if ( ! $cache_id ) { return 0; }

    if ( $$self{cached_data} ) { return 1; }

    my $data = $$self{session}->retrieve($cache_id);
    if ( !$data ) { die "Sorry, the session $cache_id expired.\n"; }

    $data = thaw($data);
    $$self{selected_strains} = $$data{selected_strains};
    $$self{selected_conseq} = $$data{selected_conseq};
    $$self{selected_svtypes} = $$data{selected_svtypes};
    $$self{chrm} = $$data{chrm};
    $$self{from} = $$data{from};
    $$self{to}   = $$data{to};
    $$self{loc}   = $$data{loc};
    $$self{display_html_params}  = $$data{display_html_params};
    $$self{display_dload_params} = $$data{display_dload_params};
    $$self{display_html_svparams}  = $$data{display_html_svparams};
    $$self{display_dload_svparams} = $$data{display_dload_svparams};

    my @rows = @{$$data{data}};
    if ( !@rows ) { return 0; }

    $$self{cached_data} = \@rows;

    $$self{cache_page} = $$self{writer}{cgi}->param('page');
    if ( !$$self{cache_page} ) { $$self{cache_page}=0; }
    $$self{icache} = 0;

    return 1; 
}

sub cache_npages
{
    my ($self) = @_;
    if ( !$$self{cached_data} ) { return 0; }
    my $npages = scalar(@{$$self{cached_data}});
    return POSIX::ceil($npages / $$self{cache_page_size});
}

sub cache_get_next
{
    my ($self) = @_;

    if ( !$$self{cached_data} ) { return 0; }
    if ( $$self{icache} >= $$self{cache_page_size} ) { return 0; }

    my $i = $$self{cache_page_size} * $$self{cache_page} + $$self{icache}++;
    my $data = $$self{cached_data};

    if ( $i >= scalar @$data ) { return 0; }

    return $$data[$i];
}

sub cache_nrows
{
    my ($self) = @_;
    if ( !$$self{cached_data} ) { return 0; }
    return scalar(@{$$self{cached_data}});
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

sub get_filter_order 
{
	my ($filter) = @_;
	my @filter_fields = split(':', $filter);
	# Get the position of the parameter in the FORMAT, e.g. GT:ATG:MQ:HCG:GQ:DP
	# as they might not necessarly always be in that order.
	# Will return an array with the position in the order GT, ATG, MQ, GQ, DP.
	my ($gt, $atg, $mq, $gq, $dp);
	
	for my $i (0 .. $#filter_fields) {
		if ($filter_fields[$i] eq 'GT') {
			$gt = $i;
		}
		elsif ($filter_fields[$i] eq 'ATG') {
			$atg = $i;
		}
		elsif ($filter_fields[$i] eq 'MQ') {
			$mq = $i;
		}
		elsif ($filter_fields[$i] eq 'GQ') {
			$gq = $i;
		}
		elsif ($filter_fields[$i] eq 'DP') {
			$dp = $i;
		}
	}
	return ($gt, $atg, $mq, $gq, $dp);	
}

sub run
{
    my ($self,@args) = @_;
    my $result;
    if ( !$self->cache_exists() )
    {
        $result = $self->run_tabix_command;
        $$self{'nrows'} = scalar %{ $result };    # Beware: this is not the number of rows in the result table!
        if ( !$$self{'nrows'} )
        {
            $self->print_no_match();
            return;
        }
    }

    $$self{'writer'}->out("<h3>Mouse Genomes Project - $$self{'title'} [$$self{'loc'}]</h3>");
    $$self{'writer'}->out(qq[<div style="font-size:x-small;">(NOTE: All SNPs,indels and SVs reported have not been experimentally validated)</div>]);
    $$self{'writer'}->out(qq[<table style="width:100%;"><tr><td><div class="resultsdiv">]);
    $self->print_header();

    my $cache_id = 0;
    my $npages   = 0;

    # We need to store the data in a cache in order to get the functionality of multiple
    #   screens per one mysql query. The data may have complex structure, how to store
    #   them? Use Storable for this.
    #
    if ( $self->cache_exists() )
    {
        while (my $pos=$self->cache_get_next())
        {
            $self->print_row($pos);
        }
        $$self{'nrows'}  = $self->cache_nrows();
        $npages   = $self->cache_npages();
    }
    else
    {
        my %store = ();

        # We are given variable number of rows for each position.
        #   The subroutine get_next_position collects all data for one position
        #   and the last item is stored in the buffer which is reused upon subsequent
        #   call.
        #
        #my $buffer = [];
        #while (my $pos=$self->get_next_position($result,$buffer))
        my @resultsort = sort { $a <=> $b } keys %{ $result };
        my $irow = 0;
        for my $row ( @resultsort )
        {
            my $pos = $$result{$row};
            
            push @{$store{data}}, $pos;

            if ( $irow<$$self{cache_page_size} )
            {
                $self->print_row($pos);
            }
            $irow++;
        }
        # These are necessary values for the cached pages display.
        $store{selected_strains} = $$self{selected_strains};
        $store{selected_conseq} = $$self{selected_conseq};
        $store{selected_svtypes} = $$self{selected_svtypes};
        $store{chrm} = $$self{chrm};
        $store{from} = $$self{from};
        $store{to}   = $$self{to};
        $store{loc}   = $$self{loc};
        $store{display_html_params}  = $$self{display_html_params};
        $store{display_dload_params} = $$self{display_dload_params};
        $store{display_html_svparams}  = $$self{display_html_svparams};
        $store{display_dload_svparams} = $$self{display_dload_svparams};       

        $$self{'session'}->append(freeze(\%store));
        $$self{'session'}->store();
        $npages = POSIX::ceil($irow/$$self{cache_page_size});
    }
    $self->print_footer();

    # Add links to next and prev result pages.
    if ( $npages )
    {
        # The desired behaviour for 10 pages and $half_win=2
        #
        #       0123456789
        #       x1234.   >
        #       0x234.   >
        #       01x34.   >
        #       012x45.  >
        #       <.23x56. >
        #       < .34x67.>
        #       <   45x789
        #       <   456x89
        #       <   4567x9
        #       <   45678x
        #
        $$self{writer}->out(q[<div style="padding-top:1em;text-align:center;">]);

        my $cache_id = $$self{session}->id();
        my $ipage = $$self{cache_page} ? $$self{cache_page} : 0;
        my $half_win = 3;
        my $pg_size = "pg_size=$$self{cache_page_size}";

        my ($ipage_from,$ipage_to);
        if ( $ipage<=$half_win )
        {
            $ipage_from = 0;
            $ipage_to   = $half_win*2;
        }
        elsif ( $npages-$half_win-1 <= $ipage )
        {
            $ipage_from = $npages - $half_win*2 -1;
            $ipage_to   = $npages - 1;
        }
        else
        {
            $ipage_from = $ipage - $half_win;
            $ipage_to   = $ipage + $half_win;
        }
        if ( $ipage_from<=1 ) { $ipage_from=0; }
        if ( $ipage_to>=$npages-2 ) { $ipage_to=$npages-1; }

        if ( $ipage_from>0 ) 
        { 
            $$self{'writer'}->out(qq[<a href="$$self{'myself'}?action=$$self{action}&amp;cache=$cache_id&amp;page=0&amp;$pg_size">&lt;&lt;</a>&nbsp;&nbsp;]);
            $$self{'writer'}->out('&nbsp;...'); 
        }
        for (my $i=$ipage_from; $i<=$ipage_to; $i++)
        {
            if ( $ipage == $i )
            {
                $$self{'writer'}->out('&nbsp;' . ($i+1));
            }
            else
            {
                $$self{'writer'}->out(qq[&nbsp;&nbsp;<a href="$$self{'myself'}?action=$$self{action}&amp;cache=$cache_id&amp;page=$i&amp;$pg_size">].($i+1).qq[</a>]);
            }
        }
        if ( $ipage_to<$npages-1 ) 
        { 
            $$self{'writer'}->out('&nbsp;...'); 
            $$self{'writer'}->out(qq[&nbsp;&nbsp;<a href="$$self{'myself'}?action=$$self{action}&amp;cache=$cache_id&amp;page=].($npages-1).qq[&amp;$pg_size">&gt;&gt;</a>]);
        }
        $$self{writer}->out(q[</div>]);
    }

    $$self{'writer'}->out("</div></td><td>");
    $self->print_legend();
    $$self{'writer'}->out("</td></tr></table>");

    if ( !$$self{'nrows'} )
    {
        $$self{writer}->print_form();
    }
    $$self{writer}->print_footer();
    return;
}

1;

