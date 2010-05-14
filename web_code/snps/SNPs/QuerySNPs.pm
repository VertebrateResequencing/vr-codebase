#
# Author:    Petr Danecek (pd3@sanger.ac.uk)    Team 145
#
#--------------- QuerySNPs ------------------------------------------
#
# Builds and executes the SQL query for SNPs.
#

package SNPs::QuerySNPs;

use strict;
use warnings;

use SNPs::Session;
use POSIX;
use Storable qw(freeze thaw);

sub new
{
    my ($class,$args) = @_;

    my $self = $args ? $args : {};
    if ( !exists($$self{'dbh'}) )
    {
        if ( !exists $$self{db_info} ) { die("Expected db_info parameter.\n"); }
        if ( !exists $$self{db_user} ) { die("Expected db_user parameter.\n"); }
        if ( !exists $$self{db_pass} ) { die("Expected db_pass parameter.\n"); }

        my $dbh = DBI->connect($$self{db_info},$$self{db_user},$$self{db_pass}, {'RaiseError' => 0, 'PrintError'=>0});
        if ( !$dbh || $DBI::err ) { die(sprintf("DB connection failed: %s\n",$DBI::errstr,"\n")); }
        $$self{'dbh'} = $dbh;
    }
    if ( !exists($$self{'writer'}) ) { die "Missing the 'writer' option.\n"; }
    bless $self, ref($class) || $class;

    if ( !$$self{'session'} ) 
    { 
        my $args = { 'sw'=>$$self{'writer'}->{'sw'}, 'prefix'=>$$self{store_key}, 'hours'=>$$self{store_hours} };

        my $cache_id = $$self{writer}{cgi}->param('cache');
        if ( $cache_id ) { $$args{id}=$cache_id; }

        $$self{'session'} = SNPs::Session->new($args); 
    }
    if ( !$self->cache_exists() ) 
    { 
        $self->validate_cgi_params(); 
    }

    $$self{title}        = 'SNPs' unless exists($$self{'title'});
    $$self{print_legend} = 1 unless exists($$self{'print_legend'});
    $$self{'action'}     = 'snps';

    # Why was this here?
    #   $$self{writer}->get_cookies('rows');
    #   $$self{'cache_page_size'} = $self->validate_int($$self{writer}->param('rows'));
    #   $$self{writer}->set_cookies('rows');

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


sub mysql_query
{
    my ($self,$query) = @_;
    my $sth = $$self{'dbh'}->prepare($query);
    if ( !$sth ) { die "$query: $DBI::errstr\n"; }

    # Unfortunately, there is no way how to set query timeout.
    #   We can kill the mysql connection, call $sth->cancel etc,
    #   but the query is still continues on the mysql server.
    #   To connect again and run SHOW PROCESSLIST followed by KILL 
    #   is unfortunately not an option, because of the load sharing 
    #   server design - there is no garantuee that we will reconnect 
    #   to the same server.
    #
    $sth->execute or die "$query: $DBI::errstr\n";
    return $sth;
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


sub chr_pos_to_chrpos
{
    my ($self,$chr,$pos) = @_;
    my $max_chr_len = 200000000;
    if ( !($chr=~/^\d+$/) )
    {
        if ( $chr eq 'X' ) { $chr=20; }
        elsif ( $chr eq 'Y' ) { $chr=21; }
        else { $chr=22; }
    }
    if ( $pos >= 200000000 ) { $pos=200000000; }

    my $chrpos = ($chr-1)*$max_chr_len + $pos;
    return $chrpos;
}


# The following locations are valid:
#   1:10000-20000
#   1 : 1,000,000 -2,000,000
#   Cops5
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
        elsif ($chrm ne 'X' && $chrm ne 'Y') { $chrm=0; }

        if ( !$chrm ) { die("Sorry, only the chromosomes 1-19,X,Y are available.\n"); }

        if ( $to < $from )
        {
            die ("Is this a typo? The end coordinate ($to) is smaller than the start coordinate ($from).\n");
        }

        if ( $to - $from > 10_000_000 )
        {
            die ("We are sorry, but the selected region is too big. The limit of the
                    web interface is 10Mb. Please contact us if you need bigger 
                    amounts of data.\n");
        }

        $$self{chrm}    = $chrm;
        $$self{from}    = $from;
        $$self{to}      = $to;
        my $chrpos_from = $self->chr_pos_to_chrpos($chrm,$from);
        my $chrpos_to   = $self->chr_pos_to_chrpos($chrm,$to);

        return "s.chrpos>=$chrpos_from AND s.chrpos<=$chrpos_to";
    }
    elsif ( !($loc=~/^[A-Za-z0-9.\-_]+$/) )
    {
        # Security, no tricks with altering the SQL commands
        die(qq[Sorry, could not parse the region '$location'.\n]);
    }

    $$self{gene} = $loc;
    return "c.gene=g.id AND g.value='$loc'";
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

    my $sql;
    if ( @selected )
    {
        $sql = q[s.strain=r.id AND r.value IN ('] . join(q[','], @selected) . q[')];
    }
    else
    {
        @selected = @{$self->{'strains'}};
    }

    my $strains={};
    my $i=0;
    for my $key (sort @selected) { $$strains{$key} = $i++; }

    return ($strains,$sql);
}


# Similar to validate_strains, but the selection can be void.
#
sub validate_conseq_type
{
    my ($self,$list) = @_;

    my $sql_filter = '';
    my @selected   = ();
    my $intergenic = 0;
    for my $type (@$list)
    {
        if ( $type eq 'INTERGENIC' ) { $intergenic=1; next; }
        # Security, no tricks with altering the SQL commands
        if ( !($type=~/^[A-Za-z0-9_-]+$/) ) { die("Could not parse the consequence type: [$type]\n"); }
        push @selected, $type;
    }

    if ( @selected )
    {
        $sql_filter = q[c.consequence IN ('] . join(q[','], @selected) . q[')];
        if ( $intergenic ) { $sql_filter = "($sql_filter OR c.consequence IS NULL)"; }
    }
    elsif ( $intergenic )
    { 
        $sql_filter = "c.consequence IS NULL";
    }
    return $sql_filter;
}


sub validate_cgi_params
{
    my ($self) = @_;

    my $cgi = $$self{'writer'}->{'cgi'};

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

    $$self{'writer'}->set_cookie('strain',join(',',@strains));
    $$self{'writer'}->set_cookie('conseqs',join(',',@conseqs));
    $$self{'writer'}->set_cookies('location','snp_qual','map_qual','cons_qual','depth_qual','action','depth_min','depth_max','rows');

    $$self{cache_page_size} = $self->validate_int($$self{writer}{cgi}->param('rows'));
    ($$self{'selected_strains'},$$self{'strains_sql'}) = $self->validate_strains(\@strains);
    $$self{'conseq_sql'} = $self->validate_conseq_type(\@conseqs);
    $$self{'pos_sql'}    = $self->validate_location($cgi->param('location'));
    $$self{'snp_qual'}   = $self->validate_int($cgi->param('snp_qual'));
    $$self{'map_qual'}   = $self->validate_int($cgi->param('map_qual'));
    $$self{'cons_qual'}  = $self->validate_int($cgi->param('cons_qual'));
    $$self{'depth_min'}  = $self->validate_int($cgi->param('depth_min'));
    $$self{'depth_max'}  = $self->validate_int($cgi->param('depth_max'));

    if ( !$$self{'pos_sql'} && !$$self{'gene_sql'} ) { die "Neither region nor gene selected!\n"; }


    # Save the parameters for display, both on the web page and the downloadable file.
    my @display = ();
    if ( $$self{'pos_sql'} ) { push @display, ['Location', $cgi->param('location')]; }
    if ( $$self{'conseq_sql'} ) 
    { 
        push @display, ['Consequence', \@conseqs]; 
    }
    if ( $$self{'snp_qual'} ) { push @display, ['SNP quality', $$self{'snp_qual'}]; }
    if ( $$self{'map_qual'} ) { push @display, ['Mapping quality', $$self{'map_qual'}]; }
    if ( $$self{'cons_qual'} ) { push @display, ['Consensus quality', $$self{'cons_qual'}]; }
    if ( $$self{'depth_min'} ) { push @display, ['Min depth', $$self{'depth_min'}]; }
    if ( $$self{'depth_max'} ) { push @display, ['Max depth', $$self{'depth_max'}]; }

    my ($dload_params,$html_params);
    for my $param (@display)
    {
        my $key   = $$param[0];
        my $value = $$param[1];
        $html_params .= qq[<tr><td>$key</td><td>...</td><td>];
        $html_params .= (ref($value) eq 'ARRAY' ) ? join('<br>', @$value) : $value;
        $html_params .= qq[</td></tr>];

        $dload_params .= qq[#\t$key .. ];
        $dload_params .= (ref($value) eq 'ARRAY' ) ? join(',', @$value) : $value;
        $dload_params .= "\n";
    }
    if ( $html_params )
    {
        $html_params  = q[<table id="params">] . $html_params . q[</table>];
        $dload_params = qq[# The filters used:\n] . $dload_params . "#\n";
        $$self{display_html_params}  = $html_params;
        $$self{display_dload_params} = $dload_params;
    }

    return;
}


sub sql_filters
{
    my ($self) = @_;

    my @filters = ();
    if ( $$self{'strains_sql'} ) { push @filters, $$self{'strains_sql'}; }
    if ( $$self{'pos_sql'} ) { push @filters,$$self{'pos_sql'}; }
    if ( $$self{'gene_sql'} ) { push @filters,$$self{'gene_sql'}; }
    if ( $$self{'snp_qual'} ) { push @filters, "(s.snp_qual>=$$self{snp_qual} OR s.snp_qual IS NULL)"; }
    if ( $$self{'map_qual'} ) { push @filters, "(s.map_qual>=$$self{map_qual} OR s.map_qual IS NULL)"; }
    if ( $$self{'cons_qual'} ) { push @filters, "(s.cons_qual>=$$self{cons_qual} OR s.cons_qual IS NULL)"; }
    if ( $$self{'depth_min'} ) { push @filters, "(s.depth>=$$self{depth_min} OR s.depth IS NULL)"; }
    if ( $$self{'depth_max'} ) { push @filters, "(s.depth<=$$self{depth_max} OR s.depth IS NULL)"; }
    if ( $$self{'conseq_sql'} ) { push @filters, $$self{'conseq_sql'}; }

    push @filters, 'r.id=s.strain';
    push @filters, 'chr.id=s.chr';

    return \@filters;
}


sub sql_query
{
    my ($self,$filters) = @_;
    
    my $query = qq[
        SELECT 
            r.value AS strain_name,
            chr.value AS chr,
            s.pos, s.ref_base, s.snp_base1, s.snp_base2, s.depth, s.snp_qual, s.map_qual, s.cons_qual, s.verified,
            g.value AS gene_name,
            c.*
        FROM 
            strains r, chroms chr, snps s
            LEFT JOIN snp_conseqs c ON s.snp_id=c.snp_id
            LEFT JOIN genes g ON c.gene=g.id
        WHERE
            ] .join("\n AND ",@$filters). q[ ORDER BY s.chr ASC, s.pos ASC];

    return $query;
}


sub print_legend
{
    my ($self) = @_;

    if ( !$$self{print_legend} ) { return; }

    # TODO: position taken from the first result on the page??
    my $html = $$self{'writer'};
    my $lookseq_params = 
        "show=$$self{chrm}:$$self{from}-$$self{from},paired_pileup" .
        "&amp;win=$$self{lseq_win}" .
        "&amp;width=$$self{lseq_width}" .
        "&amp;$$self{lseq_display}" .
        "&amp;lane=129S1_SvImJ.bam";

    $html->out(qq[<div id="legend">
            <b>Try:</b>
                <div style="margin-left:1em;"> <a href="$$self{'myself'}">New search</a> </div>
                <div style="margin-left:1em;"> <a href="http://www.ensembl.org/Mus_musculus/Location/View?db=core;r=$$self{chrm}:$$self{from}-$$self{to};">View in Ensembl</a> </div>
                <div style="margin-left:1em;"> <a href="$$self{'lookseq'}?$lookseq_params">View in LookSeq</a> </div>
                <div style="margin-left:1em;"> Click on SNPs for details</div>
            <div style="padding-top:1em;"><b>Download:</b></div>
        ]);

    my $session = $$self{session}->id();
    $html->out(qq[<div style="margin-left:1em;"> <a href="$$self{'myself'}?cache=$session&action=$$self{action}_dload">Results</a>\n </div>]);

    $html->out(qq[
        <div style="padding-top:1em;">
        <b>Legend:</b>
        <table>
        ]);
    for my $key (sort keys %{$$self{'conseqs'}})
    {
        my $conseq = $$self{'conseqs'}{$key};
        $html->out(qq[
            <tr><td><img class="$$conseq{'style'}" width="10" height="10" alt="" src="/icons/blank.gif"></td>
                <td>$$conseq{'label'}</td>\n]);
    }
    $html->out(qq[<tr><td><sup>+</sup></td><td>Multiple consequences</td></tr>]);
    $html->out(qq[</table></div>]);
    if ( $$self{display_html_params} )
    {
        $html->out(qq[
            <div style="padding-top:1em;position:relative;">
            <span class="button" style="font-weight:bold;" onclick="toggle_element('#params')">Filters</span>
            <div class="hidden" id="params" style="right:0em;">
                <div style="float:right;cursor:pointer;padding:0px;margin-top:-0.5em; margin-right:-0.5em;" onclick="hide_element('#params')">[x]</div>
                $$self{display_html_params}
            </div>
            </div>]);
    }
    $html->out(qq[</div>]);
}


sub print_no_match
{
    my ($self) = @_;
    die("Sorry, no matching SNPs.\n");
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
        # $session->append("\t$str");
    }
    $html->out("</tr></thead><tbody>\n");
    # $session->append("\n");
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

    my ($pos,$chr,$base,$gene_name,$gene_id);

    my $ncols = @$row;
    for (my $i=0; $i<$ncols; $i++)
    {
        if ( !$$row[$i] || !$$row[$i]->{'pos'} ) { next; }

        $pos  = $$row[$i]->{'pos'};
        $chr  = $$row[$i]->{'chr'};
        $base = $$row[$i]->{'ref_base'};
        if ( exists($$row[$i]->{'_conseqs'}) )
        {
            my $conseqs = $$row[$i]->{'_conseqs'};
            $gene_id   = $$conseqs[0]->{'gene_ensid'} || '';
            $gene_name = $$conseqs[0]->{'gene_name'} || $gene_id;
            if ( $gene_id && $pos ) { last; }
        }
    }

    return ($pos,$chr,$base,$gene_name,$gene_id);
}


sub print_row
{
    my ($self,$row) = @_;

    my ($pos,$chr,$base,$gene_name,$gene_id) = $self->nonzero_column_data($row);

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
    $style = $styles[$style];

    my $html = $$self{writer};
    $html->out("<tr><td $style>");
    if ( $print_gene ) { $html->out(sprintf qq[<a href="http://www.ensembl.org/Mus_musculus/Gene/Summary?g=ENSMUSG%.11d">$gene_name</a>], $gene_id); }
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
        for my $cons (@{$$row[$i]->{'_conseqs'}})
        {
            my $type = $$cons{'consequence'};

            if ( !$type || $type eq 'SPLICE_SITE' ) { next }  # ignore these - according to Dave these are rubbish

            if ( exists($$reported{$type}) ) { next } # report each type only once
                $$reported{$type} = 1;

            # UTR consequence types have lower priority, but we still want to keep it if there is nothing else
            #
            if ( $conseq_type && $type =~ /UTR$/ ) { next }
            $conseq_type = $type;
            if ( $$cons{'ref_codon'} ) { $codons .= '<tr><td>Ref codon:</td><td>' .$$cons{'ref_codon'}. '</td></tr>'; }
            if ( $$cons{'snp_codon1'} ) { $codons .= '<tr><td>SNP codon:</td><td>' .$$cons{'snp_codon1'}. '</td></tr>'; }
            if ( $$cons{'ref_aa'} ) { $codons .= '<tr><td>Ref aa:</td><td>' .$$cons{'ref_aa'}. '</td></tr>'; }
            if ( $$cons{'snp_aa1'} ) { $codons .= '<tr><td>SNP aa:</td><td>' .$$cons{'snp_aa1'}. '</td></tr>'; }
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
            $style   = $$self{'conseqs'}{$conseq_type}{'style'};
            $details = '<tr><td>Type:</td><td>' . join('<br />',@types) . '</td></tr>';
        }

        if ( $codons ) { $details .= $codons; }

        my $onclick='';
        my $div='';
        my $snp1 = '-';
        if ( $$row[$i]->{'depth'} )
        {
            # If we are here, there is a SNP, not an empty row.
            #
            $details = "<tr><td>Location:</td><td>$chr:$pos</td></tr>" . $details;
            $details = "<tr><td>Strain:</td><td>" .$$row[$i]->{strain_name}. "</td></tr>" . $details;
            $details .= '<tr><td>Depth:</td><td> ' . $$row[$i]->{'depth'} . '</td></tr>';
            $details .= '<tr><td>Quality:</td><td> ' . $$row[$i]->{'snp_qual'} .'/'. $$row[$i]->{'map_qual'} .'/'. $$row[$i]->{'cons_qual'} 
                            . '<br><span style="font-size:xx-small">(SNP/map/cons)</span></td></tr>';
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

            if ( !$$row[$i]->{'snp_base1'} ) { die "FIXME: yes, it can happen\n"; }
            $snp1 = $$row[$i]->{'snp_base1'};
            if ( $$row[$i]->{'snp_base2'} )
            {
                $snp1 .= '/' . $$row[$i]->{'snp_base2'};
            }
            if ( keys %$reported > 1 ) { $snp1 .= '<sup>+</sup>'; }

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
    $$self{chrm} = $$data{chrm};
    $$self{from} = $$data{from};
    $$self{to}   = $$data{to};
    $$self{display_html_params}  = $$data{display_html_params};
    $$self{display_dload_params} = $$data{display_dload_params};

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


# Reads the mySQL result and collects all rows, which belong to one
#   position. It assumes that the results are ordered by chrom and pos
#   in the ascendent order.
# The parameters are:
#   result     .. mysql result
#   buffer     .. contains a row left from previous call of get_next_position
#
sub get_next_position
{
    my ($self,$result,$buffer) = @_;

    my $strains = $$self{'selected_strains'};
    my ($row,$prev_pos,$prev_chr,$idx);
    my @data = ();
    while (1)
    {
        $row = @$buffer ? shift(@$buffer) : $result->fetchrow_hashref;
        if ( !$row ) { last }

        my $strain = $$row{'strain_name'};
        my $chr    = $$row{'chr'};
        my $pos    = $$row{'pos'};

        # If the search was by gene and not by region, the values of chrm:pos
        #   are not set. They will be required later when printing the legend
        #   (a link to LookSeq).
        if ( !$$self{chrm} ) { $$self{chrm}=$chr; }
        if ( !$$self{from} ) { $$self{from}=$pos; }

        if ( $prev_pos && ($chr ne $prev_chr || $pos ne $prev_pos) )
        {
            push @$buffer, $row;
            last;
        }

        $idx = $$strains{$strain};
        if ( !exists($data[$idx]) )
        {
            $data[$idx] = $row;
            $data[$idx]->{'_conseqs'} = [];
        }

        $data[$idx]->{strain_name} = $strain;
        if ( $$self{mouseinfo}{$strain}{bam} )
        {
            $data[$idx]->{bam} = $$self{mouseinfo}{$strain}{bam};
        }

        # Rather than copying the consequence-specific data, save the whole row.
        push @{$data[$idx]->{'_conseqs'}}, $row;

        $prev_pos   = $pos;
        $prev_chr   = $chr;
    }

    if ( $prev_pos ) { $$self{to} = $prev_pos; }

    if ( !@data ) { return 0; }
    return \@data;
}

sub run
{
    my ($self,@args) = @_;

    my $result;
    if ( !$self->cache_exists() )
    {
        # If the cache does not exist, build the SQL query
        my $filters = $self->sql_filters();
        my $query   = $self->sql_query($filters);

        # $$self{'writer'}->out(qq[<div><pre>$query</pre></div>]); return;

        $result = $self->mysql_query($query);
        $$self{'nrows'} = $result->rows;    # Beware: this is not the number of rows in the result table!
        if ( !$$self{'nrows'} )
        {
            $self->print_no_match();
            return;
        }
    }

    $$self{'writer'}->out(qq[<h3>Mouse Genomes Project - $$self{'title'}</h3>
        <div style="font-size:x-small;">(NOTE: All SNPs/indels reported have not been experimentally validated)</div>
    ]);
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
        my $irow = 0;
        my $buffer = [];
        while (my $pos=$self->get_next_position($result,$buffer))
        {
            push @{$store{data}}, $pos;

            if ( $irow<$$self{cache_page_size} )
            {
                $self->print_row($pos);
            }
            $irow++;
        }
        # These are necessary values for the cached pages display.
        $store{selected_strains} = $$self{selected_strains};
        $store{chrm} = $$self{chrm};
        $store{from} = $$self{from};
        $store{to}   = $$self{to};
        $store{display_html_params}  = $$self{display_html_params};
        $store{display_dload_params} = $$self{display_dload_params};

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

    #   if ( $$self{cached_data} )
    #   {
    #       use Data::Dumper;
    #       $$self{writer}->out(qq[<div><pre>] . Dumper($$self{cached_data}) . qq[</pre></div>]);
    #   }

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

