#
# Author:    Petr Danecek (pd3@sanger.ac.uk)    Team 145
# Modified:		John Maslen  (jm23@sanger.ac.uk)   Team 145
#
#--------------- Writer ----------------------------------------------------
#
# The basic HTML output with Sanger header and footer.
#

package SNPs::Writer;

use strict;
use warnings;

use SangerWeb;

sub new
{
    my ($class,$args) = @_;

    my $self = $args ? $args : {};
    if ( !exists($$self{'sw'}) )
    {
        my $opts =  {
            'title'      => q(Mouse Genomes Project),
            'banner'     => q(),
            'inifile'    => SangerWeb->document_root() . q(/modelorgs/mousegenomes/snps-header.ini),
            'stylesheet' => $$self{'css'},
            'jsfile'     => $$self{'jsfile'},
        };
        $$self{'sw'}  = SangerWeb->new($opts);
        $$self{'cgi'} = $$self{'sw'}->cgi();
    }
    $$self{strains} = [ sort keys %{$$self{mouseinfo}} ];

    bless $self, ref($class) || $class;
    return $self;
}

sub param
{
    my ($self,@args) = @_;
    return $$self{'cgi'}->param(@args);
}

sub get_cookie
{
    my ($self,$key) = @_;

    return $$self{'cgi'}->cookie($key);
}

sub set_cookie
{
    my ($self,$key,$value) = @_;

    if ( $$self{'header_printed'} ) { $self->error_exit("Header already sent, cannot set cookie $key=$value\n"); }
    $$self{'set_cookies'}{$key} = $value;
    return;
}

sub unset_cookie
{
    my ($self,$key) = @_;

    if ( $$self{'header_printed'} ) { $self->error_exit("Header already sent, cannot unset cookie $key\n"); }
    push @{$$self{'unset_cookies'}}, $key;
    return;
}

sub set_cookies
{
    my ($self,@keys) = @_;

    for my $key (@keys)
    {
        my (@values) = $self->param($key);
        if ( @values ) 
        {
            $self->set_cookie($key,join('|',@values)); 
        }
        else
        {
            $self->unset_cookie($key);
        }
    }
}


sub get_cookies
{
    my ($self,@keys) = @_;

    for my $key (@keys)
    {
        my $value = $self->param($key);
        if ( defined $value ) { next; }

        $value = $self->get_cookie($key);
        if ( defined $value ) { $self->param($key,split(/\|/,$value)); }
    }
}


sub out
{
    my ($self,@args) = @_;
    if ( !$$self{'header_printed'} ) 
    {
        my @cookies = ();
        while (my ($key, $value) = each %{$$self{'set_cookies'}})
        {
            push @cookies, $$self{'cgi'}->cookie('-name'=>$key, '-value'=>$value, '-domain'=>'/', '-domain'=>'.sanger.ac.uk');
        }
        for my $key (@{$$self{'unset_cookies'}})
        {
            push @cookies, $$self{'cgi'}->cookie('-name'=>$key, '-value'=>'', '-expire'=>'-1h', '-domain'=>'/', '-domain'=>'.sanger.ac.uk');
        }
        if ( @cookies ) { $$self{'sw'}->cookie(\@cookies); }

        print $self->header(); 
        $$self{'header_printed'} = 1;
    }
    print @args;
    return;
}

sub header
{
    my ($self,@args) = @_;
    return $$self{'sw'}->header();
}

sub print_footer
{
    my ($self,@args) = @_;
    $self->out($$self{'sw'}->footer());
    return;
}

sub print_error
{
    my ($self,@msg) = @_;

    my $msg = scalar @msg ? join('',@msg) : $self->param('error');
    if ( !$msg ) { return; }

    $self->out(qq[
        <div class="error">
        <h3>A problem occurred</h3>
        <div style="padding-top:1em;">$msg</div>
        </div>
    ]);
}

sub error_exit
{
    my ($self,@msg) = @_;

    if ( scalar @msg )
    {
        my $msg = join('',@msg);

        $self->out(qq[
            <form action="$$self{myself}" method="POST" id="error_form">
            <input type="hidden" name="error" value="$msg" />
            </form>
            <script type="text/javascript">
            <!--
            jQuery("#error_form").submit();
            //-->
            </script>
        ]);

        # This code should be run only with javascript disabled.
        $self->print_error($msg);
    }
    $self->print_footer();
    exit;
}


sub help_link
{
    my ($self,$keyword) = @_;
    #return qq[<sup><a href="javascript:toggle_element_here('$keyword')" class="help_link">(?)</a></sup>];
    return qq[<span onclick="toggle_element_here('$keyword')" class="help_link"><sup>(?)</sup></span>];
}

sub print_help_divs
{
    my ($self) = @_;

    my $divs =
    {
        'help_loc' => q[
            The location can only be given e.g. as
            <ul class="example">
            <li>1:10,000,000-10,040,000 </li>
            <li>1: 10000000 - 10040000 </li>
            </ul>
            ],
    };

    for my $key (keys %$divs)
    {
        $self->out(qq[
            <div class="help" id="$key" onmouseout="hide_element('#$key')" onmouseover="show_element('#$key')">
                <div style="float:right;cursor:pointer;padding:0px;margin-top:-0.5em; margin-right:-0.5em;" onclick="hide_element('#$key')">[x]</div>
                $$divs{$key}
            </div>]);
    }
}


sub print_form
{
    my ($self) = @_;

    $self->get_cookies('strain','location','action','rows','conseqs','svtypes');
    $self->print_help_divs();

    # Prefill the form fields with values from the last query.
    #
    my $def_snps   = 'checked="checked"';
    my $def_indels = '';
    my $def_svs = '';
    if ( $self->param('action') && $self->param('action') eq 'indels' )
    {
        $def_snps   = '';
        $def_svs = '';
        $def_indels = 'checked="checked"';       
    }
    elsif ( $self->param('action') && $self->param('action') eq 'svs' )
    {
        $def_snps   = '';
        $def_svs = 'checked="checked"';
        $def_indels = '';            	
    }
    
    my $def_loc     = $self->param('location') ? $self->param('location') : '';
    my $def_rows      = $self->param('rows') ? $self->param('rows') : 50;
    my @conseqs = $self->param('conseqs') ?  split(/,/,$self->param('conseqs')) 
        : ('3PRIME_UTR','5PRIME_UTR','ESSENTIAL_SPLICE_SITE','NON_SYNONYMOUS_CODING','SYNONYMOUS_CODING',
            'NO_CODING_TRANSCRIPTS','STOP_GAINED','STOP_LOST','COMPLEX_INDEL','FRAME_SHIFT');
    my $selected_conseqs = {};
    for my $conseq (@conseqs)
    {
        $$selected_conseqs{$conseq} = 'checked="checked"';
    }
    for my $conseq (keys %{$$self{conseqs}})
    {
        if ( !$$selected_conseqs{$conseq} ) { $$selected_conseqs{$conseq}=''; }
    }
    
    my @sv_types = $self->param('svtypes') ? split(/,/,$self->param('svtypes')) 
        : ('Deletion', 'Insertion', 'Copy_Number_Gain', 'Inversion', 'Combined_events');
    my $selected_svtypes = {};
    for my $svtype (@sv_types)
    {
        $$selected_svtypes{$svtype} = 'checked="checked"';
    }
    for my $svtype (keys %{$$self{sv_types}})
    {
        if ( !$$selected_svtypes{$svtype} ) { $$selected_svtypes{$svtype}=''; }
    }

    my @strains = $self->param('strain') ? split(/,/,$self->param('strain')) : (keys %{$$self{mouseinfo}});
    my $selected_strains = {};
    for my $strain (@strains)
    {
        $$selected_strains{$strain} = 'checked="checked"';
    }
    for my $strain (keys %{$$self{mouseinfo}})
    {
        if ( !$$selected_strains{$strain} ) { $$selected_strains{$strain}=''; }
    }


    # Now print the form
    #
    $self->out(qq[
        <table style="width:100%;"><tr><td>
        <form action="$$self{myself}" method="post" id="form">
        <h3>Mouse Genomes Project - Query SNPs, indels or SVs</h3>
        <table><tr>
        <td>

        <fieldset>
        <legend>Located</legend>
        <table>
            <tr><td>.. in the region] .$self->help_link('#help_loc') . qq[<br />e.g. 2:10000000-10500000</td>
                <td> <input type="text" class="inputbox" name="location" size="25" value="$def_loc" /> </td>
                <td style="text-align:right;"><input type="submit" class="button_submit" style="padding:0.2em 0.5em 0.1em 0.5em;margin-left:1.5em;" value="Submit" name="get" /> </td></tr>
        </table>
        </fieldset>

        <fieldset>
        <legend>Options</legend>
        <table>
            <tr><td><input type="radio" name="action" $def_snps  value="snps"  /> <b>Show SNPs</b></td><td>(<a href="ftp://ftp-mouse.sanger.ac.uk/REL-1105/" target="_ftp">REL-1005</a>)</td></tr>
            <tr><td><input type="radio" name="action" $def_indels value="indels"  /> <b>Show indels</b></td><td>(<a href="ftp://ftp-mouse.sanger.ac.uk/REL-1105/" target="_ftp">REL-1005</a>)</td></tr>
            <tr><td><input type="radio" name="action" $def_svs  value="structvar"  /> <b>Show SVs</b></td><td>(<a href="ftp://ftp-mouse.sanger.ac.uk/REL-110804-SV/" target="_ftp">REL-110804-SV</a>)</td></tr>
            <tr><td colspan="2">Show max <input type="text" class="inputbox" size="3" name="rows" value="$def_rows" /> rows on page</td></tr>
        </table>
        </fieldset>


		
		<fieldset>
        <legend>Filters</legend>
		<table style="text-align:left;">
        	<tr><td colspan="2"><b>Consequence type (SNPs/Indels)</b></td></tr>
        	<tr><td colspan = "2">
                <div style="padding-bottom:0.4em;">
                (<span class="button_sel" onclick="check_all('cnsq_',true);">All</span> / <span class="button_sel" onclick="check_all('cnsq_',false);">None</span>)
                </div>
                </td>
            </tr>        
        ]);
        
    my $counter = 0;    
    for my $key (sort keys %{$$self{'conseqs'}})
    {
        my $conseq = $$self{'conseqs'}{$key};
        if ($counter % 2) {
        	$self->out(qq[<td><input type="checkbox" $$selected_conseqs{$key} name="cnsq_$key"> $$conseq{'label'}</td></tr>\n]);
        }
        else {
        	$self->out(qq[<tr><td><input type="checkbox" $$selected_conseqs{$key} name="cnsq_$key"> $$conseq{'label'}</td>\n]);
        }
        $counter++;	
    }
    if (!$counter % 2) {
        	$self->out(qq[</tr>\n]);
        }
    $self->out(qq[
                </table>
                <hr size="2"/>
 				<table style="text-align:left;">
            	<tr><td colspan="2"><b>Structural variation type (SVs only)</b></td></tr>
            	<tr><td colspan = "2">
            		<div style="padding-bottom:0.4em;">
                		(<span class="button_sel" onclick="check_all('svt_',true);">All</span> / <span class="button_sel" onclick="check_all('svt_',false);">None</span>)
                	</div>
                </td></tr>]);
    $counter = 0;    	
    for my $key (sort keys %{$$self{'sv_types'}})
    {
        my $svtype = $$self{'sv_types'}{$key};
        if ($counter % 2) {
        	$self->out(qq[<td><input type="checkbox" $$selected_svtypes{$key} name="svt_$key"> $$svtype{'front_label'}</td></tr>\n]);
        }
        else {
        	$self->out(qq[<tr><td><input type="checkbox" $$selected_svtypes{$key} name="svt_$key"> $$svtype{'front_label'}</td>\n]);
        }
        $counter++;        
    }
    if (!$counter % 2) {
        	$self->out(qq[</tr>\n]);
    }    
    $self->out(qq[
        </table>
        </fieldset>
		</div>
        </td>        
        <td width="150px">
            <fieldset>
            <legend>From the strains</legend>
            <div style="padding-bottom:0.4em;">
            (<span class="button_sel" onclick="check_all('str_',true);">All</span> / <span class="button_sel" onclick="check_all('str_',false);">None</span>)
            </div>
        ]);
    for my $strain (sort keys %{$$self{'mouseinfo'}}) 
    { 
        my $id = $strain;
        $id =~ s{/}{_}g;
		my $display_name = $$self{'mouseinfo'}{$strain}{'display_name'};        
        $self->out(qq[<input type="checkbox" $$selected_strains{$strain} name="str_$strain"> 
            <span onmouseover="show_element('#$id')" onmouseout="hide_element('#$id')">$display_name</span><br>\n]);
    }
    $self->out(qq[
            </fieldset>
        <fieldset>
            <legend>Further information</legend>
            Full datasets are available via <a href="ftp://ftp-mouse.sanger.ac.uk/" target="_ftp">FTP</a>.<br /><br />
            For further information please visit the project <a href="http://www.sanger.ac.uk/resources/mouse/genomes/" target="_blank">homepage</a>.
        ]);
    $self->out(qq[
            </fieldset>
        </td>
        
        </tr>
        <tr>
        <td colspan="3">
        <fieldset>
        <legend>References</legend>
		<table>
        <tr><td>Keane TM, Goodstadt L, Danecek P, <i>et al.</i> (2011) <br />Mouse genomic variation and its effect on phenotypes and gene regulation, <br /><a href = "http://dx.doi.org/10.1038/nature10413" target = "_natureref"><i>Nature</i><a>, 477(7364):289-294.</td></tr>
        <tr><td>Yalcin B, Wong K, Agam A, <i>et al.</i> (2011) <br />Sequence-based characterization of structural variation in the mouse genome, <br /><a href = "http://dx.doi.org/10.1038/nature10432" target = "_natureref"><i>Nature</i></a>, 477(7364):326-329.</td></tr> 
        </table>
        </fieldset>
        </td>
        </tr>
        </table>
        </form>

        </td><td>
        <div id="mouse_info" style="float:right;">
    ]);
    
    #OLD SUBMIT BUTTON:
#             <tr><td colspan="2" style="text-align:right;">
#             <input type="submit" class="button" style="padding:0.2em 0.5em 0.1em 0.5em;margin-right:0.5em;" value="Submit" name="get" /> 
#         </td></tr>
    
    for my $strain (@{$$self{'strains'}})
    {
        my $caption = $$self{'mouseinfo'}{$strain}{'caption'};
        my $img     = $$self{'mouseinfo'}{$strain}{'img'} ? qq[<img src="$$self{imgs}/$$self{'mouseinfo'}{$strain}{'img'}" alt="$strain" />] : '';
        my $id      = $strain;
        $id =~ s{/}{_}g;
        $self->out(qq[
            <div id="$id" style="display:none;" class="mouseinfo">
                <div style="float:right;">$img</div>
                <h3>$strain</h3>
                <div style="padding-top:1em;">$caption</div>
            </div>
            ]);
    }
    $self->out(qq[</div></td></tr>
    </table>]);
}

1;

