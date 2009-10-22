#!/usr/local/bin/perl -T
# qcgrind.pl
# 
# Displays QC information for Sanger short-read sequencing to allow lanes
# to be passed/failed
#
# Author:        jws
# Maintainer:    jws
# Created:       2009-03-10
# Updated:       2009-08-07 modified to run against new tracking API/db

use strict;
use warnings;
no warnings 'uninitialized';

use SangerPaths qw(core team145);
use VRTrack::VRTrack;
use VRTrack::Project;
use VRTrack::Sample;
use VRTrack::Library;
use VRTrack::Lane;
use SangerWeb;

$ENV{PATH}= '/usr/local/bin';   # solely to stop taint from barfing

$|++;


our $SEQ_LIMIT = 10_000;
our $SIG_COLOUR = "#CCFF00";
our $SEQ_COLOUR = "#FFFF00";
our $CUT_COLOUR = "#FF9900";
our $CONT_COLOUR = "#CC99FF";
my $javascript = <<JS ;
<!-- 
function toggleAdv() {
    var block = document.getElementById("advopt");
    var link = document.getElementById("advtoggleanchor");
    if (block.style.display == "none"){
        block.style.display = "block";
	link.firstChild.nodeValue="Hide advanced options";
    } else {
        block.style.display = "none";  
	link.firstChild.nodeValue="Show advanced options";
    }
}

// -->
JS

my $css = <<CSS ;

.centerFieldset {
text-align:center;
}

.centerFieldset fieldset {
margin-left:auto;
margin-right:auto;
/* INHERITED ALIGNMENT IS CENTER. ONLY INCLUDE THIS IF YOU WANT */
/* TO CHANGE THE ALIGNMENT OF THE CONTENTS OF THE FIELDSET */
text-align:left;
}

.centerFieldset table {
    margin-left:auto;
    margin-right:auto;
}

table.summary {
    border-collapse: collapse;
    font: 0.9em Verdana, Arial, Helvetica, sans-serif;
}

table.summary td {
    white-space: nowrap;
    text-align: right;
    padding-right: 1em;
    padding-left: 1em;
    padding-top: 2px;
    padding-bottom: 2px;
    border-bottom: #83a4c3 dotted 1px;
}

table.summary tr.header th, table.summary tr.header td {
    font-weight:bold;
    border-bottom: #83a4c3 solid 1px;
    text-align: left;
    vertical-align:middle;
}

table.summary tr.level th, table.summary tr.level td {
    font-weight:bold;
    border-bottom: #83a4c3 dotted 1px;
    text-align: left;
    vertical-align:middle;
}

table.summary tr.total th, table.summary tr.total td {
    font-weight:bold;
    background-color: #CBDCED;
    border-bottom: #83a4c3 dotted 1px;
    text-align: right;
    vertical-align:middle;
}



input.btn {
  font: bold 150% 'trebuchet ms',helvetica,sans-serif;
  border:1px solid;
  border-color: #707070 #000 #000 #707070;
}

input.btnhov { 
  cursor:pointer;
  border-color: #c63 #930 #930 #c63; 
}

.clear
{
    clear: both;
    display: block;
    overflow: hidden;
    visibility: hidden;
    width: 0;
    height: 0;
}



CSS


my $sw  = SangerWeb->new({
    'title'   => q(QC Grind),
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
    'script'  => $javascript,
    'style'   => $css,
    'stylesheet' => '/Teams/Team145/hoverbox.css',
});


# ENTRY POINTS INTO SCRIPT
# 1) No params : show list of species
# 2) Lane, no index : get list of lanes, find index, show lane
# 3) submit of pass/fail/pend : update lane, get list of lanes, redirect to
# next lane
# 4) Lane && index - show lane
# The point of index is simply to be able to find the next lane (or prev)
# without having to search the list of lanes for the current lane.

# 2009-03-15 jws - new entry points:
# species, project, sample, library - give summary of stats at the appropriate
# level and link lanes back to main lane view

# 2009-03-24 jws - new entry point:
# Also need to QC a library, independent of lanes.  This is because we might
# want to QC lanes as OK, but a library as failed if e.g. the insert size is too
# small.  Data might be good, but we don't want to order more on those libs.

# 2009-03-24 jws - add filters to display pend/pass/or fail

my $cgi = $sw->cgi();

my $spp		= $cgi->param('spp');
my $lane_name	= $cgi->param('lane');
my $update	= $cgi->param('update');
my $lib_update	= $cgi->param('lib_update');
my $lane_index	= $cgi->param('i');
my $sum_proj	= $cgi->param('p');
my $sum_samp	= $cgi->param('s');
my $sum_lib	= $cgi->param('l');
my @filter	= $cgi->param('filter');
my $add_filter  = $cgi->param('filt');


my %db_for_spp = ( 'mouse'  => 'mouse_reseq_track',
		    'g1k'   => 'g1k_track',
		  );

# ENTRY - no species, or wrong species
my $db = $db_for_spp{$spp};
unless ($db){
    print $sw->header();
    show_spp_list(\%db_for_spp);
    print $sw->footer();
    exit;
}
my $vrtrack = VRTrack::VRTrack->new({ host => 'mcs4a',
                                    port => 3306,
                                    user => 'vreseq_rw',
                                    password => 't3aml3ss',
                                    database => $db,
                                   });
my $dbh;
unless ($vrtrack){
    print $sw->header();
    print_and_exit( sprintf('DB connection failed: %s',$DBI::errstr));
}

# We have a species & a database connection: dispatch to right action
my $lanes = $vrtrack->qc_filtered_lane_names(@filter);
unless ($lane_name){
    $lane_name=$lanes->[0];
}
#if ($lane_name){
    if (scalar @$lanes == 0){
        print $sw->header();
        print_and_exit("There are no lanes with a QC status of ".(join " + ", @filter));
    }

    $lane_index = undef;    # skip indexes for now - this is just problematic
    $lane_index = undef if $lane_index < 0;
    if (! defined ($lane_index)){
	$lane_index = find_index_of_lane($lane_name, $lanes);
	if (! defined $lane_index){
	    if ($add_filter){	
		# pressed button to add new filter, so just use first lane in 
		# filtered set
		$lane_index = 0;
		$lane_name = $lanes->[0];
	    }
	    else {
		print $sw->header();
		print_and_exit("Can't find $lane_name in list of lanes");
	    }
	}
    }

    # Either update lane_name & show next, or just show lane_name
    if ($update){
	# no page header so we can redirect
	update_and_redirect($update, $lane_name, $lane_index, $lanes, \@filter);
    }    
    elsif($lib_update){
	lib_update_and_redirect($lib_update, $lane_name, $lane_index, $lanes, \@filter);
    }
    else {
	# MAIN PAGE DISPLAY
	print $sw->header();
	show_lane($lane_name, $lane_index, $lanes, \@filter);
    }
#}
#else {
#    if ($sum_proj){	# summarise project
#	print $sw->header();
#	show_summary('project', $sum_proj);
#    }
#
#    elsif ($sum_samp){	# summarise sample
#	print $sw->header();
#	show_summary('sample', $sum_samp);
#    }
#
#    elsif ($sum_lib){	# summarise library
#	print $sw->header();
#	show_summary('library', $sum_lib);
#    }
#    else {		# summarise species
#	print $sw->header();
#	show_summary('species');
#    }
#
#}

# if we made it here, we haven't exited or redirected, so close page & quit
print $sw->footer();
exit;

###############################################################################
###############################################################################
###############################################################################


# updates the qc_status of the library of the given lane_name, and redirects to the
# same lane_name
sub lib_update_and_redirect {
    my ($lib_update, $lane_name, $lane_index, $lanes, $filter) = @_;
    update_lib($lib_update, $lane_name);

    # redirect back to the same lane_name
    redirect_to_lane($lane_name, $lane_index, $filter);
}


# updates the qc_status of the given lane_name, and redirects to the next lane
sub update_and_redirect {
    my ($update, $lane_name, $lane_index, $lanes, $filter) = @_;
    update_lane($update, $lane_name);

    my ($next_lane, $next_index) = get_next_lane($lane_index, $lanes);
    redirect_to_lane($next_lane, $next_index, $filter);
}


sub filterurl {
    my $filter = shift;
    my $url;
    if (@$filter){
	$url = '&amp;filter=';
	$url .= join '&amp;filter=', @$filter;
    }
    return $url;

}

sub redirect_to_lane {
    my ($lane_name, $index, $filter) = @_;
    my $location = "qcgrind.pl?spp=$spp&amp;lane=$lane_name";

    # filters change the underlying set - force index reset
    if (@$filter){
	$location .= filterurl($filter);
    }
    else {
	$location .= "&amp;i=$index";
    }

    print $cgi->redirect( -location => $location,
			  -method   => 'GET',
			  -status   => 303);
    exit;

}

# performs the sql update to change a library qc status
sub update_lib {
    my ($update, $lane_name) = @_;
    my $lane = get_lane_by_name($lane_name);
    my ($project,$sample,$library) = get_hierarchy_for_lane($lane);
    my %status_for_update = (	'Pend' => 'pending',
				'Pass' => 'passed',
				'Fail' => 'failed'
			    );
    my $new_status = $status_for_update{$update};
    unless ($new_status){
	print $sw->header();
	print_and_exit("$update is not a valid QC status");
    }
    eval {
        $library->qc_status($new_status);
        $library->update;
    };
    if ($@){
	print $sw->header();
	print_and_exit(sprintf('Cannot update library to QC status $new_status: %s',$@));
    }
}

# performs the sql update to change a lane qc status
sub update_lane {
    my ($update, $lane_name) = @_;
    my %status_for_update = (	'Pend' => 'pending',
				'Pass' => 'passed',
				'Fail' => 'failed'
			    );
    my $new_status = $status_for_update{$update};
    unless ($new_status){
	print $sw->header();
	print_and_exit("$update is not a valid QC status");
    }
    my $lane = get_lane_by_name($lane_name);
    eval {
        $lane->qc_status($new_status);
        $lane->update;
    };
    if ($@){
	print $sw->header();
	print_and_exit(sprintf('Cannot update lane $lane_name to QC status $new_status: %s', $@));
    }
}

# get the index of the given lane in the list of all lanes
sub find_index_of_lane {
    my ($lane_name, $lanes) = @_;
    my $index;
    for (my $i=0; $i<scalar(@$lanes); ++$i){
	if ($lanes->[$i] eq $lane_name){
	    $index = $i;
	    last;
	}
    }
    return $index;
}

# Main display of script
# display QC information for one lane.
sub show_lane {
    my ($lane_name, $lane_index, $lanes, $filter) = @_;
    my $lane = get_lane_by_name($lane_name);
    my ($project,$sample,$library) = get_hierarchy_for_lane($lane);
    my $project_name = $project->hierarchy_name;
    my $sample_name = $sample->hierarchy_name;
    my $library_name = $library->hierarchy_name;
    $lane_name =~ /(\d+)_\d/;
    my $run = $1;
    my $libcount = scalar @{$library->lane_ids};
    print qq(<h2 style="font: italic normal 900 1.5em arial">QC Grind 
    <a href="qcgrind.pl?spp=$spp">$spp</a> - 
    <a href="qcgrind.pl?spp=$spp&amp;p=$project_name">$project_name</a> : 
    <a href="qcgrind.pl?spp=$spp&amp;s=$sample_name">$sample_name</a> : 
    <a href="qcgrind.pl?spp=$spp&amp;l=$library_name">$library_name</a> : 
    $lane_name [<a href="http://intweb.sanger.ac.uk/perl/prodsoft/npg/npg/run/$run">npg</a>]
    </h2>);

    my $status = $lane->qc_status;
    my $lib_status = $library->qc_status;
    my $mapstats = $lane->latest_mapping;
    my $npg_qc              = $lane->npg_qc_status;

    my $imglist = [];
    my ($reads_mapped, $bases_mapped, $reads_paired, $rmdup_reads_mapped, $rmdup_bases_mapped, $clip_bases,$adapter_reads);

    my $cycles		    = $lane->read_len;
    my $reads		    = $lane->raw_reads;
    my $bases		    = $lane->raw_bases;
    my $total_bases         = $bases;
    # derived stats
    my $dupe_rate	    = "";
    my $reads_mapped_perc   = "";
    my $reads_paired_perc   = "";
    my $bases_mapped_perc   = "";
    my $rmdup_reads_mapped_perc   = "";
    my $rmdup_bases_mapped_perc   = "";
    my $clip_bases_perc     = "";
    my $error_rate          = "";
    my $gt_status           = $lane->genotype_status;
    my $gt_display          = $gt_status;
    my $gt_found;
    my $gt_ratio;
    my $adapter_perc = '';

    if ($mapstats){
        $imglist = $mapstats->images;
        
        # over-ride raw_reads and raw_bases if we have a mapping
        # This handles sampling for mapping
        $reads              = $mapstats->raw_reads;
        $bases              = $mapstats->raw_bases;
        $clip_bases         = $mapstats->clip_bases;    # bases after clipping
        $total_bases        = defined $clip_bases ? $clip_bases : $bases;
        $reads_mapped	    = $mapstats->reads_mapped;
        $bases_mapped	    = $mapstats->bases_mapped;
        $reads_paired	    = $mapstats->reads_paired;
        $rmdup_reads_mapped  = $mapstats->rmdup_reads_mapped;
        $rmdup_bases_mapped  = $mapstats->rmdup_bases_mapped;
        $error_rate	    = sprintf("%.3f",$mapstats->error_rate);
        if (defined $mapstats->adapter_reads){
            $adapter_reads = $mapstats->adapter_reads;
            $adapter_perc = sprintf("%.1f",($adapter_reads/$reads)*100);

        }

        # genotypes
        $gt_found           = $mapstats->genotype_found;
        $gt_ratio           = sprintf("%.3f",$mapstats->genotype_ratio);
        $gt_display = "$gt_status ($gt_found:$gt_ratio)";


        if ($bases_mapped){	# sometimes the mapping fails
            $dupe_rate	    = sprintf("%.4f", (1-$rmdup_reads_mapped/$reads_mapped)); #($bases_mapped - $rmdup_bases_mapped)/$bases_mapped);
            $reads_mapped_perc   = sprintf("%.1f", ($reads_mapped/$reads)*100);
            $reads_paired_perc   = sprintf("%.1f", ($reads_paired/$reads)*100);
            $bases_mapped_perc   = sprintf("%.1f", ($bases_mapped/$total_bases)*100);
            $rmdup_reads_mapped_perc   = sprintf("%.1f", ($rmdup_reads_mapped/$reads)*100);
            $rmdup_bases_mapped_perc   = sprintf("%.1f", ($rmdup_bases_mapped/$total_bases)*100);
            $clip_bases_perc    = sprintf("%.1f", ($clip_bases/$bases)*100);

            $reads_mapped = commify($reads_mapped);
            $bases_mapped = commify($bases_mapped);
            $reads_paired = commify($reads_paired);
            $rmdup_reads_mapped = commify($rmdup_reads_mapped);
            $rmdup_bases_mapped = commify($rmdup_bases_mapped);
        }

        $reads = commify($reads);
        $bases = commify($bases);
        $clip_bases = commify($clip_bases) if $clip_bases;

        # visual checks
        unless ($gt_status eq 'confirmed'){
            $gt_display = qq[<span style="color : red;">$gt_display</span>];
        }
        if ($bases_mapped_perc < 80){
            $bases_mapped_perc = qq[<span style="color : red;">$bases_mapped_perc</span>];
        }
        if ($mapstats->error_rate > 0.02){
            $error_rate = qq[<span style="color : red;">$error_rate</span>];
        }
        if ($npg_qc eq 'fail'){
            $npg_qc = qq[<span style="color : red;">$npg_qc</span>];
        }

    }

    # QC plots
     print qq[
    <br />
    <div class="centerFieldset">
    <fieldset style="width: 800px" > 
    <legend>QC plots</legend>
    ];

    my %images = map {$_->name => $_} @$imglist;

    # print gc dist then insert size, then any others smaller
    my $gcimg = $images{'gc.gif'} || $images{'gc-content.png'};
    print qq[<table><tr><td>];
    if ($gcimg){
	my $id = $gcimg->id;
	my $caption = $gcimg->caption;
        print qq[
        <img src="vrimg.pl?spp=$spp&img=$id" alt="$caption" width="300" height="300" style="border:1px dotted #83A4C3;"/></a>
        ];
        delete $images{$gcimg->name};
    }
    my $insimg = $images{'rmdup.inserts.gif'} || $images{'insert-size.png'};
    if ($insimg){
	my $id = $insimg->id;
	my $caption = $insimg->caption;
        print qq[
        <img src="vrimg.pl?spp=$spp&img=$id" alt="$caption" width="300" height="300" style="border:1px dotted #83A4C3;"/></a>
        ];
        delete $images{$insimg->name};
    }
    print qq[</td></tr></table>];

    # now the rest
    #print qq[<span class="clear"></span>
    print "<table><tr><td>";
    print qq[<ul class="hoverbox">];
    foreach (sort keys %images){
        my $img = $images{$_};
	my $id = $img->id;
	my $caption = $img->caption;
        print qq[<li>
        <a href="#"><img src="vrimg.pl?spp=$spp&img=$id" alt="$caption" /><img src="vrimg.pl?spp=$spp&img=$id" alt="$caption" class="preview" /></a>
        </li>
        ];
    }
    print qq[</ul>];
    print qq[</td></tr></table>];

    print qq[
    </fieldset>
    </div>
    ];

    # Mapping stats
    my ($adapter_str, $rmapped_str, $rpaired_str, $rmdup_rmapped_str);
    my ($clip_str, $bmapped_str, $rmdup_bmapped_str);
    
    $adapter_str = "$adapter_reads ($adapter_perc\%)" if $adapter_reads;
    $rmapped_str = "$reads_mapped ($reads_mapped_perc\%)" if $reads_mapped;
    $rpaired_str = "$reads_paired ($reads_paired_perc\%)" if $reads_paired;
    $rmdup_rmapped_str = "$rmdup_reads_mapped ($rmdup_reads_mapped_perc\%)" if $rmdup_reads_mapped;

    if (defined $clip_bases){
        $clip_str = "$clip_bases ($clip_bases_perc\% raw)";
    }
    else {  
        $clip_str = "unclipped";
    }
    $bmapped_str = "$bases_mapped ($bases_mapped_perc\%)" if $bases_mapped;
    $rmdup_bmapped_str = "$rmdup_bases_mapped ($rmdup_bases_mapped_perc\%)" if $rmdup_bases_mapped;
    print qq[
    <br />
    <div class="centerFieldset">
    <fieldset style="width: 900px">
    <legend>Mapping data</legend>
    <table width="100%">
    <tr>
    <td><table>
	<tr><td>Reads: </td><td>$reads</td></tr>
	<tr><td>Reads w/adapter</td><td>$adapter_str</td></tr>
	<tr><td>Reads mapped: </td><td>$rmapped_str</td></tr>
	<tr><td>Reads paired: </td><td>$rpaired_str</td></tr>
	<tr><td>Reads mapped (rmdup): </td><td>$rmdup_rmapped_str</td></tr>
	</table>
    </td>
    <td><table>
	<tr><td>Bases: </td><td>$bases</td></tr>
	<tr><td>Bases postclip: </td><td>$clip_str</td></tr>
	<tr><td>Bases mapped: </td><td>$bmapped_str</td></tr>
	<tr><td>&nbsp;</td><td>&nbsp;</td></tr>
	<tr><td>Bases mapped (rmdup): </td><td>$rmdup_bmapped_str</td></tr>
	</table>
    </td>
    <td><table>
	<tr><td>Cycles: </td><td>$cycles</td></tr>
	<tr><td>NPG QC: </td><td>$npg_qc</td></tr>
	<tr><td>Error rate: </td><td>$error_rate</td></tr>
	<tr><td>Genotype check: </td><td>$gt_display</td></tr>
	<tr><td>Duplication rate: </td><td>$dupe_rate</td></tr>
    </table>
    </td>
    </tr>
    </table>
    </fieldset>
    </div>
    ];


    # Lane status block:
    my $status_colour = get_colour_for_status($status);
    my $lanecount = scalar @$lanes;
    my $lanenum = $lane_index + 1;
    my ($next_lane, $next_index) = get_next_lane($lane_index, $lanes);
    my ($prev_lane, $prev_index) = get_prev_lane($lane_index, $lanes);
    my $statehtml = get_state_html($lane_name, $spp, $lane_index,$filter);
    my $filterstr = filterurl($filter);
    print qq[<br />
    <div class="centerFieldset">
    <fieldset style="background-color: $status_colour;width: 500px">
    <legend>Lane QC : $status</legend>
    <form action="qcgrind.pl">
    <p>
    $statehtml
    <table width="400px" align="center">
    <tr>
    <td align="center"><input type="submit" name="update" value="Pend" class="btn" style="background-color: #F5F5F5;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
    <td align="center"><input type="submit" name="update" value="Fail" class="btn" style="background-color: #FFC0C0;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
    <td align="center"><input type="submit" name="update" value="Pass" class="btn" style="background-color: #C0FFC0;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
    </tr>
    <tr>
	<td align="left"><a href="qcgrind.pl?spp=$spp&amp;lane=$prev_lane&amp;i=${prev_index}${filterstr}">&lt;</a></td>
	<td align="center">$lanenum / $lanecount</td>
	<td align="right"><a href="qcgrind.pl?spp=$spp&amp;lane=$next_lane&amp;i=${next_index}${filterstr}">&gt;</a></td>
    </tr>
    </table>
    </p>
    </form>
    </fieldset>
    </div>
    ];

    # Library status block
    my $lib_status_colour = get_colour_for_status($lib_status);
    print qq[<br />
    <div class="centerFieldset">
    <fieldset style="background-color: $lib_status_colour;width: 500px">
    <legend>Library QC : $lib_status</legend>
    <form action="qcgrind.pl">
    <p>
    $statehtml
    <input type="hidden" name="i" value="$lane_index" />
    <table width="400" align="center">
    <tr>
    <td align="center"><input type="submit" name="lib_update" value="Pend" class="btn" style="background-color: #F5F5F5;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
    <td align="center"><input type="submit" name="lib_update" value="Fail" class="btn" style="background-color: #FFC0C0;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
    <td align="center"><input type="submit" name="lib_update" value="Pass" class="btn" style="background-color: #C0FFC0;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
    </tr>
    <tr>
	<td align="left"></td>
	<td align="center">lanes in library: $libcount</td>
	<td align="right"></td>
    </tr>
    </table>
    </p>
    </form>
    </fieldset>
    </div>
    ];

    print qq[<br />
    <div class="centerFieldset">
    <fieldset style="width: 500px;text-align:center;">
    <legend>Filter: </legend>
    <form action="qcgrind.pl">
    <p>
    ];
    foreach my $filt( 'no_qc', 'pending', 'passed', 'failed'){
	my $checked = 'checked="yes"' if grep {$_ eq $filt} @$filter;
	print qq[$filt: <input type="checkbox" name="filter" value="$filt" $checked />&nbsp;&nbsp;];
    }
    print qq[
    <input type="submit" name="filt" value= "Filter" />
    <input type="hidden" name="lane" value="$lane_name" />
    <input type="hidden" name="spp" value="$spp" />
    </p>
    </form>
    </fieldset>
    </div>
    ];
}


sub get_next_lane {
    my ($lane_index, $lanes) = @_;

    my ($next_lane, $next_index);
    if ($lane_index == $#{$lanes}){
	$next_lane = $lanes->[$lane_index];
	$next_index = $lane_index;
    }
    else {
	$next_index = $lane_index + 1;
	$next_lane = $lanes->[$next_index];
    }

    return ($next_lane, $next_index);
}


sub get_prev_lane {
    my ($lane_index, $lanes) = @_;

    my ($prev_lane, $prev_index);
    if ($lane_index == 0){
	$prev_lane = $lanes->[$lane_index];
	$prev_index = $lane_index;
    }
    else {
	$prev_index = $lane_index - 1;
	$prev_lane = $lanes->[$prev_index];
    }

    return ($prev_lane, $prev_index);

}


# returns CSS colour for a QC status
sub get_colour_for_status {
    my $status = shift;
    my $status_colour;
    
    if ($status eq 'no_qc'){
	$status_colour="#FFFFFF";
    }
    if ($status eq 'passed'){
	$status_colour="#C0FFC0";
    }
    elsif ($status eq 'failed'){
	$status_colour="#FFC0C0";
    }
    else {
	$status_colour="#F5F5F5";
    }

    return $status_colour;
}

# return number of lanes in this library
sub get_lanecount_for_library {
    my ($project, $sample, $library) = @_;
    my $sql = qq[select count(lane.name) from lane, library, sample, project where lane.library_id = library.library_id and library.sample_id = sample.sample_id and sample.project_id = project.project_id and project.name=? and sample.name= ? and library.name=?];
    my $sth = $dbh->prepare($sql);
    my $count;
    if ($sth->execute( $project, $sample, $library )) {
	my $data = $sth->fetchall_arrayref()->[0];
	($count) = @$data;
    }
    else {
	print_and_exit(sprintf('Cannot retrieve lane count for library $library: %s', $DBI::errstr));
    }
    return $count;
}


# returns the project, sample, library objects for a given lane
sub get_hierarchy_for_lane {
    my $lane = shift;
    my ($proj, $samp, $lib);
    eval {
        $lib = VRTrack::Library->new($vrtrack,$lane->library_id);
        if ($lib){
            $samp = VRTrack::Sample->new($vrtrack,$lib->sample_id);
            if ($samp){
                $proj = VRTrack::Project->new($vrtrack,$samp->project_id);
            }
        }
    };

    unless ($proj && $samp && $lib){
	print_and_exit(sprintf('Cannot retrieve hierarchy for lane %s', $lane->hierarchy_name));
    }
    return ($proj, $samp, $lib);
}

# returns the project, sample, library objects for a given lane
sub get_hierarchy_for_lib 
{
    my $lib = shift;
    my ($proj, $samp);
    eval {
            $samp = VRTrack::Sample->new($vrtrack,$lib->sample_id);
            if ($samp){
                $proj = VRTrack::Project->new($vrtrack,$samp->project_id);
            }
        }
    };

    unless ($proj && $samp)
	{
		print_and_exit(sprintf('Cannot retrieve hierarchy for lane %s', $ib->hierarchy_name));
    }
    return ($proj, $samp);
}


# list species for selection
sub show_spp_list {
    my $spphash = shift;
    print qq(<h2 style="font: italic normal 900 1.5em arial">QC Grind</h2>);
    print qq[
	<br />
	<div class="centerFieldset">
	<fieldset style="width: 900px">
	<legend>Select dataset to QC</legend>
	];

    foreach my $spp (keys %$spphash){
	my $sppname = ucfirst($spp);
	print qq[<p><a href="qcgrind.pl?spp=$spp">$sppname</a></p>\n];
    }

    print qq[
    </fieldset>
    </div>
    ];

}


sub show_summary {
    my ($level, $datum) = @_;

}


# summarise the data at the level specified and list lanes
sub show_summary_old {
    my ($level, $datum) = @_;
    my $show_all;
    if ($level eq 'species'){
	$show_all = 1;
	$level = 'project';
    }
    my @levels = ('project', 'sample', 'library','lane');
    my( $index )= grep { $levels[$_] eq $level } 0..$#levels;
    unless(defined $index){
	print_and_exit("$level is not a hierarchy level");
    }
    
    my $join_sql;
    if ($index < $#levels){ # do we need to join tables?
	for(my $i=$index;$i < $#levels;++$i){
	    my $from = $levels[$i];
	    my $to = $levels[$i+1];
	    $join_sql .= " and $from.${from}_id = $to.${from}_id";
	}
    }
    my $select_sql = join ",", map{"$_.name as ${_}_name"} @levels[$index..$#levels]; 
    my $from_sql = join ",", @levels[$index..$#levels];
    my $order_sql = join ",", map{"$_.name"} @levels[$index..$#levels];

    my @metafields = qw(raw_reads raw_bases reads_mapped reads_paired bases_mapped rmdup_reads_mapped rmdup_bases_mapped readlen error_rate); # all lane fields, in retrieval order
    my $metafields = join ",",@metafields;
    my @sumfields = qw(raw_reads raw_bases reads_mapped reads_paired bases_mapped rmdup_reads_mapped rmdup_bases_mapped); # fields to summarise, in retrieval order

    my $sql;
    if ($show_all){
	$sql = qq[select $select_sql, $metafields from $from_sql where 1=1 $join_sql order by $order_sql];
    }
    else {
	$sql = qq[select $select_sql, $metafields from $from_sql where $level.name = ? $join_sql order by $order_sql];
    }
    #warn $sql;
    my $sth = $dbh->prepare($sql);
    my $data;
    if ($show_all){
	if ($sth->execute) {
	    $data = $sth->fetchall_arrayref();
	}
	else {
	    print_and_exit(sprintf('Cannot retrieve summary data for %s: %s', $level, $DBI::errstr));
	}
    }
    else {
	if ($sth->execute($datum)) {
	    $data = $sth->fetchall_arrayref();
	}
	else {
	    print_and_exit(sprintf('Cannot retrieve summary data for %s %s: %s', $level, $datum, $DBI::errstr));
	}
    }


    # pass through data to output each lane and summarise at each level
    my %summary_for;	# running totals for each level
    my @sum_levels = @levels[$index..$#levels-1]; # actual levels we're summarising
    my $num_levels = scalar @sum_levels;
    my $num_cols = $num_levels+1 + scalar @metafields;
    print qq[<table  class="summary">\n];
    print_table_row([map{my $foo=$_; $foo=~s/_/ /g;ucfirst($foo)} @sum_levels,'lane',@metafields], $num_cols, "header");

    foreach my $row (@$data){
	my $cells = scalar @$row;
	# columns in the row are:
	# first $num_levels columns are level names
	# next field is the lane name
	# rest of fields are the @metafields

	# check if we have changed level,
	# update level running totals,
	# output lane data

	# Key the current level name by concatenating it with its parent level
	# names, so that if a parent level (e.g. project) changes, then by
	# definition, so does the subordinate level.  This gets around issues
	# where the project changes but the sample name is the same (all
	# individuals called "1" in Mouse)

	# need to print out the totals going up the hierarchy, but the level
	# headers going down
	# i.e. Project, then sample, but summarise sample first, then project.
	# because the output is nested.
	# So, do two loops - first to output summaries (going up)
	# second to output new headers (going down)

	# Output summaries of changed levels
	my @new_levels;

	# handle output by loading everything into an array of row arrays.
	# Each row will have data in the appropriate cell for output
	my @table_rows;
	for(my $i = $num_levels - 1;$i >= 0;--$i){

	    my $this_level = $sum_levels[$i];
	    my $this_level_key = join "\t", @$row[0..$i];
	    my $this_level_name = $row->[$i];
	    if ($summary_for{$this_level}{key} ne $this_level_key){
		# level has changed
		# print summary of previous level
		if (keys %{$summary_for{$this_level}{totals}}){
		    print_summary_row($summary_for{$this_level}{totals},\@sumfields, $i, $num_levels + 1, $num_cols);
		}
		# record changed level for output later
		$new_levels[$i] = $this_level_name;

		# reset summary hash
		$summary_for{$this_level}{name} = $this_level_name;
		$summary_for{$this_level}{key} = $this_level_key;
		$summary_for{$this_level}{totals} = {};

	    }

	    # increment summaries at this level for the new row data
	    for (my $n = 0;$n < scalar @sumfields;++$n){
		my $col_idx = $n + $num_levels + 1; # offset into row columns
		$summary_for{$this_level}{totals}{$sumfields[$n]} += $row->[$col_idx];
	    }
	    $summary_for{$this_level}{totals}{count}++;

	}
	# print new level headers (going down hierarchy)
	for(my $i = 0;$i < scalar @new_levels; ++$i){
	    next unless defined $new_levels[$i];
	    my @table_row;
	    $table_row[$i] = $new_levels[$i];
	    print_table_row( \@table_row, $num_cols, 'span_level');
	}

	# Output lane data
	my @table_row;

	# lane name output
	my $lane_name = $row->[$num_levels];
	$table_row[$num_levels] = qq[<a href="qcgrind.pl?spp=$spp&amp;lane=$lane_name">$lane_name</a>];

	my $col = $num_levels;
	my $raw_reads = 0;
	my $raw_bases = 0;
	# get the base & read counts for percentages
	for (my $n = 0;$n < scalar @metafields;++$n){
	    $col++;
	    my $sf = $metafields[$n];	
	    my $col_idx = $n + $num_levels + 1; # offset into row columns
	    $raw_reads = $row->[$col_idx] if $sf eq 'raw_reads';
	    $raw_bases = $row->[$col_idx] if $sf eq 'raw_bases';
	}
	$col = $num_levels;
	for (my $n = 0;$n < scalar @metafields;++$n){
	    $col++;
	    my $sf = $metafields[$n];	
	    my $col_idx = $n + $num_levels + 1; # offset into row columns
	    $table_row[$col] = format_mapval($row->[$col_idx], $sf, $raw_reads, $raw_bases);
	}
	print_table_row(\@table_row, $num_cols);
    }
    
    # last - print summaries for stuff still in the summary buffer:
    for(my $i = $num_levels - 1;$i >= 0;--$i){
	my $this_level = $sum_levels[$i];
	if (keys %{$summary_for{$this_level}{totals}}){
	    print_summary_row($summary_for{$this_level}{totals},\@sumfields, $i, $num_levels + 1, $num_cols);
	}
    }
    print "</table>\n";

}

# print table rows for sparse arrays
sub print_table_row {
    my ($row, $num_cols, $style) = @_;
    my $span;
    if ($style eq "span_level"){    
	$style = 'level';
	$span = 1;
    }

    if ($style){
	print qq[<tr class="$style">\n]; 
    }
    else {
	print "<tr>\n";
    }
    for (my $i=0; $i<$num_cols; ++$i){
	if ($span && $row->[$i]){
	    my $span = $num_cols - $i;
	    print qq[<td colspan="$span">],$row->[$i],"</td>";
	    last;
	}
	else {
	    print "<td>",$row->[$i],"</td>";
	}
    }
    print "</tr>\n";
}

sub print_and_exit {
    my ($msg) = @_;
    print qq[<h2>A problem occurred</h2>\n];
    print qq[<p class="error1">$msg</p>\n];
    print $sw->footer();
    exit;
}

# hidden fields for maintaining state
sub get_state_html {
    my ($lane_name, $spp, $lane_index, $filter) = @_;
    my $html = qq[
    <input type="hidden" name="lane" value="$lane_name" />
    <input type="hidden" name="spp" value="$spp" />
    ];
    if (! @$filter){	# filters might change index set, so force recalc.
	$html .= qq[<input type="hidden" name="i" value="$lane_index" />];
    }

    foreach (@$filter){
	$html .= qq[<input type="hidden" name="filter" value="$_" />];
    }
    return $html;
}

# output summary totals row

sub print_summary_row {
    my ($sumhash, $fields, $label_ind, $index, $num_cols) = @_;
    # output totals
    my @table_row;
    $table_row[$label_ind] = "total";
    my $raw_reads = $sumhash->{raw_reads};
    my $raw_bases = $sumhash->{raw_bases};
    my $tot_idx = $index;
    foreach my $sf (@$fields){
	$table_row[$tot_idx++] = format_mapval($sumhash->{$sf}, $sf,$raw_reads, $raw_bases);
    }
    print_table_row(\@table_row, $num_cols,'total');

    # output means
    if ($sumhash->{count} > 1){
	@table_row=();
	$table_row[$label_ind] = "mean";
	my $mean_idx = $index;
	foreach my $sf (@$fields){
	    $table_row[$mean_idx++] = val_to_nearest_SI($sumhash->{$sf}/$sumhash->{count});
	}
	print_table_row(\@table_row, $num_cols,'total');
    }

}


# format mapping values for output
sub format_mapval { 
    my ($val, $name, $reads, $bases) = @_;
    my $valstr = val_to_nearest_SI($val);
    $name =~ s/^rmdup_//;
    if ($name =~ /^reads/){
	my $perc = sprintf "%.1f%%", ($val/$reads)*100 ;
	$valstr .= " ($perc)";
    }
    elsif ($name =~ /^bases/){
	my $perc = sprintf "%.1f%%", ($val/$bases)*100 ;
	$valstr .= " ($perc)";

    }
    return $valstr;
}

# add commas in to big numbers
sub commify { 
    my $number = shift;
    $number =  reverse $number;
    $number =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g; 
    return scalar reverse $number; 
}


# round numbers to nearest SI multiple
sub val_to_nearest_SI {
    my $val = shift;
    my @units = ('', 'K','M','G','T','P' );
     
    my $thou = int((length(sprintf("%.0f%s",$val)) - 1)/3);
    my $unit_str = $val;
    if ($thou){
	my $unit = $units[$thou];
	
	my $value = $val/( 10 ** ($thou * 3)) ;
	  
	$unit_str = sprintf( "%.2f%s", $value, " $unit" );
    }
    return $unit_str;

}

sub get_lane_by_name {
    my $lane_name = shift;
    my $lane;
    eval {
        $lane = VRTrack::Lane->new_by_hierarchy_name($vrtrack,$lane_name);
    };
    unless ($lane){
        print_and_exit(sprintf("Can't retrieve lane $lane_name from database<br>%s",$@));
    }
    return $lane;
}
