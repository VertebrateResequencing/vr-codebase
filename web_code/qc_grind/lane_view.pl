#!/usr/local/bin/perl -T
# Displays Lane QC information 

BEGIN {
    $ENV{VRTRACK_HOST} = 'mcs4a';
    $ENV{VRTRACK_PORT} = 3306;
    $ENV{VRTRACK_RO_USER} = 'vreseq_ro';
    $ENV{VRTRACK_RW_USER} = 'vreseq_rw';
    $ENV{VRTRACK_PASSWORD} = 't3aml3ss';
};

use strict;
use warnings;
use URI;

use SangerPaths qw(core team145);
use SangerWeb;
use VRTrack::Project;
use VertRes::Utils::VRTrackFactory;
use VertRes::QCGrind::Util;

my $QC_GRIND_SCRIPT = "./qc_grind.pl";
my $STUDY_LANES_SCRIPT = "./study_lanes.pl";
my $pending_view = "../pending-view/pending_view.pl";

my $sw  = SangerWeb->new({
    'title'   => q(QC Grind Lane View),
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
});
my $utl = VertRes::QCGrind::Util->new();

my $USER = $sw->username();
my $cgi = $sw->cgi();
my $SCRIPT_NAME = $cgi->url(-relative=>1);

print $sw->header();
#use Data::Dumper;print "<pre>parms:", Dumper($cgi->param()), "</pre>"; ## DEBUG

my $db = $cgi->param('db');
$utl->displayError( "No database ID",$sw ) unless $db;
my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $db, mode => 'rw');

my $laneID = $cgi->param('lane_id');
$utl->displayError( "No Lane ID",$sw ) unless $laneID;

displayLane( $cgi, $vrtrack, $db, $laneID);
print $sw->footer();
exit;

#####

sub displayLane
{
	my ( $cgi, $vrtrack, $database, $laneID) = @_;

	my $lane;
	if( $laneID == -1 )
	{
		#pick a random library
		my $lanes = $vrtrack->qc_filtered_lane_names();
		$lane = VRTrack::Lane->new_by_hierarchy_name($vrtrack,$lanes->[0]);

		displayError("Null lane returned: $lane->[0]" ) unless defined( $lane );
	}
	else
	{
		$lane = VRTrack::Lane->new( $vrtrack, $laneID );
	}

	unless ($lane){
		displayError("Can't retrieve lane $laneID\n");
	}

	my $name = $lane->name;
	my ($project,$sample,$library) = get_hierarchy_for_lane($vrtrack, $lane);

	my $pname = $project->name;
	my $pid = $project->id();
	my $sname = $sample->name;
	my $sid = $sample->id();
	my $libname = $library->name;
	my $lid = $library->id();
	$name =~ /(\d+)_(\d+)/;
	my $run = $1;

	my $SPECIES_VIEW = $utl->{MODES}{SPECIES_VIEW};
	my $PROJ_VIEW = $utl->{MODES}{PROJ_VIEW};
	my $LIB_VIEW = $utl->{MODES}{LIB_VIEW};
	print qq[
		<h2 align="center" style="font: normal 900 1.5em arial">QC Grind</h2>
		<h3 style="font: normal 700 1.5em arial">
		<a href="$SCRIPT_NAME?mode=$SPECIES_VIEW&amp;db=$database">].ucfirst($database).qq[</a> : 
		<a href="$SCRIPT_NAME?mode=$PROJ_VIEW&amp;proj_id=$pid&amp;db=$database">$pname</a> :
		<a href="$SCRIPT_NAME?mode=$PROJ_VIEW&amp;proj_id=$pid&amp;db=$database">$sname</a> :
		<a href="$SCRIPT_NAME?mode=$LIB_VIEW&amp;lib_id=$lid&amp;db=$database">$libname</a> :
		$name 
		];

	if( $run !~ /SRR|ERR/ )
	{
		print qq{[<a href="http://intweb.sanger.ac.uk/perl/prodsoft/npg/npg/run/$run">npg</a>]};
	}

	print "</h3>";

	my $status = $lane->qc_status;
	my $auto_qc_status = $lane->auto_qc_status();

#work out which mapstats has the QC data (i.e. check for images in the mapstats)
	my @mappings = @{ $lane->mappings() };
	my $mapstats;
	foreach( sort {$a->row_id() <=> $b->row_id()} @mappings )
	{
		my $map = $_;
		my $im = $map->images();
		if( @{$im} > 0 ){$mapstats = $map;}
	}

	my $npg_qc = $lane->npg_qc_status;
	my ($error_rate, $adapter_perc, $reads_mapped, $bases_mapped, $reads_paired, $rmdup_reads_mapped, $rmdup_bases_mapped, $clip_bases,$adapter_reads,
			$bait_near_bases_mapped, $target_near_bases_mapped, $bait_bases_mapped, $mean_bait_coverage, $bait_coverage_sd, $off_bait_bases,
			$reads_on_bait, $reads_on_bait_near, $reads_on_target, $reads_on_target_near, $target_bases_mapped, $mean_target_coverage,
			$target_coverage_sd, $targets, $target_bases_1X, $target_bases_2X, $target_bases_5X, $target_bases_10X, $target_bases_20X, $target_bases_50X, $target_bases_100X);
	my $cycles = $lane->read_len;
	my $raw_reads = $lane->raw_reads;
	my $raw_bases = $lane->raw_bases;
	my $reads = $raw_reads; # if no mapping
		my $bases = $raw_bases; # if no mapping
		my $gt_status = $lane->genotype_status;
	my $total_bases = $bases;

   	my $status_colour = $utl->get_colour_for_status($lane->qc_status);
	my $printed = 0;
	if ($mapstats)
	{
		my $images = $mapstats->images;
		my $exome_design = $mapstats->exome_design();
		my ($target_bases, $bait_bases, $bait_near_bases_mapped_perc, $target_near_bases_mapped_perc, $bait_bases_mapped_perc, $off_bait_bases_perc, $reads_on_bait_perc, $reads_on_bait_near_perc, $reads_on_target_perc, $reads_on_target_near_perc, $target_bases_mapped_perc);

#my $bait_bases = $exome_design->bait_bases();
#my $target_bases = $exome_design->target_bases();

# over-ride raw_reads and raw_bases if we have a mapping
# This handles sampling for mapping
		$reads = $mapstats->raw_reads;
		$bases = $mapstats->raw_bases;
		$clip_bases = $mapstats->clip_bases;    # bases after clipping
			$total_bases = defined $clip_bases ? $clip_bases : $bases;
		$reads_mapped = $mapstats->reads_mapped;
		$bases_mapped = $mapstats->bases_mapped;
		$reads_paired = $mapstats->reads_paired;
		$rmdup_reads_mapped = $mapstats->rmdup_reads_mapped;
		$rmdup_bases_mapped = $mapstats->rmdup_bases_mapped;
		$error_rate = sprintf("%.3f",$mapstats->error_rate);
		if (defined $mapstats->adapter_reads)
		{
			$adapter_reads = $mapstats->adapter_reads;
			$adapter_perc = sprintf("%.1f",($adapter_reads/$reads)*100);
		}

#get the aligner
		my $mapper = $mapstats->mapper();
		my $mapper_string = qq[None found];
		if( $mapper )
		{
			$mapper_string = $mapper->name().qq[ ].$mapper->version();
		}

# genotypes
		my $gt_found = $mapstats->genotype_found;
		my $gt_ratio = sprintf("%.3f",$mapstats->genotype_ratio);
		my $gt_display = "$gt_status ($gt_found:$gt_ratio)";

# exome stats
		if (defined $exome_design) {
			$target_bases = $utl->commify($exome_design->target_bases());
			$bait_bases = $utl->commify($exome_design->bait_bases());
			$bait_near_bases_mapped = $mapstats->bait_near_bases_mapped;
			$target_near_bases_mapped = $mapstats->target_near_bases_mapped;
			$bait_bases_mapped = $mapstats->bait_bases_mapped;
			$mean_bait_coverage = $mapstats->mean_bait_coverage;
			$bait_coverage_sd = $mapstats->bait_coverage_sd;
			$off_bait_bases = $mapstats->off_bait_bases;
			$reads_on_bait = $mapstats->reads_on_bait;
			$reads_on_bait_near = $mapstats->reads_on_bait_near;
			$reads_on_target = $mapstats->reads_on_target;
			$reads_on_target_near = $mapstats->reads_on_target_near;
			$target_bases_mapped = $mapstats->target_bases_mapped;
			$mean_target_coverage = $mapstats->mean_target_coverage;
			$target_coverage_sd = $mapstats->target_coverage_sd;
			$target_bases_1X   = sprintf("%.1f", 100 * $mapstats->target_bases_1X);
			$target_bases_2X   = sprintf("%.1f", 100 * $mapstats->target_bases_2X);
			$target_bases_5X   = sprintf("%.1f", 100 * $mapstats->target_bases_5X);
			$target_bases_10X  = sprintf("%.1f", 100 * $mapstats->target_bases_10X);
			$target_bases_20X  = sprintf("%.1f", 100 * $mapstats->target_bases_20X);
			$target_bases_50X  = sprintf("%.1f", 100 * $mapstats->target_bases_50X);
			$target_bases_100X = sprintf("%.1f", 100 * $mapstats->target_bases_100X);
		}

		if ($bases_mapped)
		{   # sometimes the mapping fails
			my $dupe_rate_perc       = sprintf("%.2f", (1-($rmdup_reads_mapped/$reads_mapped))*100); #($bases_mapped - $rmdup_bases_mapped)/$bases_mapped);
			my $reads_mapped_perc   = sprintf("%.1f", ($reads_mapped/$reads)*100);
			my $reads_paired_perc   = sprintf("%.1f", ($reads_paired/$reads)*100);
			my $bases_mapped_perc   = sprintf("%.1f", ($bases_mapped/$total_bases)*100);
			my $rmdup_reads_mapped_perc   = sprintf("%.1f", ($rmdup_reads_mapped/$reads)*100);
			my $rmdup_bases_mapped_perc   = sprintf("%.1f", ($rmdup_bases_mapped/$total_bases)*100);
			my $clip_bases_perc    = sprintf("%.1f", ($clip_bases/$bases)*100);
			my $error_rate_perc    = sprintf("%.2f", ($error_rate*100));
			$reads_mapped = $utl->commify($reads_mapped);
			$bases_mapped = $utl->commify($bases_mapped);
			$reads_paired = $utl->commify($reads_paired);
			$rmdup_reads_mapped = $utl->commify($rmdup_reads_mapped);
			$rmdup_bases_mapped = $utl->commify($rmdup_bases_mapped);
			$bases = $utl->commify($bases);
			$raw_reads = $utl->commify($raw_reads);
			$raw_bases = $utl->bp_to_nearest_unit($raw_bases);
			$clip_bases = $utl->commify($clip_bases) if $clip_bases;

			if (defined $exome_design) {
				$bait_near_bases_mapped_perc   = sprintf("%.1f", 100 * $bait_near_bases_mapped / $total_bases);
				$target_near_bases_mapped_perc = sprintf("%.1f", 100 * $target_near_bases_mapped / $total_bases);
				$bait_bases_mapped_perc        = sprintf("%.1f", 100 * $bait_bases_mapped / $total_bases);
				$off_bait_bases_perc           = sprintf("%.1f", 100 * $off_bait_bases / $total_bases);
				$reads_on_bait_perc            = sprintf("%.1f", 100 * $reads_on_bait / $reads);
				$reads_on_bait_near_perc       = sprintf("%.1f", 100 * $reads_on_bait_near / $reads);
				$reads_on_target_perc          = sprintf("%.1f", 100 * $reads_on_target / $reads);
				$reads_on_target_near_perc     = sprintf("%.1f", 100 * $reads_on_target_near / $reads);
				$target_bases_mapped_perc      = sprintf("%.1f", 100 * $target_bases_mapped / $total_bases);

				$mean_target_coverage = sprintf("%.1f", $mean_target_coverage);
				$mean_bait_coverage = sprintf("%.1f", $mean_bait_coverage);
				$bait_near_bases_mapped = $utl->commify($bait_near_bases_mapped);
				$target_near_bases_mapped = $utl->commify($target_near_bases_mapped);
				$bait_bases_mapped = $utl->commify($bait_bases_mapped);
				$mean_bait_coverage = $utl->commify($mean_bait_coverage);
				$bait_coverage_sd = $utl->commify($bait_coverage_sd);
				$off_bait_bases = $utl->commify($off_bait_bases);
				$reads_on_bait = $utl->commify($reads_on_bait);
				$reads_on_bait_near = $utl->commify($reads_on_bait_near);
				$reads_on_target = $utl->commify($reads_on_target);
				$reads_on_target_near = $utl->commify($reads_on_target_near);
				$target_bases_mapped = $utl->commify($target_bases_mapped);
				$mean_target_coverage = $utl->commify($mean_target_coverage);
				$target_coverage_sd = $utl->commify($target_coverage_sd);
			}

			$reads = $utl->commify($reads);

			print qq[
				<div class="centerFieldset">
				<fieldset style="width: 800px">
				<legend>QC Plots</legend>
				<table width="80%">
				<tr>
				];

#map the image names to the object
			my $imglist = $mapstats->images;
			my %images = map {$_->name => $_} @$imglist;
			my @big_images;
			if ($exome_design) {
				@big_images = qw(exomeQC.gc_mapped.png exomeQC.gc_mapped.png exomeQC.insert_size.png);
			}
			else {
				@big_images = qw(gc-content.png insert-size.png);
			}

			foreach (@big_images) {
				my $img = $images{$_};
				if ($img) {
					print qq[<td>];printFullImage($img);print qq[</td>];
					delete $images{$_};
				}
			} 

			print qq[
				</tr>
				</table>
				<table>
				<tr>
				];

			my @i;

			if ($exome_design){
				@i = qw(exomeQC.quality_scores_1.png
						exomeQC.quality_scores_2.png
						exomeQC.gc_unmapped.png
						exomeQC.target_gc_vs_cvg.scaled.png
						exomeQC.cumulative_coverage.png
						exomeQC.coverage_per_base.png);
			} 
			else{
				@i = keys %images;
#@i = sort( { $a->name() cmp $b->name() } @{$images} );
			}

			foreach( @i )
			{
				print qq[<td>];printPreviewImage($images{ $_ });print qq[</td>];
			}

			print qq[
				</tr>
				</table>
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

			if (defined $clip_bases)
			{
				$clip_str = "$clip_bases ($clip_bases_perc\% raw)";
			}
			else 
			{  
				$clip_str = "unclipped";
			}
			$bmapped_str = "$bases_mapped ($bases_mapped_perc\%)" if $bases_mapped;
			$rmdup_bmapped_str = "$rmdup_bases_mapped ($rmdup_bases_mapped_perc\%)" if $rmdup_bases_mapped;


			my $mapping_section_title = $exome_design ? 'General mapping data' : 'Mapping data';

			print qq[
				<br/>
				<div class="centerFieldset">
				<fieldset style="width: 900px">
				<legend>$mapping_section_title</legend>
				<table width="100%">
				<tr>
				<td><table>
				<tr><td>Total Reads: </td><td>$raw_reads</td></tr>
				<tr><td>QC Reads: </td><td>$reads</td></tr>
				<tr><td>Reads w/adapter</td><td>$adapter_str</td></tr>
				<tr><td>Reads mapped: </td><td>$rmapped_str</td></tr>
				<tr><td>Reads paired: </td><td>$rpaired_str</td></tr>
				<tr><td>Reads mapped (rmdup): </td><td>$rmdup_rmapped_str</td></tr>
												</table>
														</td>
															 <td><table>
																		<tr><td>Total Bases: </td><td>$raw_bases</td></tr>
																						  <tr><td>QC Bases: </td><td>$bases</td></tr>
																									  <tr><td>Bases postclip: </td><td>$clip_str</td></tr>
																											   <tr><td>Bases mapped: </td><td>$bmapped_str</td></tr>
																													 <tr><td>&nbsp;</td><td>&nbsp;</td></tr>
																														 <tr><td>Bases mapped (rmdup): </td><td>$rmdup_bmapped_str</td></tr>
																																						 </table>
																																								 </td>
																																									  <td><table>
																																												 <tr><td>Mapper: </td><td>$mapper_string</td></tr>
																																																   <tr><td>Cycles: </td><td>$cycles</td></tr>
																																																			   <tr><td>NPG QC: </td><td>$npg_qc</td></tr>
																																																						<tr><td>Error rate: </td><td>$error_rate_perc%</td></tr>
																																																								 <tr><td>Genotype check: </td><td>$gt_display</td></tr>
																																																									   <tr><td>Duplication rate: </td><td>$dupe_rate_perc%</td></tr>
																																																										  </table>
																																																											 </td>
																																																												</tr>
																																																												   </table>
																																																												   </fieldset>
																																																												   </div>
																																																												   ];


			if ($exome_design) {
				print qq[
					<br/>
					<div class="centerFieldset">
					<fieldset style="width: 900px">
					<legend>Exome-specific data</legend>
					<table width="100%">
					<tr>
					<td><table>
					<tr><td>% target bases >= 1X:</td><td>$target_bases_1X</td></tr>
					<tr><td>% target bases >= 2X:</td><td>$target_bases_2X</td></tr>
					<tr><td>% target bases >= 5X:</td><td>$target_bases_5X</td></tr>
					<tr><td>% target bases >= 10X:</td><td>$target_bases_10X</td></tr>
					<tr><td>% target bases >= 20X:</td><td>$target_bases_20X</td></tr>
					<tr><td>% target bases >= 50X:</td><td>$target_bases_50X</td></tr>
					<tr><td>% target bases >= 100X:</td><td>$target_bases_100X</td></tr>
					<tr><td>Mean target coverage:</td><td>$mean_target_coverage</td></tr>
					<tr><td>Mean bait coverage:</td><td>$mean_bait_coverage</td></tr>
					<tr><td>Off bait bases:</td><td>$off_bait_bases ($off_bait_bases_perc%)</td></tr>
					</table></td>
					<td><table>
					<tr><td>Bases mapped on bait:</td><td>$bait_bases_mapped ($bait_bases_mapped_perc%)</td></tr>
					<tr><td>Bases mapped near bait:</td><td>$bait_near_bases_mapped ($bait_near_bases_mapped_perc%)</td></tr>
					<tr><td>Bases mapped on target:</td><td>$target_bases_mapped ($target_bases_mapped_perc%)</td></tr>
					<tr><td>Bases mapped near target:</td><td>$target_near_bases_mapped ($target_near_bases_mapped_perc%)</td></tr>
					<tr><td>Reads on bait:</td><td>$reads_on_bait ($reads_on_bait_perc%)</td></tr>
					<tr><td>Reads near bait:</td><td>$reads_on_bait_near ($reads_on_bait_near_perc%)</td></tr>
					<tr><td>Reads on target:</td><td>$reads_on_target ($reads_on_target_perc%)</td></tr>
					<tr><td>Reads near target:</td><td>$reads_on_target_near ($reads_on_target_near_perc%)</td></tr>
					</table></td>
					<td><table>
					<tr><td>Target bases:</t><td>$target_bases</td></tr>
					<tr><td>Bait bases:</t><td>$bait_bases</td></tr>
					</table></td>
					</tr>
					</table>
					</fieldset>
					</div>
					];
			}
		}
		print qq[
			<br/>
			<div class="centerFieldset">
			<fieldset style="background-color: $status_colour;width: 500px">
			<legend>Lane QC : $status &nbsp;&nbsp;&nbsp;(Auto QC : $auto_qc_status)</legend>
			];

		if ($utl->{AUTH_USERS}{$USER}) {   
			{
				my $laneID = $lane->id();
				my $LANE_UPDATE = $utl->{MODES}{LANE_UPDATE};
				print qq[
					<form action="$SCRIPT_NAME">
					<input type="hidden" name="mode" value="$LANE_UPDATE">
					<input type="hidden" name="db" value="$database">
					<input type="hidden" name="lane_id" value="$laneID">
					];

				my $PENDING = $utl->{STATES}{PENDING};
				my $GT_PENDING = $utl->{STATES}{GT_PENDING};
				my $INVESTIGATE = $utl->{STATES}{INVESTIGATE};
				my $FAILED = $utl->{STATES}{FAILED};
				my $PASSED = $utl->{STATES}{PASSED};
				print qq[
					<p>
					<table width="400" align="center">
					<tr>
					<td align="center"><input type="submit" name="lane_update" value="$PENDING" class="btn" style="background-color: #F5F5F5;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
					<td align="center"><input type="submit" name="lane_update" value="$GT_PENDING" class="btn" style="background-color: #F5F5F5;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
					<td align="center"><input type="submit" name="lane_update" value="$INVESTIGATE" class="btn" style="background-color: #F5F5F5;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
					<td align="center"><input type="submit" name="lane_update" value="$FAILED" class="btn" style="background-color: #FFC0C0;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
					<td align="center"><input type="submit" name="lane_update" value="$PASSED" class="btn" style="background-color: #C0FFC0;" onmouseover="this.className='btn btnhov'" onmouseout="this.className='btn'" /></td>
					</tr>
					</table>
					</p>
					</form>
					];
			}
			print qq[
				</fieldset>
				</div>
				];
		}
	}
}

# returns the project, sample, library objects for a given lane
sub get_hierarchy_for_lane 
{
    my $vrtrack = shift;
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

    unless ($proj && $samp && $lib)
    {
        print_and_exit(sprintf('Cannot retrieve hierarchy for lane %s', $lane->hierarchy_name));
    }
    return ($proj, $samp, $lib);
}

sub printFullImage
{
    my $imageObj = shift;
    
    $imageObj->name =~ /\.(\w+)$/;
    my $type = $1;
    my $uri = URI->new("data:");
    $uri->media_type("image/$type");
    $uri->data($imageObj->image);
    my $caption = $imageObj->caption;
    
    print qq[
        <img src="$uri" alt="$caption" width="400" height="400" style="border:1px dotted #83A4C3;">
    ];
}

sub printPreviewImage
{
    my $imageObj = shift;
    
    $imageObj->name =~ /\.(\w+)$/;
    my $type = $1;
    my $uri = URI->new("data:");
    $uri->media_type("image/$type");
    $uri->data($imageObj->image);
    my $caption = $imageObj->caption;
    
    print qq[
        <a class="thumbnail" href="#">
        <img src="$uri" width="100" height="100" alt="$caption" style="border:1px dotted #83A4C3;">
        <span><img src="$uri" alt="$caption" style="border:1px dotted #83A4C3;"/></span></a>
    ];
}
