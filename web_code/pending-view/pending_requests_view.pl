#!/usr/local/bin/perl -T

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
use VRTrack::VRTrack;
use VRTrack::Project;
use VRTrack::Sample;
use VRTrack::Library;
use VRTrack::Lane;
use VRTrack::Multiplex_pool;
use VertRes::Utils::VRTrackFactory;
use VertRes::QCGrind::ViewUtil;

my $utl = VertRes::QCGrind::ViewUtil->new();

my $title = 'Pending View';

my $sw  = SangerWeb->new({
    'title'   => $title,
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
    'style'   => $utl->{CSS},
});

my $cgi = $sw->cgi();

my $db = $cgi->param('db');
unless ($db) {
    print $sw->header();
	$utl->displayError( "No database ID",$sw );
}
my $vrtrack = $utl->connectToDatabase($db);
unless ( defined $vrtrack ) {
	print $sw->header();
	$utl->displayError( "No database connection to $db",$sw );
}

print $sw->header();
displayPendingRequestsPage($cgi,$vrtrack,$db);
print $sw->footer();
exit;

#######


sub displayPendingRequestsPage 
{
    my $cgi = shift;
    my $vrtrack = shift;
    my $db = shift;

    my $init_script = $utl->{SCRIPTS}{PENDING_VIEW};
    
	my @projects = sort {$a->name cmp $b->name} @{$vrtrack->projects()};
	
	my @libpendings;
	my @libstarted;
	my @seqpendings;
	my @seqstarted;
	my @mplexpendings;
	my @mplexstarted;
	my $lib_type = "library";
	my $seq_type = "sequence";
	my $plex_type = "multiplex sequence";

	print qq[
    	<h2 align="center" style="font: normal 900 1.5em arial"><a href="$init_script">Pending Requests</a></h2>
    	
		<div class="centerFieldset">
    ];

    foreach my $project (@projects)
    {
        my $projectID = $project->id();
    
        my $samples = $project->samples();
        displayError( "Cant get samples for project: $projectID" ) unless $samples;
     
        my $pname = $project->name;

        my $pendingflag = 'pending';
        my $startedflag = 'started';

        foreach( @{$samples} ) 
        {	
		    my $sample = $_;
	    	my $sname = $sample->name();
		    my $library_requests = $sample->library_requests();
		    foreach ( @{$library_requests}) {
			    my $librequest = $_;
			    my $lib_status = $librequest->prep_status();
			    my $lib_date = $librequest->changed();
			    my $ssid = $librequest->ssid();
			    if ($lib_status eq $pendingflag) { push @libpendings, [$projectID, $pname, $sname, $lib_status, $lib_date, $ssid]; }
			    elsif ($lib_status eq $startedflag) { push @libstarted, [$projectID, $pname, $sname, $lib_status, $lib_date, $ssid]; }
		    }
		    my $libraries = $sample->libraries();
		    foreach ( @{$libraries}) {
			    my $library = $_;
			    my $seq_requests = $library->seq_requests();
			    foreach ( @{$seq_requests} ) {
				    my $seqrequest = $_;
				    my $seq_status = $seqrequest->seq_status();
				    my $seq_date = $seqrequest->changed();
				    my $ssid = $seqrequest->ssid();				
				    if ($seq_status eq $pendingflag) { push @seqpendings, [$projectID, $pname, $sname, $seq_status, $seq_date, $ssid];}
				    elsif ($seq_status eq $startedflag) { push @seqstarted, [$projectID, $pname, $sname, $seq_status, $seq_date, $ssid];}
			    }
			    foreach ( @{$library->library_multiplex_pools}){
				    my $mplex = VRTrack::Multiplex_pool->new($vrtrack, $_->multiplex_pool_id);
                    foreach ( @{ $mplex->seq_requests } ){
                	    my $seqrequest = $_;
                	    my $seq_status = $seqrequest->seq_status();
                	    my $seq_date = $seqrequest->changed();
					    my $ssid = $seqrequest->ssid();
					    if ($seq_status eq $pendingflag) { push @mplexpendings, [$projectID, $pname, $sname, $seq_status, $seq_date, $ssid]; }
					    if ($seq_status eq $startedflag) { push @mplexstarted, [$projectID, $pname, $sname, $seq_status, $seq_date, $ssid]; }
                    }
                }
                @mplexpendings = sort {$a->[5] <=> $b->[5]} @mplexpendings;
                @mplexstarted = sort {$a->[5] <=> $b->[5]} @mplexstarted;
		    }
        }
	}
	my @libtotal = (@libpendings, @libstarted);
	my @seqtotal = (@seqpendings, @seqstarted);
	my @mplextotal = (@mplexpendings, @mplexstarted);
	printPendingRows(\@libtotal, $lib_type, $db, $cgi);
	printPendingRows(\@seqtotal, $seq_type, $db, $cgi);
	printPendingRows(\@mplextotal, $plex_type, $db, $cgi);
    print qq[
    	</div>
    ];
}

sub printPendingRows
{
 	my ($pendinglist, $req_type, $db, $cgi) = @_;
    my $qcgrind_script = $utl->{SCRIPTS}{QCGRIND_SAMPLES};
    print qq[
    <fieldset > 
    <legend>Pending $req_type requests</legend>
    <table width="80%">
    <tr>
    <th>QC Link</th>
    <th>Project</th>
    <th>Sample Name</th>
	<th>Status</th>
	<th>Date changed</th>
	<th>SSID</th>
	</tr>
    ];
	my @pendings = @$pendinglist;
	my $current_project = '';
  	if (@pendings) {
		for my $i ( 0 .. $#pendings ) {
			my @pending = @{$pendings[$i]};
			print qq[<tr>];
			if ($current_project ne $pending[1]) {
				print qq[<td><a href="$qcgrind_script?db=$db&amp;proj_id=$pending[0]">QC_grind</a></td>];
				$current_project = $pending[1];
			}
			else {
				print qq[<td></td>];
			}			
			for my $j ( 1 .. ($#pending-1) ) {
				print qq[<td>$pending[$j]</td>];
			}
			print qq[<td><a href="http://psd-production.internal.sanger.ac.uk:6600/requests/$pending[$#pending]">$pending[$#pending]</a></td></tr>];
		}
	}
  	else {
		print qq[<tr><td colspan=5>No pending $req_type requests were found in the $db database.</td></tr>];
	}	
    print qq[
        </table>
        </fieldset><br/>
    ];
}
