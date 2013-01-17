#!/usr/local/bin/perl -T

BEGIN {
    $ENV{VRTRACK_HOST} = 'mcs10';
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
use VertRes::QCGrind::ViewUtil;

my $utl = VertRes::QCGrind::ViewUtil->new();

my $title = 'Pending View';

my $sw  = SangerWeb->new({
    'title'   => $title,
    'banner'  => q(),
    'inifile' => SangerWeb->document_root() . q(/Info/header.ini),
    'jsfile'  => ['http://jsdev.sanger.ac.uk/prototype.js','http://jsdev.sanger.ac.uk/toggle.js','ttp://jsdev.sanger.ac.uk/scriptaculous/scriptaculous.js','http://jsdev.sanger.ac.uk/sidebar.js','http://jsdev.sanger.ac.uk/urchin.js','http://jsdev.sanger.ac.uk/zebra.js','http://js.sanger.ac.uk/sorttable_v2.js',],  
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

my $dbpend = $cgi->param('dbpend');
print $sw->header();
displayPendingRequestsPage($cgi,$vrtrack,$dbpend);
print $sw->footer();
exit;

#######

sub displayPendingRequestsPage 
{

    my $cgi = shift;
    my $vrtrack = shift;
    my $db = shift;
    
    my $lib_type = "library";
    my $seq_type = "sequence";
    my $plex_type = "multiplex";

    my $init_script = $utl->{SCRIPTS}{PENDING_VIEW};
    my $index = $utl->{SCRIPTS}{INDEX_PAGE};
   
    my @libpendings = fetchRequests($vrtrack, $utl->{SQL}{PENDING_LIB}, $db, $lib_type);
    my @seqpendings = fetchRequests($vrtrack, $utl->{SQL}{PENDING_SEQ}, $db, $seq_type);
    my @mplexpendings = fetchRequests($vrtrack, $utl->{SQL}{PENDING_MULTIPLEX_SEQ}, $db, $plex_type);
    
    print qq[
        <h4 align="center" style="font: arial"><i><a href="$index">Team 145</a> : <a href="$init_script">$title</a></i> : $db</h4>
		<div class="centerFieldset">
    ];
	printPendingRows(\@libpendings, $lib_type, $db, $cgi);
	printPendingRows(\@seqpendings, $seq_type, $db, $cgi);
	printPendingRows(\@mplexpendings, $plex_type, $db, $cgi);
    print qq[
    	</div>
    ];
}

sub fetchRequests
{
  my $vrtrack = shift;
  my $sql_fetch_pending = shift;
  my $db = shift;
  my $type = shift;
  my @pending;
  my @started;
  my $sth = $vrtrack->{_dbh}->prepare($sql_fetch_pending);
  if ($sth->execute($db)) {
    my ($project_name, $sample_name, $seq_status, $changed, $ssid); 
    $sth->bind_columns(\($project_name, $sample_name, $seq_status, $changed, $ssid));
    while ($sth->fetch) {
		$seq_status eq 'pending' ? push @pending, [$project_name, $sample_name, $seq_status, $changed, $ssid]: push @started, [$project_name, $sample_name, $seq_status, $changed, $ssid];
    }
  } 
  if ($type eq 'multiplex') {
	@pending = sort {$a->[4] <=> $b->[4]} @pending;
	@started = sort {$a->[4] <=> $b->[4]} @started;
  }
  return (@pending, @started);
}

sub printPendingRows
{
 	my ($pendinglist, $req_type, $db, $cgi) = @_;
    my $qcgrind_script = $utl->{SCRIPTS}{QCGRIND_SAMPLES};
    print qq[
    <fieldset > 
    <legend>Pending $req_type requests</legend>
    <table class='sortable' width="80%">
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
			if ($current_project ne $pending[0]) {
				print qq[<td><a href="$qcgrind_script?db=$db&amp;proj_id=$pending[0]">QC_grind</a></td>];
				$current_project = $pending[0];
			}
			else {
				print qq[<td></td>];
			}			
			for my $j ( 0 .. ($#pending-1) ) {
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
