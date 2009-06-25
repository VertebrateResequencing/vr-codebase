#!/usr/local/bin/perl

use strict;
use warnings;
no warnings 'uninitialized';
use Getopt::Long;
use DBI;
use constant DBI_DUPLICATE => '1062';


my ($help,$imglist,$metadata, $spp);

GetOptions(
    'h|help'	    =>  \$help,
    'i|images=s'    =>  \$imglist,
    'm|meta=s'	    =>  \$metadata,
    's|species=s'   =>  \$spp,
    );

(-f $imglist && -f $metadata && $spp && !$help) or die <<USAGE;
    Usage: $0 
                --i    <find list of images to load>
                --s    <species - g1k or mouse>
                --m    <metadata file>
		--help <this message>

Loads QC images into qc tracking database.  

Image list must be a hierarchy path (starting at the project) down to the qc
image, e.g. 
./129S1_SvImJ_Mouse_Genome/1/SLX/129S1_SvIMJ_SLX_200_NOPCR_2/2024_3/qc/gc.gif

Images must also be rsynced over to the web image server:
webdatasrv:/mitocheck/images/qc for display

Metadata is Thomas Keane's mapping data: csv with fields:
Status,Project,Individual,Tech,Library,Lane,Fastq0,ReadLength,Fastq1,ReadLength,Fastq2,ReadLength,#Reads,#Bases,#RawReadsMapped,#RawBasesMapped,#RawReadsPaired,#RmdupReadsMapped,#RmdupBasesMapped,ErrorRate

2009-04-02 changed to:
Status,Center,Project,Individual,Tech,Library,Lane,Fastq0,ReadLength,Fastq1,ReadLength,Fastq2,ReadLength,#Reads,#Bases,#RawReadsMapped,#RawBasesMapped,#ReadsPaired,ErrorRate

USAGE

my %db = (  'g1k'   => 'g1k_track',
	    'mouse' => 'mouse_reseq_track',
	);

my $db = $db{$spp};
$db or die "Can't find database for $spp\n";

my $dbh = DBI->connect("DBI:mysql:host=mcs4a;database=$db",
			"vreseq_rw", "t3aml3ss",  
                      {'RaiseError' => 0, 'PrintError'=>0}
		      );
if ($DBI::err){
      die(sprintf('DB connection failed: %s', $DBI::errstr));
}

# hash the metadata for including with each image
open (my $META, $metadata)or die "Can't open metadata file: $!\n"; 

my %metaidx = ( 'status'	=> 0,
		'center'	=> 1,
		'project'	=> 2,
		'sample'	=> 3,
		'library'	=> 5,
		'lane'		=> 6,
		'readlen0'	=> 8,
		'readlen1'	=> 10,
		'raw_reads'	=> 13,
		'raw_bases'	=> 14,
		'reads_mapped'	=> 15,
		'bases_mapped'	=> 16,
		'reads_paired'	=> 17,
		'rmdup_reads_mapped'	=> 18,
		'rmdup_bases_mapped'	=> 19,
		'error_rate'	=> 20,
	    );

my %meta;
while (<$META>){
    chomp;
    my @data = split /, */, $_;
    next unless $data[$metaidx{'status'}] eq 'MAPPED';
    my $project = $data[$metaidx{'project'}];
    my $sample = $data[$metaidx{'sample'}];
    my $library = $data[$metaidx{'library'}];
    my $lane = $data[$metaidx{'lane'}];
    # load hash with all the fields from @data specified in %metaidx
    my %hash = map { $_ => $data[$metaidx{$_}] } keys %metaidx;

    # single-end runs will have a readlen for fastq0 (readlen0), while paired will
    # have a readlen1.  Pick one for 'readlen':
    $hash{'readlen'} = $hash{'readlen0'} || $hash{'readlen1'};

    # store that hash under the hierarchy key:
    $meta{$project}{$sample}{$library}{$lane} = \%hash;
}

close ($META);


open (my $IMG, $imglist) or die "Can't open image list: $!\n";

my $imgcount = 0;
while (<$IMG>){
    chomp;
    s|^./||;	# trim leading dot directory from find command
    my ($project, $sample, $tech, $library, $lane, $qc, $img) = split "/",$_;
    die "Can't parse $_ \n" unless $img;

    # HACK FOR MOUSE SAMPLE NAMES jws 20090624
    # mouse sample names in hierarchy are all '1'.  Set to project name, minus
    # the '_Mouse_Genome' suffix.

    if ($spp eq 'mouse'){
        $sample = $project;
        $sample =~ s/_Mouse_Genome$//;
    }
    # EOHACK
    
    my $lanemeta = $meta{$project}{$sample}{$library}{$lane};
    unless ($lanemeta){
	warn "No metadata for $lane\n";
	next;
    }
    # query/update database for new data
    # NB only data with images get added to the database
    my $proj_id = fetch_project_id($project);
    my $sample_id = fetch_sample_id($proj_id, $sample);
    my $library_id = fetch_library_id($sample_id, $library);
    my $lane_id = fetch_lane_id($library_id, $lane, $lanemeta);

    my $caption;
    if ($img =~/gc/){
	$caption = "Percentage GC";
    }
    elsif ($img =~ /insert/){
	$caption = "Insert Size Distribution";
    }
    else {
	$caption = "";
    }

    my $path = join "/", ($project, $sample, $tech, $library, $lane, $qc, $img);
    $imgcount += add_image($lane_id, $img, $path, $caption);
}

close($IMG);

print "$imgcount images loaded\n";

###############################################################################
###############################################################################
###############################################################################

sub add_image {
    my ($lane_id, $name, $path, $caption) = @_;

    my $sql = qq[INSERT INTO img (lane_id, name, caption, path) VALUES (?,?,?,?)];
    my $sth = $dbh->prepare($sql);
    my $loaded = 0;
    if ($sth->execute( $lane_id, $name, $caption,$path )) {
	# success
	$loaded = 1;
    }
    else {
	if ($DBI::err eq DBI_DUPLICATE) {
	    #warn "Image $name is already present in the database\n";
	}
	else {
	    die( sprintf('DB load insert failed: %s %s', $name, $DBI::errstr));
	}
    }
    return $loaded;

}


# returns existing lane id or adds lane and returns new id
# 2009-04-28 updates existing data if it has changed.
sub fetch_lane_id {
    my ($library_id, $name, $meta) = @_;
    die "No metadata for $name" unless $meta;	# this is trapped earlier anyway

    my @metafields =qw(readlen raw_reads raw_bases reads_mapped reads_paired bases_mapped rmdup_reads_mapped rmdup_bases_mapped error_rate);

    my $sql = qq[select * from lane where library_id = ? and name=?];
    my $id_ref = $dbh->selectrow_hashref($sql, undef, ($library_id,$name));
    my $id;
    if ($id_ref){
	$id = $id_ref->{lane_id};
	my @changes;
	foreach my $field(@metafields){
	    next unless $meta->{$field}; # don't blow away existing data if the 
					 #summary  has lost it for some reason.
	    if ($id_ref->{$field} != $meta->{$field} ){
		warn "$name.$field has changed. Was ".$id_ref->{$field}." now ".$meta->{$field}.".  Updating.\n";
		push @changes, "$field=".$meta->{$field};
	    }
	    if (@changes){
		my $updatesql = qq[UPDATE lane set changed=now(), ].join ", ", @changes;
		$updatesql .= qq[ where library_id = ? and name=?];

		my $sth = $dbh->prepare($updatesql);
		if ($sth->execute( $library_id, $name)) {
		    warn "Updated\n";
		}
		else {
		    die( sprintf('DB update failed: %s %s', $name, $DBI::errstr));
		}

	    }

	}
    }
    else {
	my $metafields = join (", ", @metafields);
	my $metaplace = "?," x scalar @metafields;
	$sql = qq[INSERT INTO lane (lane_id, library_id, name, $metafields, qc_status, changed) VALUES (NULL,?,?,$metaplace 'pending',now())];
	my $sth = $dbh->prepare($sql);
	if ($sth->execute( $library_id, $name, @$meta{@metafields} )) {
	    $id = $dbh->{'mysql_insertid'}
	}
	else {
	    if ($DBI::err eq DBI_DUPLICATE) {
		die "Lane $name is already present in the database\n";
	    }
	    else {
		die( sprintf('DB load insert failed: %s %s', $name, $DBI::errstr));
	    }
	}
    }
    return $id;
}


# returns existing lib id or adds lib and returns new id
sub fetch_library_id {
    my ($sample_id, $name) = @_;

    my $sql = qq[select library_id from library where sample_id = ? and name=?];
    my $id_ref = $dbh->selectrow_hashref($sql, undef, ($sample_id,$name));
    my $id;
    if ($id_ref){
	$id = $id_ref->{library_id};
    }
    else {

	$sql = qq[INSERT INTO library (library_id, sample_id, name, qc_status) VALUES (NULL,?,?,'pending')];
	my $sth = $dbh->prepare($sql);
	if ($sth->execute( $sample_id, $name )) {
	    $id = $dbh->{'mysql_insertid'}
	}
	else {
	    if ($DBI::err eq DBI_DUPLICATE) {
		die "INSERT ERROR: Library $name is already present in the database\n";
	    }
	    else {
		die( sprintf('DB load insert failed: %s %s', $name, $DBI::errstr));
	    }
	}
    }
    return $id;
}


# returns existing sample id or adds sample and returns new id
sub fetch_sample_id {
    my ($proj_id, $name) = @_;

    my $sql = qq[select sample_id from sample where project_id = ? and name=?];
    my $id_ref = $dbh->selectrow_hashref($sql, undef, ($proj_id,$name));
    my $id;
    if ($id_ref){
	$id = $id_ref->{sample_id};
    }
    else {
	$sql = qq[INSERT INTO sample (sample_id, project_id, name) VALUES (NULL,?,?)];
	my $sth = $dbh->prepare($sql);
	if ($sth->execute( $proj_id, $name )) {
	    $id = $dbh->{'mysql_insertid'}
	}
	else {
	    if ($DBI::err eq DBI_DUPLICATE) {
		die "INSERT ERROR: Sample $name is already present in the database\n";
	    }
	    else {
		die( sprintf('DB load insert failed: %s %s', $name, $DBI::errstr));
	    }
	}
    }
    return $id;
}



# returns existing project id or adds project and returns new id
sub fetch_project_id {
    my $name = shift;

    my $sql = qq[select project_id from project where name=?];
    my $id_ref = $dbh->selectrow_hashref($sql, undef, ($name));
    my $proj_id;
    
    if ($id_ref){
	$proj_id = $id_ref->{project_id};
    }
    else {
	$sql = qq[INSERT INTO project (project_id, name) VALUES (NULL,?)];
	my $sth = $dbh->prepare($sql);
	if ($sth->execute( $name )) {
	    $proj_id = $dbh->{'mysql_insertid'}
	}
	else {
	    if ($DBI::err eq DBI_DUPLICATE) {
		die "INSERT ERROR: Project $name is already present in the database\n";
	    }
	    else {
		die( sprintf('DB load insert failed: %s %s', $name, $DBI::errstr));
	    }
	}
    }
    return $proj_id;
}

