#!/usr/bin/env perl

=head1 NAME

pro_load.pl - load a proserver mysql database from a GFF file

=head1 DESCRIPTION

Loads a proserver database, either via DBI or bulk to a local db, from a GFF
file.  

Largely copied from the ldas_bulk_load.pl script by Lincoln Stein.

=head1 AUTHOR

Jim Stalker (jws@sanger.ac.uk) (with thanks to Lincoln Stein)

=cut

use strict;

use DBI;
use IO::File;
use Getopt::Long;

use constant FEATURETABLE   => 'feature';
use constant GROUPTABLE	    => 'fgroup';

# could do something clever to generate the schema from this, or at least check
# that it matches the schema in <DATA> at EOF
my @FEATURE_COLUMNS = qw(id label segment start end score orient phase type_id type_category method group_id target_id target_start target_end link_url link_text note);

my @GROUP_COLUMNS = qw(group_id group_label group_type group_link_url group_link_text group_note);

my ($DBNAME,$HOST,$PORT,$USER,$PASS,$CREATE,$BULK, $MYSQL,$WIPE);

GetOptions ('database:s'    => \$DBNAME,
	    'host:s'	    => \$HOST,
	    'port:s'	    => \$PORT,
	    'user:s'        => \$USER,
	    'password:s'    => \$PASS,
	    'create'        => \$CREATE,
	    'wipe'          => \$WIPE,
	    'bulk'	    => \$BULK,
	    'mysql:s'	    => \$MYSQL,
	   );
	   
	   
$DBNAME or die <<USAGE;
Usage: $0 [options] <gff file 1> <gff file 2> ...
Load a DAS database from GFF.

 Options:
     --database		   MySQL database name
     --host		   MySQL database host
     --port		   MySQL database port
     --user                MySQL username
     --password            MySQL password
     --wipe                Wipe database entries first.
     --create              Create database (or empty existing database)
     --bulk		   Do bulk load (local db only)
     --mysql		   Location of mysql client binary (if not in PATH)

NOTE: If no arguments are provided, then the input is taken from
standard input. Compressed files (.gz, .Z, .bz2) are automatically
uncompressed.

The nature of the bulk load requires that the database be on the local
machine and that the indicated user have the "file" privilege to load
the tables and have enough room in /tmp (or whatever is specified
by the \$TMPDIR environment variable), to hold the tables transiently.

USAGE
;

$MYSQL ||= 'mysql';

foreach (@ARGV) {
      $_ = "gunzip -c $_ |" if /\.gz$/;
      $_ = "uncompress -c $_ |" if /\.Z$/;
      $_ = "bunzip2 -c $_ |" if /\.bz2$/;
}

my %FH;
my $DBH;

# Check we can connect 
my $params = '';
$params .= " -u$USER"   if defined $USER;
$params .= " -p$PASS"   if defined $PASS;
unless ($BULK){
    # Can only bulk load on localhost
    $params .= " -h$HOST"   if defined $HOST;
    $params .= " -P$PORT"   if defined $PORT;
};

my $command = qq($MYSQL $params -e ";");
if (system $command){
    die "Can't connect to database server\n";
}

# Create new DB/clean out existing DB
if ($CREATE) {
    # check if we have db first
    $command = qq($MYSQL $params -e "use $DBNAME" > /dev/null 2>&1);
    if (system $command){
	# we got an error, so database doesn't exist
	warn "\nCreating database $DBNAME.\n";
	$command = qq($MYSQL $params -f -e "create database $DBNAME");
	system $command;
    }
    else {
	# db exists, so drop tables
	warn "\nCleaning database $DBNAME.\n";
	$command = qq($MYSQL $params -f -e "drop table ${\FEATURETABLE}" $DBNAME);
	system $command;

	$command = qq($MYSQL $params -f -e "drop table ${\GROUPTABLE}" $DBNAME);
	system $command;
    }

    # load db schema
    open (M,"| $MYSQL $params $DBNAME -f -B") or die "Couldn't open Mysql: $!";
    while (<DATA>) {
	print M;
    }
    close M;

}


my $tmpdir = $ENV{TMPDIR} || $ENV{TMP} || '/tmp';

if ($BULK){
    die "Option --wipe doesn't work with --bulk yet...." if $WIPE;
    foreach (FEATURETABLE, GROUPTABLE) {
	$FH{$_} = IO::File->new("$tmpdir/$_.$$",">") or die $_,": $!";
	$FH{$_}->autoflush;
    }
}
else {
    my $dsn = "DBI:mysql:database=$DBNAME;";
    $dsn .= "host=$HOST;" if $HOST;
    $dsn .= "port=$PORT;" if $PORT;
    $DBH = DBI->connect($dsn, $USER, $PASS, {
			PrintError => 0,   ### Don't report errors via warn(  )
			RaiseError => 1    ### Do report errors via die(  )
		  });
    if ($WIPE){
      warn "Deleting all current data.\n";
      $DBH->do("DELETE FROM feature") or warn "Could not delete from feature table";
      $DBH->do("DELETE FROM fgroup") or warn "Could not delete from fgroup table";
    }
}

my %fid_done;
my %gid_done;

my $total_count	    = 0;
my $feat_count	    = 0;
my $featload_count  = 0;
my $group_count	    = 0;
my $groupload_count = 0;

while (<>) {
    chomp;
    next if /^\#/;
    $total_count++;
    my ($segment,$method,$type,$start,$end,$score,$orient,$phase,$attribs) = split "\t";
   
    # Quick hacky check to see if this looks OK-ish.
    unless ($start =~ /^\d+$/ or $end =~ /^\d+$/){
	die "Can't find start & end - this doesn't look like a valid GFF input file (line $total_count)\n";
    }

    $score  = undef if $score  eq '.';
    $orient = undef if $orient eq '.';
    $phase  = '-' if $phase  eq '.';
    $orient = '+' if $orient eq '1';
    $orient = '-' if $orient eq '-1';

    my $feature = {	
		    id		    => undef,
		    label	    => undef,
		    segment	    => $segment,
		    start	    => $start,
		    end		    => $end,
		    score	    => $score,
		    orient	    => $orient,
		    phase	    => $phase,
		    type_id	    => $type,
		    type_category   => undef,
		    method	    => $method,
		    group_id	    => undef,
		    group_label	    => undef,
		    group_type	    => undef,
		    group_note	    => undef,
		    group_link_url  => undef,
		    group_link_text => undef,
		    target_id	    => undef,
		    target_start    => undef,
		    target_end	    => undef,
		    link_url	    => undef,
		    link_text	    => undef,
		    note	    => undef,
		    };


    # handle attribs parsing
    # Far too simplistic, unfortunately - needs complex regex to do properly.
    #$attribs =~ s/(\"[^\"]*);([^\"]*\")/$1$;$2/g;  # protect embedded semicolons in the attribs
    my @attribs = split(/\s+;\s+/,$attribs);

    foreach (@attribs) { 
	#s/$;/;/g;
	my ($tag , $value) = /^(\S+)\s*(.*)/;
	$tag ||= 'note';
	$tag = lc($tag);
	$value =~ s/\\t/\t/g;
	$value =~ s/\\r/\r/g;
	$value =~ s/\\n/\n/g;
	$value =~ s/^"//;
	$value =~ s/"$//;
	if ($feature->{$tag}){
	    #warn "$tag value override from attribs\n" unless $feature->{$tag} eq $value;
	}
	$feature->{$tag} = $value;
    }

    # Fill in default values
    
    # ID

    # get id specified in attribs, if any
    my $id = $feature->{'id'};
    $feature->{'group_id'} ||= $id;
   
    # if no id specified, make one up.
    $id ||= "$segment:$start-$end";
    my $orig_id = $id;

    # if necessary, uniquify id with a suffix.
    while ($fid_done{$id}++){
	$id = "$id.".$fid_done{$id};
    }
    $feature->{'id'} = $id;
    
    # LABEL
    $feature->{'label'} ||= $orig_id;	# set label to possibly non-unique id

    # GROUP
    # if no group id or feature id specified in attribs, fallback to grouping
    # on unique id.  It's a group of one, obviously...
    $feature->{'group_id'} ||= $id;

    if ( $total_count % 1000 == 0) {
      print STDERR "$total_count features parsed...\n";
    }
    
    if ($BULK){
	print_feature($feature);
	print_group($feature) unless $gid_done{$feature->{'group_id'}}++;
    }
    else {
	dbwrite_feature($feature);
	dbwrite_group($feature) unless $gid_done{$feature->{'group_id'}}++;
    }

} # end loop over infiles

# Do bulk loading
if ($BULK && $total_count){
    $_->close foreach values %FH;
    warn "Bulk loading...\n";

    my $AUTH = '';
    $AUTH .= " -u$USER"     if defined $USER;
    $AUTH .= " -p$PASS" if defined $PASS;

    foreach (FEATURETABLE, GROUPTABLE) {
	my $command = qq($MYSQL $AUTH -e "lock tables $_ write; delete from $_; load data infile '$tmpdir/$_.$$' replace into table $_; unlock tables" $DBNAME);
	  $command =~ s/\n/ /g;
	  system $command;
	  unlink "$tmpdir/$_";
    }
}

warn "$total_count features parsed\n";
unless ($BULK){
    warn "$featload_count out of $feat_count features loaded\n";
    warn "$groupload_count out of $group_count groups loaded\n";
}

#flush tables - give us new data for our new selects please
warn "\nFlushing tables in RDBMS.\n";
$command = qq($MYSQL $params -f -e "flush tables" $DBNAME);
system $command;
#warn "\nFlushing tables in database $DBNAME.\n";
#$command = qq($MYSQL $params -f -e "flush table ${\FEATURETABLE}" $DBNAME);
#system $command;
#$command = qq($MYSQL $params -f -e "flush table ${\GROUPTABLE}" $DBNAME);
#system $command;
    

###############################################################################

sub print_feature {
    my $feature = shift;
    $FH{ FEATURETABLE()}->print( join("\t", map ($feature->{$_}, @FEATURE_COLUMNS)));
    $FH{ FEATURETABLE()}->print("\n");
}

sub print_group {
    my $feature = shift;
    $FH{ GROUPTABLE()}->print( join("\t", map ($feature->{$_}, @GROUP_COLUMNS)));
    $FH{ GROUPTABLE()}->print("\n");
}

sub dbwrite_feature {
    my $feature = shift;
    $feat_count++;
    my $sql = "insert into ".FEATURETABLE." (";
    $sql .= join (",", @FEATURE_COLUMNS);
    $sql .= ") values (";
    $sql .= join (",",split("","?" x scalar(@FEATURE_COLUMNS)));
    $sql .= ")";

    my $sth = $DBH->prepare ($sql);
    eval {
	$sth->execute(map ($feature->{$_}, @FEATURE_COLUMNS));
    };
    if($@){
	warn "Error inserting ".$feature->{'id'}. " into ${\FEATURETABLE} table: ".$sth->errstr."\n";
    }
    else {
	$featload_count++;
    }
    $sth->finish;
}

sub dbwrite_group {
    my $feature = shift;
    $group_count++;
    my $sql = "insert into ".GROUPTABLE." (";
    $sql .= join (",", @GROUP_COLUMNS);
    $sql .= ") values (";
    $sql .= join (",",split("","?" x scalar(@GROUP_COLUMNS)));
    $sql .= ")";

    my $sth = $DBH->prepare ($sql);
    eval {
	$sth->execute(map ($feature->{$_}, @GROUP_COLUMNS));
    };
    if($@){
	warn "Error inserting ".$feature->{'id'}. " into ${\GROUPTABLE} table: ".$sth->errstr."\n";
    }
    else {
	$groupload_count++;
    }
    $sth->finish;
}

###############################################################################



__END__

CREATE TABLE feature (
    id		    varchar(60) NOT NULL,
    label	    varchar(60) default NULL,
    segment	    varchar(60) NOT NULL,
    start	    int(11) unsigned NOT NULL default '0',
    end		    int(11) unsigned NOT NULL default '0',
    score	    float default NULL,
    orient	    enum('0','+','-') default '0',
    phase	    enum('0','1','2','-') default '-',
    type_id	    varchar(240) NOT NULL,
    type_category   varchar(60) default NULL,
    method	    varchar(60) NOT NULL,
    group_id	    varchar(60) default NULL,
    target_id	    varchar(60) default NULL,
    target_start    int(11) unsigned default NULL,
    target_end	    int(11) unsigned default NULL,
    link_url	    varchar(255) default NULL,
    link_text	    varchar(60) default NULL,
    note	    text,
    UNIQUE KEY	    id_key (id,segment,start,end),
    KEY segment_key (segment,start)
);

CREATE TABLE fgroup (
    group_id	    varchar(60) NOT NULL,
    group_label	    varchar(60) default NULL,
    group_type	    varchar(60) default NULL,
    group_link_url  varchar(255) default NULL,
    group_link_text varchar(60) default NULL,
    group_note	    text,
    PRIMARY KEY	    group_id_key (group_id)
);

