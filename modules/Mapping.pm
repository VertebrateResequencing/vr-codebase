package Mapping;
use strict;
use Carp;
use File::Basename;
use File::Spec;
use Cwd qw(getcwd abs_path);

use Mapping_SLX_Maq;
use Mapping_454_ssaha;
use HierarchyUtilities;

our $LOCAL_CACHE_DIR = '/localcache';

=pod

=head1 NAME

MapHierarchy

=head1 SYNOPSIS

A module that is that is used to map on the vert reseq hierarchy

=head1 REQUIRES

Perl5.8.8

=head1 DESCRIPTION

A module for carrying out mapping on the Vertebrate Resequencing Informatics 
mapping hierarchy. The hierarchy can consist of multiple different sequencing 
technologies.

=head1 METHODS

=head2 map

	Arg [1]    : root of the mapping hierarchy
	Arg [2]    : LSF queue for jobs
	Arg [3]    : file of lanes to be mapped
	Arg [3]    : index file (empty string if none)
	Example    : map( '$G1K/MOUSE/MAPPING', 'normal', 'lanes.fofn', 'sequence.index');
	Description: Maps the lanes specified in the lanes.fofn file
	Returntype : none
=cut

sub map
{
	croak "Usage: rootDir lsf_queue lanes_fofn index_file" unless @_ == 4;
	
	my $mRoot = shift;
	my $lsf_queue = shift;
	my $lanes_fofn = shift;
	my $indexF = shift;
	
	croak "Cant find the lanes_fofn file: $lanes_fofn\n" unless -f $lanes_fofn;
	
	if( ! -f $indexF )
	{
		print "Cant find index file: $indexF - ignoring index file\n";
		$indexF = '';
	}
	
	my %lanes;
	
	my $cwd = getcwd;

	open( LANES, $lanes_fofn ) or die "Cant open lanes fofn file: $!\n";
	while( <LANES> )
	{
		chomp;
		
		my $laneAbsPath = $_;
		if( $_ !~ /^$mRoot.*/ )
		{
			$laneAbsPath = $mRoot.'/'.$_;
		}
		
		chdir($cwd); # mapLane can change our dir and not change us back
		$laneAbsPath = abs_path($laneAbsPath);

		if( ! -d $laneAbsPath )
		{
			print "Cant find lane directory: $laneAbsPath\n";
			next;
		}
		
		if( $laneAbsPath =~ /\/SLX\// )
		{
			if( ! Mapping_SLX_Maq::isLaneMapped( $laneAbsPath ) && ! Mapping_SLX_Maq::mappingInProgress( $laneAbsPath ) )
			{
				print "Mapping Lane: $laneAbsPath\n";
				Mapping_SLX_Maq::mapLane( $laneAbsPath, $lsf_queue, $indexF );
			}
		}
		elsif( $laneAbsPath =~ /\/454\// )
		{
			if( ! Mapping_454_ssaha::isLaneMapped( $laneAbsPath ) && ! Mapping_454_ssaha::mappingInProgress( $laneAbsPath ) )
			{
				print "Mapping Lane: $laneAbsPath\n";
				Mapping_454_ssaha::mapLane( $laneAbsPath, $lsf_queue, $indexF );
			}
		}
	}
	close( LANES );
}

sub makeBam
{
	croak "Usage: rootDir lsf_queue lanes_fofn index_file male_fai female_fai" unless @_ == 6;
	
	my $mRoot = shift;
	my $lsf_queue = shift;
	my $lanes_fofn = shift;
	my $indexF = shift;
	my $male_fai = shift;
	my $female_fai = shift;
	
	croak "Cant find the lanes_fofn file: $lanes_fofn\n" unless -f $lanes_fofn;
	croak "Cant find the male fai file: $male_fai\n" unless -f $male_fai;
	croak "Cant find the female fai file: $female_fai\n" unless -f $female_fai;
	
	if( ! -f $indexF )
	{
		print "Cant find index file: $indexF - ignoring index file\n";
		$indexF = '';
	}
	
	my %lanes;
	
	open( LANES, $lanes_fofn ) or die "Cant open lanes fofn file: $!\n";
	while( <LANES> )
	{
		chomp;
		
		my $laneAbsPath = $_;
		if( $_ !~ /^$mRoot.*/ )
		{
			$laneAbsPath = $mRoot.'/'.$_;
		}
		
		if( ! -d $laneAbsPath )
		{
			print "Cant find lane directory: $laneAbsPath\n";
			next;
		}
		
		if( $laneAbsPath =~ /\/SLX\// )
		{
			if( Mapping_SLX_Maq::isLaneMapped( $laneAbsPath ) )
			{
				print "Making bam for Lane: $laneAbsPath\n";
				Mapping_SLX_Maq::laneToBAM( $laneAbsPath, $lsf_queue, $indexF, $male_fai, $female_fai );
			}
		}
		elsif( $laneAbsPath =~ /\/454\// )
		{
			if( Mapping_454_ssaha::isLaneMapped( $laneAbsPath ) )
			{
				print "Making bam for Lane: $laneAbsPath\n";
				Mapping_454_ssaha::laneToBAM( $laneAbsPath, $lsf_queue, $indexF);
			}
		}
		else
		{
			print "Unknown sequencing technology: $laneAbsPath\n"
		}
	}
	close( LANES );
}

sub checkBam
{
	croak "Usage: rootDir root_dir lanes_fofn" unless @_ == 2;
	
	my $mRoot = shift;
	my $lanes_fofn = shift;
	
	my %lanes;
	
	open( LANES, $lanes_fofn ) or die "Cant open lanes fofn file: $!\n";
	while( <LANES> )
	{
		chomp;
		
		my $laneAbsPath = $_;
		if( $_ !~ /^$mRoot.*/ )
		{
			$laneAbsPath = $mRoot.'/'.$_;
		}
		
		if( ! -d $laneAbsPath )
		{
			print "Cant find lane directory: $laneAbsPath\n";
			next;
		}
		
		if( $laneAbsPath =~ /\/SLX\// )
		{
			if( Mapping_SLX_Maq::isLaneMapped( $laneAbsPath ) )
			{
				Mapping_SLX_Maq::verifyBamFile( $laneAbsPath );
			}
		}
		elsif( $laneAbsPath =~ /\/454\// )
		{
			if( Mapping_454_ssaha::isLaneMapped( $laneAbsPath ) )
			{
				Mapping_454_ssaha::verifyBamFile( $laneAbsPath );
			}
		}
		else
		{
			print "Unknown sequencing technology: $laneAbsPath\n"
		}
	}
	close( LANES );
}

=head2 checkDownloadLocalFile

	Arg [1]    : remote filename to be downloaded
	Arg [2]    : md5 string of remote file
=cut
sub checkDownloadLocalFile
{
	croak "Usage: downloadLocalFile remote_file remote_md5" unless @_ == 2;
	
	my $remoteF = shift;
	my $remote_md5 = shift;
	
	croak "Cant find remote file: $remoteF\n" unless -f $remoteF;
	croak "Cant find remote md5 file: $remote_md5\n" unless -f $remote_md5;
	
	my $remoteFilename = basename( $remoteF );
	if( ! -f $LOCAL_CACHE_DIR.'/'.$remoteFilename || `md5sum $LOCAL_CACHE_DIR/$remoteFilename | awk '{print \$1}'` ne `awk '{print \$1}' $remote_md5` )
	{
		#download
		system( "cp $remoteF $LOCAL_CACHE_DIR" );
		
		my $numRetries = 0;
		while( `md5sum $LOCAL_CACHE_DIR/$remoteFilename | awk '{print \$1}'` ne `awk '{print \$1}' $remote_md5` && $numRetries < 10 )
		{
			print "ERROR: Failed to download local copy of $remoteFilename Trying again....\n";
			system( "cp $remoteF $LOCAL_CACHE_DIR" );
			$numRetries ++;
		}
		
		croak "ERROR: Failed to download local copy of $remoteFilename after 10 retries\n" unless $numRetries < 10;
	}
}

sub cleanPartiallyMappedLanes
{
	croak "Usage: cleanPartiallyMappedLanes rootDir lanes_fofn" unless @_ == 2;
	
	my $mRoot = shift;
	my $lanes_fofn = shift;
	
	croak "Cant find the lanes_fofn file: $lanes_fofn\n" unless -f $lanes_fofn;
	
	my %lanes;
	
	open( LANES, $lanes_fofn ) or die "Cant open lanes fofn file: $!\n";
	while( <LANES> )
	{
		chomp;
		
		my $laneAbsPath = $_;
		if( $_ !~ /^$mRoot.*/ )
		{
			$laneAbsPath = $mRoot.'/'.$_;
		}
		
		if( ! -d $laneAbsPath )
		{
			print "Cant find lane directory: $laneAbsPath\n";
			next;
		}
		
		my ($project, $sample, $platform, $lib, $lane ) = split "/", $_;
		
		if( $platform eq "SLX" )
		{
			if( ! Mapping_SLX_Maq::isLaneMapped( $laneAbsPath ) )
			{
				print "Cleaning Lane: $laneAbsPath\n";
				Mapping_SLX_Maq::removeIntermediateMappingFiles( $laneAbsPath );
			}
		}
		elsif( $platform eq "454" )
		{
			if( ! Mapping_454_ssaha::isLaneMapped( $laneAbsPath ) )
			{
				print "Cleaning Lane: $laneAbsPath\n";
				#Mapping_454_ssaha::removeIntermediateMappingFiles( $laneAbsPath );
			}
		}
	}
	close( LANES );
}

sub calculateClipPoints
{
	croak "Usage: calculateClipPoints rootDir lanes_fofn" unless @_ == 2;
	
	my $mRoot = shift;
	my $lanes_fofn = shift;
	
	croak "Cant find the lanes_fofn file: $lanes_fofn\n" unless -f $lanes_fofn;
	
	my %lanes;
	
	open( LANES, $lanes_fofn ) or die "Cant open lanes fofn file: $!\n";
	while( <LANES> )
	{
		chomp;
		
		my $laneAbsPath = $_;
		if( $_ !~ /^$mRoot.*/ )
		{
			$laneAbsPath = $mRoot.'/'.$_;
		}
		
		if( ! -d $laneAbsPath )
		{
			print "Cant find lane directory: $laneAbsPath\n";
			next;
		}
		
		my ($project, $sample, $platform, $lib, $lane ) = split "/", $_;
		
		if( $platform eq "SLX" )
		{
			Mapping_SLX_Maq::calculateClipPointsLane( $laneAbsPath );
		}
		elsif( $platform eq "454" )
		{
			#Mapping_454_ssaha::calculateClipPointsLane( $laneAbsPath );
		}
	}
	close( LANES );
}

sub mappingHierarchyReport
{
	croak "Usage: mappingHierarchyReport rootDir lanes_fofn output_csv" unless @_ == 3;
	
	my $mRoot = shift;
	my $lanes_fofn = shift;
	my $output_csv = shift;
	
	croak "Cant find the lanes_fofn file: $lanes_fofn\n" unless -f $lanes_fofn;
	
	open( OUT, ">$output_csv" ) or die "Cannot create output file: $output_csv\n";
	print OUT "Status,Center,Project,Individual,Tech,Library,Lane,Fastq0,ReadLength,Fastq1,ReadLength,Fastq2,ReadLength,#Reads,#Bases,#RawReadsMapped,#RawBasesMapped,#RawReadsPaired,#RmdupReadsMapped,#RmdupBasesMapped,ErrorRate\n";
	open( LANES, "$lanes_fofn" ) or die "Cannot open in file: $lanes_fofn\n";
	while( <LANES> )
	{
		chomp;
	
        # This should replace the code below for QC-mapped lines.
        #
        #   my $line = $_;
        #   if ( -s "$line/qc-sample/_stats" ) 
        #   {
        #       open(my $fh,'<',"$line/qc-sample/_stats") or Utils::error("$line/qc-sample/_stats: $!");
        #       my @stats = <$fh>;
        #       close $fh;
        #   
        #       if ( scalar @stats != 1 ) { Utils::error("Expected different content in $line/qc-sample/_stats\n") }
        #       print OUT $stats[0];
        #   
        #       next;
        #   }
		
		my $laneAbsPath = $_;
		if( $_ !~ /^$mRoot.*/ )
		{
			$laneAbsPath = $mRoot.'/'.$_;
		}
		
		if( ! -d $laneAbsPath )
		{
			print "Cant find lane directory: $laneAbsPath\n";
			next;
		}
		
		#get the fastq names
		my @fastq = @{ HierarchyUtilities::getFastqInfo( $laneAbsPath ) };
		
		my $numBases = 0;
		my $numReads = 0;
		for( my $i = 0; $i < @fastq; $i ++ )
		{
			$numReads += $fastq[ $i ][ 2 ];
			$numBases += $fastq[ $i ][ 3 ];
		}
		
		my ($project, $sample, $platform, $lib, $lane ) = (split "/", $_)[-5..-1];
		
		if( $laneAbsPath =~ /\/SLX\// )
		{
			if( Mapping_SLX_Maq::isLaneMapped( $laneAbsPath ) )
			{
				#get the mapping stats
				my $t = Mapping_SLX_Maq::getMappingStats( $laneAbsPath );
				
				if( defined( $t ) )
				{
					print OUT "MAPPED,,$project,$sample,$platform,$lib,$lane,$fastq[ 0 ][ 0 ],$fastq[ 0 ][ 1 ],$fastq[ 1 ][ 0 ],$fastq[ 1 ][ 1 ],$fastq[ 2 ][ 0 ],$fastq[ 2 ][ 1 ],$numReads,$numBases";
					foreach( @$t )
					{
						print OUT ",$_";
					}
					print OUT "\n";
				}
				else
				{
					print OUT "1UNMAPPED,,$project,$sample,$platform,$lib,$lane,$fastq[ 0 ][ 0 ],$fastq[ 0 ][ 1 ],$fastq[ 1 ][ 0 ],$fastq[ 1 ][ 1 ],$fastq[ 2 ][ 0 ],$fastq[ 2 ][ 1 ],$numReads,$numBases\n";
				}
			}
			else
			{
				print OUT "2UNMAPPED,,$project,$sample,$platform,$lib,$lane,$fastq[ 0 ][ 0 ],$fastq[ 0 ][ 1 ],$fastq[ 1 ][ 0 ],$fastq[ 1 ][ 1 ],$fastq[ 2 ][ 0 ],$fastq[ 2 ][ 1 ],$numReads,$numBases\n";
			}
		}
		elsif( $laneAbsPath =~ /\/454\// )
		{
			if( Mapping_454_ssaha::isLaneMapped( $laneAbsPath ) )
			{
				my $t = Mapping_454_ssaha::getMappingStats( $laneAbsPath );
				
				if( defined( $t ) )
				{
					print OUT "MAPPED,,$project,$sample,$platform,$lib,$lane,$fastq[ 0 ][ 0 ],$fastq[ 0 ][ 1 ],$fastq[ 1 ][ 0 ],$fastq[ 1 ][ 1 ],$fastq[ 2 ][ 0 ],$fastq[ 2 ][ 1 ],$numReads,$numBases";
					foreach( @$t )
					{
						print OUT ",$_";
					}
					print OUT "\n";
				}
				else
				{
					print OUT "1UNMAPPED,,$project,$sample,$platform,$lib,$lane,$fastq[ 0 ][ 0 ],$fastq[ 0 ][ 1 ],$fastq[ 1 ][ 0 ],$fastq[ 1 ][ 1 ],$fastq[ 2 ][ 0 ],$fastq[ 2 ][ 1 ],$numReads,$numBases\n";
				}
			}
			else
			{
				print OUT "UNMAPPED,,$project,$sample,$platform,$lib,$lane,$fastq[ 0 ][ 0 ],$fastq[ 0 ][ 1 ],$fastq[ 1 ][ 0 ],$fastq[ 1 ][ 1 ],$fastq[ 2 ][ 0 ],$fastq[ 2 ][ 1 ],$numReads,$numBases\n";
			}
		}
	}
	close( LANES );
	close( OUT );
}


sub unmapped2Bam
{
	croak "Usage: rootDir lsf_queue lanes_fofn index_file" unless @_ == 4;
	
	my $mRoot = shift;
	my $lsf_queue = shift;
	my $lanes_fofn = shift;
	my $indexF = shift;
	
	croak "Cant find the lanes_fofn file: $lanes_fofn\n" unless -f $lanes_fofn;
	
	if( ! -f $indexF )
	{
		print "Cant find index file: $indexF - ignoring index file\n";
		$indexF = '';
	}
	
	my %lanes;
	
	my $cwd = getcwd;

	open( LANES, $lanes_fofn ) or die "Cant open lanes fofn file: $!\n";
	while( <LANES> )
	{
		chomp;
		
		my $laneAbsPath = $_;
		if( $_ !~ /^$mRoot.*/ )
		{
			$laneAbsPath = $mRoot.'/'.$_;
		}
		
		chdir($cwd); # mapLane can change our dir and not change us back
		$laneAbsPath = abs_path($laneAbsPath);
		
		if( ! -d $laneAbsPath )
		{
			print "Cant find lane directory: $laneAbsPath\n";
			next;
		}

		if( $laneAbsPath =~ /\/SLX\// )
		{
			if( ! -f $laneAbsPath."/raw.unmapped.bam" && ! Mapping_SLX_Maq::mappingInProgress( $laneAbsPath ) )
			{
				print "Making unmapped bam for lane: $laneAbsPath\n";
				my $cmd = qq[bsub -q $lsf_queue -o $laneAbsPath/unmapped.o -e $laneAbsPath/unmapped.e "perl -w -e use Mapping_SLX_Maq;Mapping_SLX_Maq::unmapped2Bam( '$laneAbsPath', '$indexF');"];
				system( $cmd );
				#print qq[$cmd\n];
			}
		}
		elsif( $laneAbsPath =~ /\/454\// )
		{
			if( ! Mapping_454_ssaha::isLaneMapped( $laneAbsPath ) && ! Mapping_454_ssaha::mappingInProgress( $laneAbsPath ) )
			{
				#print "Mapping Lane: $laneAbsPath\n";
				#Mapping_454_ssaha::mapLane( $laneAbsPath, $lsf_queue, $indexF );
			}
		}
	}
	close( LANES );
}

1;
